"""Experimental compact decoder using inverse-product error magnitudes.

The state transform and inverse are emitted with the attributed upstream field
arithmetic. A native-field runtime-mask sweep is inserted between them. Common
monomial padding is removed by coefficient rotations, avoiding the fixed-degree
Wronskian transform. No score is claimed before native full-size verification.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from blocked_euclid_sequence import frames
from blocked_monic_euclid import apply as block,schedule
from monic_euclid_local import first_step
from monomial_serial_euclid import apply as serial
from components.trimmed_middle import TrimmedArithmetic as BoundedArithmetic
from compact_monomial_unpad import apply as unpad
from inverse_product_kernel import apply as phase_kernel
from published_dqi import isomorphism
from tower_inverse import map_ops,apply_map,linear_xor
from fast_dqi import FastBuilder

ROOT=Path(__file__).resolve().parents[3]/'4/codex'


def build(N,width,serial_after,M,output,linear_padding=False,virtual_zero=False,virtual_lead=False,low_inverse=False,lazy_flags=False):
    begin=time.time();m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)]
    if low_inverse:
        from components.budget_inverse import bind
        m.inv_word=bind(m)
    alloc=m.Alloc(b)
    flags=list(alloc.open('first_flag',1)) if lazy_flags else b.qubits(N//2)
    E=m.Em(b,*b.qubits(3),alloc);math=BoundedArithmetic(m,E)
    math.virtual_zero=virtual_zero
    stages,final=frames(words,width);bp,bb,_,p=final
    def scale(inv=False):
        with m.inv_word(E,p) as pi,m.opened(E,'scan_scale_product',12) as temp:
            for word in bp:m.mon_slot(E,word,pi if inv else p,p if inv else pi,temp)
    with m.eq_test(E,words[-1],0) as zero:E.cx(zero,flags[0])
    E.cx(flags[0],words[-1][0]);first_step(m,E,math,words)
    for k,count,oldp,old,c,pin in stages:
        if lazy_flags:flags.extend(alloc.open(f'flags_{k}',count))
        if k>=serial_after:
            for j in range(count):oldp,old,c,pin=serial(m,E,math,oldp,old,c,pin,flags[k+j],linear_padding=linear_padding)
        else:block(m,E,math,oldp,old,c,pin,count,flags=flags[k:k+count],monomial=True,suffix=True,linear_padding=linear_padding)
    scale();hook=len(b._items)
    scale(True)
    for k,count,oldp,old,c,pin in reversed(stages):
        if k>=serial_after:
            Q=list(reversed(c))[:2*count];ff,_=schedule(oldp,old,Q)
            for j in reversed(range(count)):
                prev,cur,be,al=ff[j];serial(m,E,math,prev,cur,c[:len(c)-2*j],pin,flags[k+j],True,linear_padding=linear_padding)
        else:block(m,E,math,oldp,old,c,pin,count,True,flags=flags[k:k+count],monomial=True,suffix=True,linear_padding=linear_padding)
        if lazy_flags:alloc.close(f'flags_{k}',flags[k:k+count])
    first_step(m,E,math,words,True);E.cx(flags[0],words[-1][0])
    with m.eq_test(E,words[-1],0) as zero:E.cx(zero,flags[0])
    if lazy_flags:alloc.close('first_flag',flags[:1])
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tape=output.with_suffix('.npy');np.save(tape,rows)
    print(json.dumps(dict(stage='core',N=N,Q=b.next_qubit,ops=len(rows),hook=hook,seconds=time.time()-begin)),flush=True)
    forward,back,transform=isomorphism(m)
    kb=FastBuilder(km.CircuitBuilder());syn=[kb.create_quantum_register(12,f'S{i}') for i in range(N)];masks=[kb.create_classical_register(12,f'E{i}') for i in range(M)]
    allq=[q for word in syn for q in word]+list(kb.alloc_qubits(b.next_qubit-12*N))
    fops=map_ops(forward)
    for word in syn:apply_map(kb,word,fops)
    kb.info_tag('COMPACT_SCAN_START')
    # All core-private scratch is clean at the insertion point. The syndrome
    # words and completion flags remain live and keep their physical numbers.
    flagset=set(flags)
    kb.free(km.array([q for i,q in enumerate(allq) if i>=12*N and i not in flagset]))
    BP=[km.array([allq[q] for q in word]) for word in bp]
    B=[km.array([allq[q] for q in word]) for word in bb]
    P=km.array([allq[q] for q in p]);F=km.array([allq[q] for q in flags])
    # The explicit leading coefficient of scaled Bprev equals pinv. Move
    # that existing word through the coefficient bank, restoring it before
    # the core inverse, rather than keeping a second copy during the scan.
    leadp=P;leadb=None if virtual_lead else kb.alloc_qubits(12)
    if leadb is not None:kb.x(leadb[0]);B.append(leadb)
    BP.append(leadp);counter=kb.alloc_qubits((N//2).bit_length())
    unpad(kb,BP,B,F,counter)
    root=kb.alloc_qubits(12);derivative=kb.alloc_qubits(12);prev=kb.alloc_qubits(12)
    gf=km.GF2Field(12)
    def lin(word,target,c):
        cols=[transform(m.gmul(c,1<<i),back) for i in range(12)]
        linear_xor(kb,word,target,cols)
    for j,mask in enumerate(masks):
        x=transform(gf.pow(2,(-j)%4095),forward)
        def evals():
            for i,word in enumerate(B):
                lin(word,root,m.gpow(x,i))
                if i%2:lin(word,derivative,m.gpow(x,i-1))
            for i,word in enumerate(BP):lin(word,prev,m.gpow(x,i))
            if virtual_lead:
                maximum=N//2
                def constant(value,target):
                    for bit in range(12):
                        if value>>bit&1:kb.x(target[bit])
                def controlled(value,target,control):
                    for bit in range(12):
                        if value>>bit&1:kb.cx(control,target[bit])
                constant(transform(m.gpow(x,maximum),back),root)
                if maximum%2:constant(transform(m.gpow(x,maximum-1),back),derivative)
                for k,flag in enumerate(F):
                    controlled(transform(m.gpow(x,k)^m.gpow(x,k+1),back),root,flag)
                    value=(m.gpow(x,k-1) if k%2 else 0)^(m.gpow(x,k) if (k+1)%2 else 0)
                    controlled(transform(value,back),derivative,flag)
        evals()
        c=gf.pow(gf.pow(2,(-j)%4095),N)
        cols=[gf.mul(c,1<<i) for i in range(12)]
        trans=[sum(((cols[jj]>>i)&1)<<jj for jj in range(12)) for i in range(12)]
        ops=map_ops(trans);apply_map(kb,mask,ops)
        phase_kernel(kb,root,prev,derivative,mask)
        apply_map(kb,mask,ops,True);evals()
        if j%256==0:print(json.dumps(dict(stage='scan',position=j,seconds=time.time()-begin)),flush=True)
    kb.free(root);kb.free(derivative);kb.free(prev)
    unpad(kb,BP,B,F,counter,True);kb.free(counter)
    if leadb is not None:kb.x(leadb[0]);kb.free(leadb)
    kb.info_tag('COMPACT_SCAN_END')
    for word in syn:apply_map(kb,word,fops,True)
    wrapper=output.with_suffix('.wrapper.kmx');circuit=kb.finish_circuit();circuit.to_file(str(wrapper),format='kmx')
    parts=[output.with_suffix('.'+x+'.kmx') for x in ['pre','scan','post']];handles=[p.open('w') for p in parts];part=0
    with wrapper.open() as f:
        for line in f:
            if 'COMPACT_SCAN_START' in line:part=1
            elif 'COMPACT_SCAN_END' in line:part=2
            else:handles[part].write(line)
    for h in handles:h.close()
    empty=output.with_suffix('.marks');empty.write_text('');full=output.with_suffix('.kmx')
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(hook),str(circuit.num_bits),*map(str,parts),str(full),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(full),format='kmx');c.to_file(str(output),format='kmb')
    report=dict(status='BUILT',scope='experimental complete phase decoder; requires native verification',N=N,M=M,width=width,serial_after=serial_after,linear_padding=linear_padding,virtual_zero=virtual_zero,virtual_lead=virtual_lead,low_inverse=low_inverse,lazy_flags=lazy_flags,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-begin)
    Path('results/'+output.stem+'_build.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    return c


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=14);p.add_argument('--block',type=int,default=2);p.add_argument('--serial-after',type=int,default=5);p.add_argument('--m',type=int,default=63);p.add_argument('--output',default='candidates/compact_product_n14_m63.kmb');p.add_argument('--linear-padding',action='store_true');p.add_argument('--virtual-zero',action='store_true');p.add_argument('--virtual-lead',action='store_true');p.add_argument('--low-inverse',action='store_true');p.add_argument('--lazy-flags',action='store_true');a=p.parse_args();build(a.n,a.block,a.serial_after,a.m,Path(a.output),a.linear_padding,a.virtual_zero,a.virtual_lead,a.low_inverse,a.lazy_flags)
