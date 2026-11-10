"""46T inverse-only phase body and in-place Wronskian preparation.

The Wronskian identity and the proposed 46T cubic adaptation were already
described in round-4 WRONSKIAN_RESEARCH.md. This file implements the body;
sylvester_circuit.py supplies a different Euclidean factorization.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from tower_inverse import map_ops,apply_map,linear_xor,init_inverse
from cubic_inverse_product_kernel import ClassicalBuilder,classical_product_xor,trace_inverse
from cubic_trace_field import data
from import_published import load
from low_workspace_inverse import bind
from sylvester_circuit import apply as sylvester,output_layout,generic_chart
from algebra import Field,trim

HERE=Path(__file__).resolve().parent


def phase(b,root,h,mask):
    active=b.alloc_bits(1)[0];b.bit_store(active,True,control=mask);b.push_condition(active)
    ops=map_ops(data()['dual_columns']);apply_map(b,mask,ops)
    ik=b.alloc_bits(12);sc=init_inverse(ClassicalBuilder(b),mask,ik)
    scaled=b.alloc_qubits(12);classical_product_xor(b,h,ik,scaled)
    trace_inverse(b,root,scaled)
    classical_product_xor(b,h,ik,scaled);b.free(scaled)
    b.bit_store(ik,False);b.bit_store(sc,False);b.free(ik);b.free(sc)
    apply_map(b,mask,ops,True);b.pop_condition();b.bit_store(active,False);b.free(km.array([active]))


def binary_map(m,E,word,exponent,inverse=False):
    ops=map_ops([m.gpow(1<<j,exponent) for j in range(12)])
    for op,i,j in (reversed(ops) if inverse else ops):
        if op=='cx':E.cx(word[i],word[j])
        else:E.swaps(None,[(word[i],word[j])])


def prepare(m,E,A,B,inverse=False):
    """Monic B, deg A < deg B. Replace A by sqrt((A B)').

    B is restored. Result order is a returned wire list, with no physical
    output permutation needed. Generic E/O Euclidean charts only.
    """
    assert len(B)==len(A)+1 and len(A)>=1
    if len(A)==1:
        binary_map(m,E,A[0],2048,inverse)
        return A
    t=len(A)//2
    C,D=A[::2],A[1::2];even,odd=B[::2],B[1::2]
    if len(A)%2:
        c,d,extra,offset,other,monic=D,C[:-1],C[-1],t,odd,True
    else:
        c,d,extra,offset,other,monic=D[:-1],C,D[-1],t-1,even,False
    H=output_layout(len(even)-1,len(odd)-1,c,d)+[extra]
    def cross():
        for j in range(t):E.xor(m.prod_rt(H[offset+j],extra,other[j]))
    if inverse:
        for w in B:binary_map(m,E,w,2048)
        cross();sylvester(m,E,even,odd,c,d,True,g_one=monic)
        for w in A:binary_map(m,E,w,2048,True)
        for w in B:binary_map(m,E,w,2048,True)
    else:
        for w in A+B:binary_map(m,E,w,2048)
        sylvester(m,E,even,odd,c,d,g_one=monic);cross()
        for w in B:binary_map(m,E,w,2048,True)
    return H


def point_test():
    b=FastBuilder(km.CircuitBuilder());root=b.create_quantum_register(12,'root');h=b.create_quantum_register(12,'h');mask=b.create_classical_register(12,'mask')
    phase(b,root,h,mask);c=b.finish_circuit();c.to_file(str(HERE/'candidates/inverse_only_point.kmb'),format='kmb')
    F=Field();values=np.arange(4096);inv=np.array([F.inv(int(v)) if v else 0 for v in values]);rng=np.random.default_rng(620193)
    bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    public=set(root)|set(h);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    cases=[(np.zeros(4096,dtype=int),np.full(4096,1<<j)) for j in range(12)]
    cases += [(np.zeros(4096,dtype=int),values),(rng.integers(1,4096,4096),values),(np.zeros(4096,dtype=int),np.zeros(4096,dtype=int))]
    for i,(rv,mv) in enumerate(cases):
        sim=km.Simulator(batch_size=4096,seed=51026+i);sim.use_same_registers_as(c)
        for reg,v in [(root,rv),(h,values),(mask,mv)]:sim.write_across_shots(reg,bits(v))
        sim.do(c)
        want=np.array([r==0 and (int(v)&int(mm)).bit_count()%2 for r,v,mm in zip(rv,inv,mv)],dtype=bool)
        assert np.array_equal(sim.read_phase_flipped_across_shots(),want)
        for reg,v in [(root,rv),(h,values),(mask,mv)]:assert np.array_equal(sim.read_across_shots(reg),bits(v))
        assert not sim.read_across_shots(scratch).any()
    r=dict(status='PASS',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=4096*len(cases),
           scope='all divisors against every mask basis, all runtime masks, nonroots, zero masks; input/mask/scratch/phase',
           sha256=hashlib.sha256((HERE/'candidates/inverse_only_point.kmb').read_bytes()).hexdigest())
    assert r['T']==46
    (HERE/'results/inverse_only_point.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)


def prep_test(degree,shots,salt,profile=False,points=1):
    m=load();m.inv_word=bind(m);b=m.Builder();A=[b.qubits(12) for _ in range(degree)];B=[b.qubits(12) for _ in range(degree)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));lead=E.A.open('h_monic_lead',12);E.x(lead[0]);BB=B+[lead]
    H=prepare(m,E,A,BB);forward_T=b._cur_stats['toffoli'];hook=len(b._items)
    prepare(m,E,A,BB,True);E.x(lead[0]);E.A.close('h_monic_lead',lead)
    r=dict(degree=degree,Q=b.peak_qubits,T=b._cur_stats['toffoli'],forward_T=forward_T,
           public_Q=24*degree,scratch_Q=b.peak_qubits-24*degree,
           scope='full H preparation and erasure, monic fixed degree; generic Euclidean chart')
    if profile:
        (HERE/f'results/h_degree{degree}_profile.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);return
    tag=f'h_degree{degree}'+(f'_m{points}' if points>1 else '')+'_phase'
    stem=HERE/'candidates'/tag
    tape=stem.with_suffix('.npy')
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
    kb=FastBuilder(km.CircuitBuilder())
    for key,words in [('A',A),('B',B)]:
        for i in range(degree):kb.create_quantum_register(12,f'{key}{i}')
    masks=[kb.create_classical_register(12,'M' if points==1 else f'M{j}') for j in range(points)];q=list(km.q(i) for i in range(b.next_qubit))
    # The private lead is LIVE and equals one at the hook; reserve it.
    kb.alloc_qubits(b.next_qubit-24*degree)
    live=set(lead)
    kb.free(km.array([q[i] for i in range(24*degree,b.next_qubit) if i not in live]))
    root=kb.alloc_qubits(12);value=kb.alloc_qubits(12)
    from published_dqi import isomorphism
    fw,bk,tr=isomorphism(m)
    host=km.GF2Field(12)
    xs=[]
    for j,mask in enumerate(masks):
        beta=host.pow(2,(-j)%4095);x=tr(beta,fw);xs.append(x)
        def evaluate():
            lead_value=tr(m.gpow(x,degree),bk)
            for bit in range(12):
                if lead_value>>bit&1:kb.x(root[bit])
            for i,w in enumerate(B):linear_xor(kb,km.array([q[z] for z in w]),root,[tr(m.gmul(m.gpow(x,i),1<<bit),bk) for bit in range(12)])
            for i,w in enumerate(H):linear_xor(kb,km.array([q[z] for z in w]),value,[tr(m.gmul(m.gpow(x,i),1<<bit),bk) for bit in range(12)])
        # Absorb e=beta^70/H^2 into this shot's classical coordinate mask.
        constant=host.pow(beta,70)
        columns=[sum(((host.mul(constant,host.mul(1<<z,1<<z))>>i)&1)<<z for z in range(12)) for i in range(12)]
        maskops=map_ops(columns)
        evaluate();apply_map(kb,mask,maskops);phase(kb,root,value,mask);apply_map(kb,mask,maskops,True);evaluate()
    kb.free(root);kb.free(value)
    scan=stem.with_suffix('.scan.kmx');kb.finish_circuit().to_file(str(scan),format='kmx')
    # Split register declarations (pre) from the phase body; retain allocations
    # as actual wire references, which the native parser counts globally.
    lines=scan.read_text().splitlines(True);pre=[];body=[]
    for line in lines:
        (pre if line.startswith(('REGISTER','APPEND_TO_REGISTER')) else body).append(line)
    prepath=stem.with_suffix('.pre.kmx');prepath.write_text(''.join(pre));scan.write_text(''.join(body))
    empty=HERE/'parameters/empty';empty.write_text('');out=stem.with_suffix('.kmx')
    subprocess.run([str(HERE/'../../4/codex/translate_prefix'),str(tape),str(hook),str(12*points),str(prepath),str(scan),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    circ=km.Circuit.from_file(str(out),format='kmx');circ.to_file(str(stem.with_suffix('.kmb')),format='kmb')
    F=Field(12,0x10eb);seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:8],'little');rng=np.random.default_rng(seed)
    arrays=[];expect=[]
    while len(arrays)<shots:
        # Every second locator has one of the scanned roots; all other roots
        # are outside this small scan, so all point tests are exercised.
        available=np.array([v for v in range(1,4096) if v not in xs])
        roots=rng.choice(available,degree-(len(arrays)%2==0),replace=False).tolist()
        if len(arrays)%2==0:roots.append(xs[(len(arrays)//2)%points])
        poly=[1]
        for x in roots:poly=F.product(poly,[x,1])
        aa=rng.integers(0,4096,degree).tolist()
        if not generic_chart(F,[F.sqrt(x) for x in poly[::2]],[F.sqrt(x) for x in poly[1::2]]):continue
        mm=rng.integers(0,4096,points).tolist();parity=0
        for x,maskvalue in zip(xs,mm):
            wv=F.mul(F.evaluate(aa,x),F.evaluate(F.derivative(poly),x))
            ev=F.mul(F.pow(x,70),F.inv(wv)) if wv else 0
            ev=tr(ev,bk)
            if F.evaluate(poly,x)==0:parity^=(ev&maskvalue).bit_count()%2
        arrays.append((aa,poly[:-1],mm));expect.append(parity)
    bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(circ)
    for key,j in [('A',0),('B',1)]:
        for i in range(degree):sim.write_across_shots(circ.register_data[f'{key}{i}'],bits([a[j][i] for a in arrays]))
    for j in range(points):sim.write_across_shots(circ.register_data['M' if points==1 else f'M{j}'],bits([a[2][j] for a in arrays]))
    sim.do(circ)
    assert np.array_equal(sim.read_phase_flipped_across_shots(),np.array(expect,dtype=bool)), 'H phase'
    for key,j in [('A',0),('B',1)]:
        for i in range(degree):assert np.array_equal(sim.read_across_shots(circ.register_data[f'{key}{i}']),bits([a[j][i] for a in arrays]))
    public={q for reg in circ.register_data.values() for q in reg if isinstance(q,km.q)};scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any()
    for j in range(points):assert np.array_equal(sim.read_across_shots(circ.register_data['M' if points==1 else f'M{j}']),bits([a[2][j] for a in arrays]))
    sim.do(circ);assert not sim.read_phase_flipped_across_shots().any()
    assert not sim.read_across_shots(scratch).any()
    r.update(status='PASS',shots=shots,salt=salt,points=points,phase_sandwich_Q=circ.num_qubits,phase_sandwich_T=circ.operation_statistics().num_toffoli_like,
             sha256=hashlib.sha256(stem.with_suffix('.kmb').read_bytes()).hexdigest())
    (HERE/f'results/h_degree{degree}{"_m"+str(points) if points>1 else ""}_{salt}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--point',action='store_true');p.add_argument('--degree',type=int,default=7);p.add_argument('--shots',type=int,default=512);p.add_argument('--salt',default='core-phase-H-a');p.add_argument('--profile',action='store_true');p.add_argument('--m',type=int,default=1);a=p.parse_args()
    point_test() if a.point else prep_test(a.degree,a.shots,a.salt,a.profile,a.m)
