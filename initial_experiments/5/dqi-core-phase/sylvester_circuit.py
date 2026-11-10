"""In-place Euclidean Sylvester transform, with literal inverse.

Uses read-only Apache-2.0 Gu field MACs (via round-3 import_published)
and the round-4 low-workspace inverse. Generic Euclidean degree chart;
zero pivots are armed, so restoration is unconditional but the Sylvester
formula on such paths is deliberately NOT asserted.
"""
from contextlib import contextmanager
from pathlib import Path
import argparse, hashlib, json, subprocess
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from algebra import Field

HERE = Path(__file__).resolve().parent
TRANSLATOR = HERE / '../../4/codex/translate_prefix'


def output_layout(m, n, c, d):
    if not n: return list(d)
    k = n-1
    return output_layout(n, k, d[:k], c) + list(d[k:])


@contextmanager
def factor(m, E, source, one=False):
    if one:
        yield None
    else:
        with m.inv_word(E, source) as si: yield si


def scale(m, E, words, source, inverse=False, one=False):
    if one or not words: return
    with factor(m,E,source) as si, m.opened(E,'sylv_scale_product',12) as t:
        for w in words: m.mon_slot(E,w,si if inverse else source,source if inverse else si,t)


def divide(m,E,f,g,inverse=False,one=False):
    """(f,g) -> (remainder || quotient,g), in the original f slots."""
    n = len(g)-1; mm = len(f)-1
    with factor(m,E,g[-1],one) as gi, m.opened(E,'sylv_div_product',12) as t:
        for i in (range(n,mm+1) if inverse else range(mm,n-1,-1)):
            if not inverse and not one: m.mon_slot(E,f[i],gi,g[-1],t)
            for j in range(n): E.xor(m.prod_rt(f[i-n+j],f[i],g[j]))
            if inverse and not one: m.mon_slot(E,f[i],g[-1],gi,t)


def triangular(m,E,g,high,inverse=False,one=False):
    n=len(g)-1; L=len(high)
    with factor(m,E,g[-1],one) as gi, m.opened(E,'sylv_tri_product',12) as t:
        for j in (reversed(range(L)) if inverse else range(L)):
            if not inverse and not one: m.mon_slot(E,high[j],g[-1],gi,t)
            for ell in range(j+1,L):
                idx=n+j-ell
                if idx>=0: E.xor(m.prod_rt(high[j],g[idx],high[ell]))
            if inverse and not one: m.mon_slot(E,high[j],gi,g[-1],t)


def apply(m,E,f,g,c,d,inverse=False,g_one=False):
    """In place c*f+d*g, f deg m>=g deg n, c n words, d m words.

    f and g are restored at EACH recursion return. The output coefficient
    order is output_layout(m,n,c,d). Public and scratch restoration by the
    complete forward/inverse sandwich holds even on armed exceptional paths.
    """
    mm,n=len(f)-1,len(g)-1
    assert mm>=n>=0 and len(c)==n and len(d)==mm
    # Retain ONE zero bit, not eq_test's whole conjunction ladder.
    # Recompute the predicate after restoring g to clear that retained bit.
    with m.opened(E,'sylv_zero_flag',1) as zw:
        zero=zw[0]
        with m.eq_test(E,g[-1],0) as z:E.cx(z,zero)
        E.cx(zero,g[-1][0])
        if not n:
            scale(m,E,d,g[0],inverse,g_one)
        else:
            k=n-1; q=f[n:]; r=f[:n]
            low=output_layout(n,k,d[:k],c); high=d[k:]
            def shear(): E.xor(m.polyprod_rt(q,c,[[w] for w in d]))
            def cross():
                targets=[[low[i+k]] if i<n else [] for i in range(len(g)+len(high)-1)]
                E.xor(m.polyprod_rt(g,high,targets))
            divide(m,E,f,g,one=g_one)
            if inverse:
                triangular(m,E,g,high,True,g_one);cross()
                apply(m,E,g,r,d[:k],c,True)
                shear()
            else:
                shear();apply(m,E,g,r,d[:k],c)
                cross();triangular(m,E,g,high,False,g_one)
            divide(m,E,f,g,True,g_one)
        E.cx(zero,g[-1][0])
        with m.eq_test(E,g[-1],0) as z:E.cx(z,zero)
    return output_layout(mm,n,c,d)


def emitter(t):
    m=load();m.inv_word=bind(m);b=m.Builder()
    # f has degree t, g is monic degree t, and c,d each have t words.
    f=[b.qubits(12) for _ in range(t+1)]
    g=[b.qubits(12) for _ in range(t)]
    c=[b.qubits(12) for _ in range(t)]
    d=[b.qubits(12) for _ in range(t)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    lead=E.A.open('sylv_monic_lead',12);E.x(lead[0]);gfull=g+[lead]
    return m,b,E,f,g,c,d,gfull


def build(t,which,profile_only=False):
    m,b,E,f,g,c,d,gfull=emitter(t)
    before=b._cur_stats['toffoli']
    directions={'both':[False,True],'inverse_both':[True,False],'forward':[False],'inverse':[True]}
    for inverse in directions[which]:
        apply(m,E,f,gfull,c,d,inverse,g_one=True)
    E.x(gfull[-1][0]);E.A.close('sylv_monic_lead',gfull[-1])
    stats=dict(t=t,which=which,Q=b.peak_qubits,T=b._cur_stats['toffoli']-before,
               public_Q=12*(4*t+1),scratch_Q=b.peak_qubits-12*(4*t+1),
               final_live=b.active_qubits,layout=output_layout(t,t,c,d))
    if profile_only:return stats
    stem=HERE/f'candidates/sylvester_t{t}_{which}'
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tape=HERE/f'parameters/sylvester_t{t}_{which}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for key,words in [('f',f),('g',g),('c',c),('d',d)]:
        for i in range(len(words)):kb.create_quantum_register(12,f'{key}{i}')
    pre=stem.with_suffix('.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=stem.with_suffix('.kmx')
    subprocess.run([str(TRANSLATOR),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    circuit=km.Circuit.from_file(str(out),format='kmx');circuit.to_file(str(stem.with_suffix('.kmb')),format='kmb')
    assert circuit.num_qubits==stats['Q']
    assert circuit.operation_statistics().num_toffoli_like==stats['T']
    stats['sha256']=hashlib.sha256(stem.with_suffix('.kmb').read_bytes()).hexdigest()
    stats['layout']=[[int(x) for x in w] for w in stats['layout']]
    return circuit,stats


def generic_chart(F,f,g):
    if not g[-1]:return False
    while len(g)>1:
        _,r=F.divide(f,g)
        if len(r)!=len(g)-1:return False
        f,g=g,r
    return bool(g[0])


def test(t,shots,salt):
    F=Field(12,0x10eb);seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:8],'little')
    rng=np.random.default_rng(seed);data=[];accepted=0
    while accepted<shots:
        f=rng.integers(0,4096,t+1).tolist();g=rng.integers(0,4096,t).tolist()+[1]
        if not generic_chart(F,f,g):continue
        c=rng.integers(0,4096,t).tolist();d=rng.integers(0,4096,t).tolist()
        want=F.add(F.product(c,f),F.product(d,g));want += [0]*(2*t-len(want))
        data.append((f,g[:-1],c,d,want));accepted+=1
    values={key:np.array([case[j] for case in data],dtype=np.uint16) for j,key in enumerate(['f','g','c','d','out'])}
    bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    reports=[]
    for which in ['forward','inverse','both']:
        circ,stats=build(t,which);sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(circ)
        for key in ['f','g','c','d']:
            for i in range(values[key].shape[1]):sim.write_across_shots(circ.register_data[f'{key}{i}'],bits(values[key][:,i]))
        layout=stats['layout']
        if which=='inverse':
            for i,word in enumerate(layout):sim.write_across_shots(km.array([km.q(q) for q in word]),bits(values['out'][:,i]))
        sim.do(circ)
        for key in ['f','g']:
            for i in range(values[key].shape[1]):assert np.array_equal(sim.read_across_shots(circ.register_data[f'{key}{i}']),bits(values[key][:,i])),(t,which,key,i)
        if which=='forward':
            for i,word in enumerate(layout):assert np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in word])),bits(values['out'][:,i])),(t,which,'out',i)
        else:
            for key in ['c','d']:
                for i in range(t):assert np.array_equal(sim.read_across_shots(circ.register_data[f'{key}{i}']),bits(values[key][:,i])),(t,which,key,i)
        public={q for r in circ.register_data.values() for q in r};scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(t,which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(t,which,'phase')
        stats.update(status='PASS',shots=shots,salt=salt);reports.append(stats)
        print(json.dumps({k:v for k,v in stats.items() if k!='layout'}),flush=True)
    (HERE/f'results/sylvester_t{t}_{salt}.json').write_text(json.dumps(reports,indent=2)+'\n')


def restoration_test(t,shots,salt):
    seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:8],'little');rng=np.random.default_rng(seed)
    values={key:rng.integers(0,4096,(shots,t+(key=='f')),dtype=np.uint16) for key in ['f','g','c','d']}
    values['f'][:shots//4]=0;values['g'][shots//4:shots//2]=0
    values['f'][shots//2:3*shots//4]=0;values['g'][shots//2:3*shots//4]=0
    bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    reports=[]
    for which in ['both','inverse_both']:
        circ,stats=build(t,which);sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(circ)
        for key,arr in values.items():
            for i in range(arr.shape[1]):sim.write_across_shots(circ.register_data[f'{key}{i}'],bits(arr[:,i]))
        sim.do(circ)
        for key,arr in values.items():
            for i in range(arr.shape[1]):assert np.array_equal(sim.read_across_shots(circ.register_data[f'{key}{i}']),bits(arr[:,i])),(which,key,i)
        public={q for reg in circ.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any()
        assert not sim.read_phase_flipped_across_shots().any()
        stats.pop('layout');stats.update(status='PASS',shots=shots,salt=salt,scope='arbitrary operands, forced zero remainder chains, both composition orders; restoration only on exceptional charts');reports.append(stats)
        print(json.dumps(stats),flush=True)
    (HERE/f'results/sylvester_t{t}_restoration_{salt}.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--t',type=int,default=3);p.add_argument('--shots',type=int,default=512);p.add_argument('--salt',default='core-phase-sylvester-a');p.add_argument('--profile',action='store_true');p.add_argument('--restoration',action='store_true');a=p.parse_args()
    if a.profile:
        r=build(a.t,'both',True);r.pop('layout');(HERE/f'results/sylvester_t{a.t}_profile.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
    elif a.restoration:restoration_test(a.t,a.shots,a.salt)
    else:test(a.t,a.shots,a.salt)
