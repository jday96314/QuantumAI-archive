"""Exhaustive integer and native component certificates, both top architectures."""
import hashlib,json,math,time
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt.verify_submission_invsqrt import half_turns_to_units
from exponent_kernels import borrowed_phase,walsh_phase,walsh_coefficients,circular_phase
HERE=Path(__file__).resolve().parent
N=3133452
TABLE=np.load(HERE/'candidates/margin_curve8_t5.9e-08.npy')


def product(a,z,e):
    out=np.zeros_like(z)
    bit=lambda j:0 if j<0 else (z>>min(j,13))&1
    for g in range(7):
        digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
        for mag in (1,2):
            ee=e+2*g+(mag==2)
            v=a<<ee if ee>=0 else (a>>-ee)+((a>>(-ee-1))&1)
            out+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
    return out


def evaluate(xs):
    k=np.zeros_like(xs)
    for j in range(21):k+=xs<(1<<(21-j))
    m=xs<<k;s=k//2
    c0,c1,c2=TABLE[((m>>13)&255)+((k&1)<<8)].T
    z=(m&16383)^8192
    slope=c1-product(c2,z,-11);p=c0-product(slope,z,-12)
    return k,m,s,p


def integer_check(coefficients):
    start=time.monotonic();target_file=HERE.parent/'inv-phase/reference/exact_target.npz'
    # Prior artifacts are read-only snapshots. Reference computed by host integer
    # formula there, independently spot-checked here including worst/margin x.
    with np.load(target_file) as f:target=f['target'];tol=f['tol']
    for x in [1,2,3,4,958102,2965859,N]:
        r=math.isqrt((1<<130)//x)
        assert int(target[x-1])==((r+1)//2)%(1<<64) and int(tol[x-1])==r>>25
    maxratio=0.;margin=1<<64;worst=0;amin=0
    for first in range(1,N+1,131072):
        xs=np.arange(first,min(N+1,first+131072),dtype=np.int64)
        k,m,s,p=evaluate(xs)
        d=np.zeros_like(xs)
        for j in range(10):d|=((m>>j)&1)<<(9-j)
        word=(p|(d<<27)).astype(np.uint64)
        rot=((word<<s.astype(np.uint64))|(word>>(37-s).astype(np.uint64)))&np.uint64((1<<37)-1)
        mask=(np.uint64(1)<<(27+s).astype(np.uint64))-np.uint64(1)
        borrowed=rot&mask
        assert np.array_equal(borrowed,(p<<s).astype(np.uint64))
        invrot=((rot>>s.astype(np.uint64))|(rot<<(37-s).astype(np.uint64)))&np.uint64((1<<37)-1)
        assert np.array_equal(invrot,word)
        # Independent integer identity for the refined circular phase kernel.
        pc=p.astype(np.uint64);rc=((pc<<s.astype(np.uint64))|(pc>>(27-s).astype(np.uint64)))&np.uint64((1<<27)-1)
        cmask=(np.uint64(1)<<s.astype(np.uint64))-np.uint64(1)
        circular=rc+((rc&cmask)*np.uint64((1<<27)-1))
        assert np.array_equal(circular,(p<<s).astype(np.uint64))
        packed=np.where(k==0,m,m-(1<<21)+(np.int64(1)<<np.maximum(0,k-1)))
        decoded=np.zeros_like(k)
        low=packed&-packed
        for j in range(21):decoded+=((packed<(1<<21))&(low==(1<<j)))*(j+1)
        assert np.array_equal(decoded,k)
        restored=np.where(decoded==0,packed,packed-(np.int64(1)<<np.maximum(0,decoded-1))+(1<<21))
        assert np.array_equal(restored,m)
        actual=borrowed<<np.uint64(27)
        tt=target[first-1:first-1+len(xs)];to=tol[first-1:first-1+len(xs)]
        used=np.abs((actual-tt).view(np.int64)).astype(np.uint64)
        ratios=used/to;j=int(ratios.argmax())
        if ratios[j]>maxratio:maxratio=float(ratios[j]);worst=int(xs[j])
        mm=to.astype(np.int64)-used.astype(np.int64);j=int(mm.argmin())
        if int(mm[j])<margin:margin=int(mm[j]);amin=int(xs[j])
    for s in range(11):assert sum(v*(-1)**((s&m).bit_count()%2) for m,v in coefficients.items())==1<<s
    return dict(status='PASS' if margin>0 else 'FAIL',inputs=N,borrowed_phase_identity='exact',
        inverse_rotation='exact',circular_phase_identity='exact',packed_exponent_identity='exact',
        min_margin_units=margin,min_margin_turns=margin/2**64,margin_x=amin,
        worst_ratio=maxratio,worst_x=worst,seconds=time.monotonic()-start)


def component(mode,ns=27,nx=22,head=10,denom=37,coefficients=None):
    b=km.CircuitBuilder();x=b.create_quantum_register(nx,'x');p=b.create_quantum_register(ns,'p')
    flags=b.create_quantum_register(5,'flags')
    if mode=='borrowed':borrowed_phase(b,p,x,flags,head_bits=head,denominator=denom)
    elif mode=='circular':circular_phase(b,p,flags,denominator=denom)
    else:walsh_phase(b,p,flags[1:],coefficients,denominator=denom)
    return b.finish_circuit(),x,p,flags


def native_check(mode,coefficients):
    # Exhaustive small-word semantics for every permitted input, output word,
    # and exponent. Extra full-width samples cover actual polynomial values.
    c,x,p,f=component(mode,ns=4,nx=6,head=3,denom=12,coefficients=coefficients)
    small=[(xx,pp,ss) for ss in range(4) for xx in range(64) if xx% (1<<(2*ss))==0 for pp in range(16)]
    cases=np.array(small,dtype=np.int64)
    seeds=['census-components-1','census-components-2']
    def run(c,x,p,f,cases,denom,salt):
        batch=256;sim=km.Simulator(batch_size=batch,seed=int(hashlib.sha256(salt.encode()).hexdigest(),16),counting_operations=True)
        sim.use_same_registers_as(c);public=set(x)|set(p)|set(f)
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool)
        for first in range(0,len(cases),batch):
            values=cases[first:first+batch];n=len(values)
            if n<batch:values=np.resize(values,(batch,3))
            xx,pp,ss=values.T;sim.clear_for_shot()
            for reg,v in ((x,xx),(p,pp),(f,ss*2)):sim.write_across_shots(reg,bits(v,len(reg)))
            sim.do(c)
            for reg,v in ((x,xx),(p,pp),(f,ss*2)):assert np.array_equal(sim.read_across_shots(reg),bits(v,len(reg))),(mode,first,'register')
            assert not sim.read_across_shots(scratch).any(),(mode,first,'scratch')
            actual=np.array([half_turns_to_units(sim.read_shot_phase(j)) for j in range(batch)],dtype=np.uint64)
            expected=(pp.astype(np.uint64)<<ss.astype(np.uint64))<<np.uint64(64-denom)
            assert np.array_equal(actual,expected),(mode,first,'phase')
        return math.ceil(len(cases)/batch)*batch
    shots=sum(run(c,x,p,f,cases,12,s) for s in seeds)
    full,x,p,f=component(mode,coefficients=coefficients)
    rng=np.random.default_rng(20261024)
    xs=np.unique(np.r_[np.arange(1,4097),rng.integers(1,N+1,8192),[958102,2965859,N]])
    k,m,s,pp=evaluate(xs);fullcases=np.column_stack((m,pp,s))
    shotsfull=run(full,x,p,f,fullcases,37,mode+'-full-width-1')
    full.to_file(str(HERE/f'candidates/{mode}_component.kmb'),format='kmb')
    st=full.operation_statistics()
    return dict(status='NATIVE_PASS',scope='component, not accepted inverse',exhaustive_small_cases=len(cases),
        small_shots=shots,full_width_cases=len(xs),full_width_shots=shotsfull,
        component_Q=full.num_qubits,Toffolis=st.num_toffoli_like,
        rotations=st.num_z_pow(min_clifford_hierarchy_level=4),direct_T=st.num_t_like)


def main():
    co,search=walsh_coefficients()
    report=dict(walsh_search=search,walsh_coefficients={str(k):str(v) for k,v in co.items()})
    (HERE/'results/walsh_search.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    report['integer']=integer_check(co)
    (HERE/'results/prototypes.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['integer']),flush=True)
    for mode in ('borrowed','walsh','circular'):
        report[mode]=native_check(mode,co)
        (HERE/'results/prototypes.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report[mode]),flush=True)

if __name__=='__main__':main()
