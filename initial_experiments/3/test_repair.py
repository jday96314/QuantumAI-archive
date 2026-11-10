"""Quantum simulation checks including deliberately vanishing initial pivots."""
import sys,time,os,importlib
import numpy as np
import kickmix as km
from dqi.verify_submission_dqi import ProblemSpecDQI,FieldTables,sample_shots,verify
build_circuit=importlib.import_module(os.environ.get('DQI_MODULE','repair_dqi')).build_circuit

gf=km.GF2Field(12)
tables=FieldTables(gf)
rng=np.random.default_rng(88)


def forced_shots(m,n,N):
    syn=np.zeros((N,n),np.int64)
    masks=rng.integers(0,4096,(N,m),dtype=np.int64)
    phases=np.zeros(N,bool)
    for s in range(N):
        w=int(rng.integers(2,n//2+1))
        pos=rng.choice(m,w,replace=False)
        val=rng.integers(1,4096,w,dtype=np.int64)
        roots=[gf.pow(gf.primitive_element,int(p)) for p in pos]
        top=0
        for e,r in zip(val[:-1],roots[:-1]):top ^= gf.mul(int(e),gf.pow(r,n))
        val[-1]=gf.mul(top,gf.pow(gf.pow(roots[-1],n),4094))
        for j in range(n):
            for e,r in zip(val,roots):syn[s,j] ^= gf.mul(int(e),gf.pow(r,j+1))
        phases[s]=sum(int(e).bit_count() for e in val&masks[s,pos])%2==1
    return syn,masks,phases


for n in ([int(sys.argv[1])] if len(sys.argv)>1 else [6,10,20]):
    m=max(31,n+1)
    start=time.time();c=build_circuit(m=m,n=n)
    print('built',n,c.num_qubits,c.operation_statistics().num_toffoli_like,time.time()-start,flush=True)
    spec=ProblemSpecDQI.benchmark_instance(m=m,n=n,num_shots=1024)
    for forced in [False,True]:
        shots=forced_shots(m,n,1024) if forced else sample_shots(spec,tables,rng,1024)
        result=verify(c,spec,shots=shots,max_workers=4,log2_delta=-16)
        print('test',n,'forced',forced,{k:result[k] for k in ('value_errors','phase_errors','toffoli_count')},flush=True)
        assert result['value_errors']==0 and result['phase_errors']<=2
