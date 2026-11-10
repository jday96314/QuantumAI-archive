"""Diagnose the borrowed Psi component outside its squarefree-locator promise.

No exception repair is attempted. Failures here disqualify any claim that
this experimental half-degree pass is a total replacement on all core paths.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from published_dqi import isomorphism
from oracle import product


def main():
    m=load();fw,bk,transform=isomorphism(m)
    layout=json.loads(Path('parameters/psi_layout.json').read_text())
    c=km.Circuit.from_file('candidates/psi_roundtrip.kmb',format='kmb')
    rng=np.random.default_rng(536193);N=256;labels=[];sigma=[];shifts=[]
    for s in range(N):
        typ=s%8
        if typ==0:label,p='zero polynomial',[]
        elif typ==1:label,p='constant',[1]
        elif typ==2:label,p='repeated zero root',[0]*35+[1]
        elif typ==3:label,p='zero derivative',[1]+[0]*33+[1]
        elif typ==4:label,p='repeated nonzero root',product(product([1,2],[1,2]),[1,3])
        elif typ==5:label,p='degree-35 arbitrary monic',rng.integers(0,4096,35).tolist()+[1]
        elif typ==6:label,p='degree-34 arbitrary monic',rng.integers(0,4096,34).tolist()+[1]
        else:label,p='squarefree zero root',[0,1]
        labels.append(label);sigma.append(p+[0]*(36-len(p)));shifts.append((35-max(0,len(p)-1))//2)
    sigma=np.array(sigma).T;prior=rng.integers(0,4096,(35,N));flags=rng.integers(0,2,(35,N)).astype(bool);shifts=np.array(shifts)
    bits=lambda a,n:((np.asarray(a)[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=N,seed=194519);sim.use_same_registers_as(c)
    for prefix,arr in [('sigma',sigma),('prior',prior)]:
        for j,word in enumerate(arr):sim.write_across_shots(c.register_data[f'{prefix}{j}'],bits([transform(int(v),fw) for v in word],12))
    sim.write_across_shots(c.register_data['flags'],flags);sim.write_across_shots(c.register_data['counter'],bits(shifts,6));sim.do(c)
    value_bad=np.zeros(N,dtype=bool)
    for prefix,arr in [('sigma',sigma),('prior',prior)]:
        for j,word in enumerate(arr):value_bad|=(sim.read_across_shots(c.register_data[f'{prefix}{j}'])!=bits([transform(int(v),fw) for v in word],12)).any(axis=0)
    value_bad|=(sim.read_across_shots(c.register_data['flags'])!=flags).any(axis=0)
    value_bad|=(sim.read_across_shots(c.register_data['counter'])!=bits(shifts,6)).any(axis=0)
    public=set(q for reg in c.register_data.values() for q in reg)
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    junk=sim.read_across_shots(scratch).any(axis=0);phase=sim.read_phase_flipped_across_shots()
    reports=[]
    for label in dict.fromkeys(labels):
        ids=np.array([i for i,l in enumerate(labels) if l==label])
        reports.append(dict(case=label,shots=len(ids),value_errors=int(value_bad[ids].sum()),scratch_errors=int(junk[ids].sum()),phase_errors=int(phase[ids].sum())))
    out=dict(status='DIAGNOSTIC',scope='forced exceptional polynomials outside valid split-locator promise; no repair or native full-decoder acceptance',
             seed=536193,Q=c.num_qubits,reports=reports,all_sources_restored=not value_bad.any(),all_scratch_clean=not junk.any(),all_phases_clean=not phase.any())
    Path('results/psi_exceptions.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
