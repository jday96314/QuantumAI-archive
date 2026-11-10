"""Rebuild and check the direct shear against independent native-field series."""
import argparse,json
from pathlib import Path
import numpy as np,kickmix as km
import mobius_series as old
from direct_mobius import apply
from classical_repair import mul,inv
from check_gap_block import sp,sd
from published_dqi import isomorphism
from import_published import load


def run(n,k,shots):
    # Reuse the prior component's build plumbing with all output paths local.
    old.ROOT=Path.cwd();old.apply=apply
    rng=np.random.default_rng(983351+n+k);m=load();fw,bw,tr=isomorphism(m)
    aa=rng.integers(0,4096,(k+1,shots));dd=rng.integers(0,4096,(k,shots))
    aa[0]=dd[0]=1;pp=rng.integers(1,4096,(1,shots));tt=rng.integers(0,4096,(n,shots));tt[:,::11]=0;tt[0,::3]=0
    want=np.array([sd(sp(tt[:,s].tolist(),aa[:,s].tolist(),n),[int(pp[0,s]),0]+sp(tt[:,s].tolist(),dd[:,s].tolist(),max(0,n-2)),n) for s in range(shots)]).T
    conv=lambda a:np.array([[tr(int(v),fw) for v in row] for row in a])
    aa,dd,pp,tt,want=map(conv,[aa,dd,pp,tt,want]);bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    reports=[]
    for which in ['forward','inverse','both']:
        c=old.build(n,k,which);sim=km.Simulator(batch_size=shots,seed=19161);sim.use_same_registers_as(c)
        for name,values in [('A',aa),('D',dd),('p',pp),('T',want if which=='inverse' else tt)]:
            for i,row in enumerate(values):sim.write_across_shots(c.register_data[f'{name}{i}'],bits(row))
        sim.do(c);errors={}
        for name,values in [('A',aa),('D',dd),('p',pp),('T',want if which=='forward' else tt)]:
            for i,row in enumerate(values):
                bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'{name}{i}'])!=bits(row)).any(0))
                if len(bad):errors[f'{name}{i}']=bad[:8].tolist()
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(which=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,junk=junk,phase=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    result=dict(n=n,k=k,shots=shots,status='PASS' if all(x['status']=='PASS' for x in reports) else 'FAIL',circuits=reports)
    Path(f'results/direct_mobius_n{n}_k{k}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=1024);a=p.parse_args()
    for n in [2,4,6]:
        for k in [1,3,10]:run(n,k,a.shots)
