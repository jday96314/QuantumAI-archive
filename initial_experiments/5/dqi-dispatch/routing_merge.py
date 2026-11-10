"""Native check of the shared zero-word realignment for carry modes 1/4/5."""
import json
from pathlib import Path
import numpy as np,kickmix as km
from compact_monomial_unpad import rotate

def build(r,which):
    b=km.CircuitBuilder();words=[b.create_quantum_register(12,f'W{i}') for i in range(2*r-1)];flags=b.create_quantum_register(3,'flags')
    ops=[(words[r-2:],flags[0],-1),(words[r-3:],flags[1],-2),(words[r-1:],flags[2],2)]
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):
        for bank,flag,amount in (reversed(ops) if back else ops):rotate(b,bank,flag,amount,back)
    c=b.finish_circuit();c.to_file(f'candidates/routing_merge_r{r}_{which}.kmb',format='kmb');return c

def main(r=19,shots=4096):
    rng=np.random.default_rng(2026102531+r);source=[];target=[];flags=[]
    for i in range(shots):
        mode=i%7
        alen=[r-1,r-2,r-3,r-1,r-2,r-2,r+1][mode];blen=r if mode<3 else r-2 if mode==6 else r-1
        A=rng.integers(0,4096,alen).tolist();B=rng.integers(0,4096,blen).tolist();zeros=[0]*(2*r-1-alen-blen)
        src=A+B+zeros
        if mode==6:dst=A[:r-1]+B+A[r-1:]
        else:dst=A+[0]*(r-1-len(A))+B+[0]*(r-len(B))
        assert len(src)==len(dst)==2*r-1;source.append(src);target.append(dst);flags.append([mode in [1,4,5],mode==2,mode==6])
    source=np.array(source).T;target=np.array(target).T;flags=np.array(flags,dtype=bool).T;bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=build(r,which);sim=km.Simulator(batch_size=shots,seed=881927);sim.use_same_registers_as(c);src=target if which in ['inverse','inverse_both'] else source;dst=target if which in ['forward','inverse_both'] else source
        for i,row in enumerate(src):sim.write_across_shots(c.register_data[f'W{i}'],bits(row))
        sim.write_across_shots(c.register_data['flags'],flags);sim.do(c)
        bad=np.zeros(shots,dtype=bool)
        for i,row in enumerate(dst):bad|=(sim.read_across_shots(c.register_data[f'W{i}'])!=bits(row)).any(0)
        bad|=(sim.read_across_shots(c.register_data['flags'])!=flags).any(0);phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not bad.any() and not phase else 'FAIL',value_errors=int(bad.sum()),phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(x['status']=='PASS' for x in reports) else 'FAIL',r=r,shots=shots,circuits=reports,scope=__doc__)
    Path(f'results/routing_merge_r{r}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
