"""Native unit reciprocal, ignored constant, prefix erasure, and exact phase."""
import argparse,json,time,hashlib,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
import copied_prefix_inverse as copied
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent


def build(n,zero1=False,recursive=False):
    m=load();m.inv_word=bind(m);b=m.Builder();F=[b.qubits(12) for _ in range(n)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=copied.CopiedArithmetic(m,E,limit=10000,recursive=recursive)
    if zero1:math.known_zero.add(tuple(F[1]))
    math.inverse(F);stats=dict(Q=b.peak_qubits,T=b._cur_stats['toffoli'],copied_calls=math.copied_calls)
    tag=f'copied_inverse_n{n}_z{int(zero1)}_r{int(recursive)}';tape=HERE/'parameters'/f'{tag}.npy'
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
    kb=km.CircuitBuilder()
    for j in range(n):kb.create_quantum_register(12,f'F{j}')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb')
    stats.update(native_Q=c.num_qubits,native_T=c.operation_statistics().num_toffoli_like,sha256=hashlib.sha256(out.with_suffix('.kmb').read_bytes()).hexdigest())
    return c,stats


def test(n,zero1=False,recursive=False,shots=256):
    started=time.time();m=load();rng=np.random.default_rng(2026102512+n);a=rng.integers(0,4096,(shots,n))
    if zero1:a[:,1]=0
    expected=[]
    for f in a:
        g=[1]
        for i in range(1,n):
            v=0
            for j in range(1,i+1):v^=m.gmul(int(f[j]),g[i-j])
            g.append(v)
        g[0]=int(f[0]);expected.append(g)
    c,r=build(n,zero1,recursive);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for seed in (2026102513,2026102514):
        sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
        for i in range(n):sim.write_across_shots(c.register_data[f'F{i}'],bits(a[:,i]))
        sim.do(c)
        for i in range(n):assert np.array_equal(sim.read_across_shots(c.register_data[f'F{i}']),bits([row[i] for row in expected])),(n,zero1,recursive,i)
        assert not sim.read_phase_flipped_across_shots().any(),(n,'phase')
        scratch=km.array([km.q(i) for i in range(12*n,c.num_qubits)])
        assert not sim.read_across_shots(scratch).any(),(n,'scratch')
        sim.do(c)
        for i in range(n):assert np.array_equal(sim.read_across_shots(c.register_data[f'F{i}']),bits(a[:,i])),(n,'twice',i)
        assert not sim.read_phase_flipped_across_shots().any(),(n,'twice phase')
        assert not sim.read_across_shots(scratch).any(),(n,'twice scratch')
    r.update(status='PASS_COMPONENT',n=n,zero1=zero1,recursive=recursive,shots=2*shots,seconds=time.time()-started,
      scope='independent exact reciprocal; ignored arbitrary constant word; phase, scratch and involution')
    print(json.dumps(r),flush=True);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degrees',nargs='+',type=int,default=[6,14,34,70]);a=p.parse_args();rows=[]
    for n in a.degrees:
        for recursive in (False,True):
            for z in (False,True):
                rows.append(test(n,z,recursive));(HERE/'results/copied_prefix_inverse_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
