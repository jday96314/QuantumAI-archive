"""Native continuant synthesis and inverse, exact rational identity on regular heads."""
from pathlib import Path
import argparse,json,subprocess,time,hashlib
import numpy as np
import kickmix as km
from import_published import load
import projective_polynomial_inverse as projective
import quotient_jacobi_synthesis as qj
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent


def build(n,both=False):
    m=load();m.inv_word=projective.bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(2*n)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));c=E.A.open('quotient_jacobi_terminal',12);E.x(c[0])
    F,G=qj.synthesize(m,E,W,c)
    if both:qj.synthesize(m,E,W,c,True);E.x(c[0]);E.A.close('quotient_jacobi_terminal',c)
    tag=f'quotient_jacobi_n{n}'+('_both' if both else '');tape=HERE/'parameters'/f'{tag}.npy'
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
    kb=km.CircuitBuilder()
    for j in range(2*n):kb.create_quantum_register(12,f'W{j}')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    native=km.Circuit.from_file(str(out),format='kmx');native.to_file(str(out.with_suffix('.kmb')),format='kmb')
    return native,F,G,c


def test(n,shots=256):
    m=load();rng=np.random.default_rng(2026102515+n);started=time.time();reports=[]
    for both in (False,True):
        values=rng.integers(0 if both else 1,4096,(shots,2*n));c,F,G,terminal=build(n,both)
        expected=[qj.model(m,row) for row in values] if not both else None
        # Independent Jacobi recurrence cross-check of each regular rational.
        if not both:
            def product(a,b):
                out=[0]*(len(a)+len(b)-1)
                for i,x in enumerate(a):
                    for j,y in enumerate(b):out[i+j]^=m.gmul(x,y)
                return out
            for raw,(ff,gg) in zip(values,expected):
                P,Q=[],[1]
                for a,b in reversed(list(zip(raw[::2],raw[1::2]))):
                    nq=product([1,int(b)],Q)
                    for i,v in enumerate(P):nq[i+2]^=v
                    P,Q=[m.gmul(int(a),v) for v in Q],nq
                assert product(list(reversed(gg)),Q)==product(P,list(reversed(ff))),('identity',n)
        bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
        for seed in (2026102516,2026102517):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for j in range(2*n):sim.write_across_shots(c.register_data[f'W{j}'],bits(values[:,j]))
            sim.do(c)
            if both:
                for j in range(2*n):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(values[:,j])),(n,'inverse',j)
                live=set(range(24*n))
            else:
                for k,bank in enumerate((F,G)):
                    for j,word in enumerate(bank):assert np.array_equal(sim.read_across_shots(km.array([km.q(i) for i in word])),bits([row[k][j] for row in expected])),(n,k,j)
                live={i for bank in (F,G) for word in bank for i in word}
            scratch=km.array([km.q(i) for i in range(c.num_qubits) if i not in live])
            assert not sim.read_across_shots(scratch).any(),(n,both,'scratch')
            assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
        path=HERE/'candidates'/f'quotient_jacobi_n{n}{"_both" if both else ""}.kmb'
        reports.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    r=dict(status='PASS_COMPONENT',n=n,circuits=reports,seconds=time.time()-started,
           scope='generic Jacobi rational identity and continuant model; paired inverse arbitrary zero pivots; phase and scratch')
    print(json.dumps(r),flush=True);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=[3,7,17,35]);a=p.parse_args();rows=[]
    for n in a.degrees:
        rows.append(test(n));(HERE/'results/quotient_jacobi_synthesis_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
