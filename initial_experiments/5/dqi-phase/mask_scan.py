"""Native runtime mask interpolation + prepared-polynomial phase scan.

Gu's interpolation identity is applied to per-shot classical masks, without
build-time parity choices. Exactly 70 positions can always be skipped. All
other zero residuals are tested at runtime; a 71st skip is not hard-coded.
This verifies scan and syndrome restoration, not polynomial preparation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from tower_inverse import apply_map, map_ops, linear_xor
from oracle import (EXP, MASK_TO_DUAL, DUAL_TO_MASK, mul, power, trace,
                    interpolation_inverse, inverse_mod, derivative,
                    product, divide, known_error_polynomials,mask_residual)
from point_kernels import f_phase, bilinear_phase


def classical_scale_xor(b,source,target,coefficient):
    linear_xor(b,source,target,[mul(coefficient,1<<j) for j in range(12)])


def prepare_masks(b,masks,inverse):
    dual_ops=map_ops([int(MASK_TO_DUAL[1<<i]) for i in range(12)])
    for mask in masks:apply_map(b,mask,dual_ops)
    u=[b.alloc_bits(12) for _ in range(70)]
    for i in range(70):
        for j in range(70):classical_scale_xor(b,masks[j],u[i],int(inverse[i,j]))
    # Use ascending/descending ordered XORs with distinct source/target
    # registers; the same operation order reversed restores arbitrary masks.
    for j,mask in enumerate(masks):
        for i,word in enumerate(u):classical_scale_xor(b,word,mask,power(2,(i+1)*j))
    def cleanup():
        for j in reversed(range(len(masks))):
            for i in reversed(range(70)):classical_scale_xor(b,u[i],masks[j],power(2,(i+1)*j))
        for i in reversed(range(70)):
            for j in reversed(range(70)):classical_scale_xor(b,masks[j],u[i],int(inverse[i,j]))
        for word in u:b.free(word)
        for mask in masks:apply_map(b,mask,dual_ops,True)
    return u,cleanup


def zlinear(b,syndromes,u):
    for syn,word in zip(syndromes,u):
        for i,control in enumerate(word):
            targets=[syn[j] for j in range(12) if trace(mul(1<<i,1<<j))]
            if targets:
                b.push_condition(control);b.z(km.array(targets));b.pop_condition()


def build(kind,M):
    assert M>=70
    b=FastBuilder(km.CircuitBuilder());syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    sig=[b.create_quantum_register(12,f'sigma{i}') for i in range(36)]
    payload=[b.create_quantum_register(12,f'payload{i}') for i in range(35)]
    psi=[b.create_quantum_register(12,f'psi{i}') for i in range(35)] if kind=='Psi' else []
    masks=[b.create_classical_register(12,f'E{j}') for j in range(M)]
    u,cleanup=prepare_masks(b,masks,interpolation_inverse());zlinear(b,syn,u)
    mask_ops=map_ops([int(DUAL_TO_MASK[1<<i]) for i in range(12)])
    root=b.alloc_qubits(12);value=b.alloc_qubits(12);other=b.alloc_qubits(12) if psi else None
    for j in range(70,M):
        x=power(2,-j)
        def evals():
            for i,w in enumerate(sig):classical_scale_xor(b,w,root,power(x,i))
            for i,w in enumerate(payload):classical_scale_xor(b,w,value,power(x,i))
            for i,w in enumerate(psi):classical_scale_xor(b,w,other,power(x,i))
        evals();apply_map(b,masks[j],mask_ops)
        if psi:bilinear_phase(b,root,value,other,masks[j])
        else:f_phase(b,root,value,masks[j])
        apply_map(b,masks[j],mask_ops,True);evals()
    if other is not None:b.free(other)
    b.free(value);b.free(root);cleanup()
    return b.finish_circuit()


def test(M=128,shots=512,reuse=False):
    rng=np.random.default_rng(532921);cases=[]
    for s in range(shots):
        L=s%36;pos=rng.choice(M,L,replace=False).tolist();vals=rng.integers(1,4096,L).tolist()
        sig,om,syn=known_error_polynomials(pos,vals);ps=inverse_mod(derivative(sig),sig);F=divide(product(om,ps),sig)[1]
        masks=rng.integers(0,4096,M)
        if s%31==0:masks[:]=0
        if s%31==1:masks[:]=4095
        phase=sum((e&int(masks[j])).bit_count() for j,e in zip(pos,vals))&1
        cases.append(dict(sigma=sig,payload=om,F=F,psi=ps,S=syn,E=masks,want=phase))
    bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    vi=interpolation_inverse()
    active_count=sum(np.count_nonzero(mask_residual(row['E'],vi)[1][70:]) for row in cases)
    reports=[]
    for kind in ['F','Psi']:
        path=f'candidates/runtime_mask_{kind}_m{M}.kmb'
        if reuse:c=km.Circuit.from_file(path,format='kmb')
        else:c=build(kind,M);c.to_file(path,format='kmb')
        sim=km.Simulator(batch_size=shots,seed=9412,counting_operations=True);sim.use_same_registers_as(c);sources={}
        for prefix,count in [('sigma',36),('payload',35),*([('psi',35)] if kind=='Psi' else []),('S',70),('E',M)]:
            key='F' if prefix=='payload' and kind=='F' else prefix
            for i in range(count):
                values=np.array([row[key][i] if i<len(row[key]) else 0 for row in cases]);sources[f'{prefix}{i}']=values
                sim.write_across_shots(c.register_data[f'{prefix}{i}'],bits(values))
        sim.do(c)
        assert np.array_equal(sim.read_phase_flipped_across_shots(),np.array([r['want'] for r in cases],dtype=bool)),kind
        for name,v in sources.items():assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),(kind,name)
        public=set(q for name,reg in c.register_data.items() if not name.startswith('E') for q in reg)
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),kind
        body=11 if kind=='F' else 23
        assert c.operation_statistics().num_toffoli_like==body*(M-70)
        stats=sim.operation_statistics()
        assert stats.num_toffoli_like==body*active_count,(kind,stats.num_toffoli_like,body*active_count)
        reports.append(dict(status='PASS',kind=kind,shots=shots,M=M,all_locator_degrees=list(range(36)),
            Q=c.num_qubits,public_Q=len(public),scratch_Q=c.num_qubits-len(public),static_T=body*(M-70),guaranteed_skips=70,
            native_mean_T=stats.num_toffoli_like/shots,native_active_point_shots=int(active_count),native_num_bits=c.num_bits,
            scope='prepared-polynomial scan, runtime classical interpolation, linear syndrome phase, full 70-syndrome/mask/source restoration; excludes polynomial preparation'))
        print(json.dumps(reports[-1]),flush=True)
    Path('results/mask_scan.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--m',type=int,default=128);p.add_argument('--shots',type=int,default=512);p.add_argument('--reuse',action='store_true');a=p.parse_args();test(a.m,a.shots,a.reuse)
