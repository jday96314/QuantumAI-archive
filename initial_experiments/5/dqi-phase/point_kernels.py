"""Point-body measurements; F/Psi preparation is deliberately separate.

Uses the read-only round-4 field and measurement-repair primitives. The
division control is the runtime-mask recursive Forney body from round 3,
derived there from the attributed Apache-2.0 Gu architecture.
"""
from contextlib import contextmanager
import json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from forney_dqi import product_terms
from recursive_forney import RecursiveForney
from cubic_inverse_product_kernel import apply as incumbent
from joint_cubic_inverse_product_kernel import apply as joint
from inverse_product_kernel import apply as old_tower
from oracle import mul, inv


@contextmanager
def root_flag(b, root):
    ladder = b.alloc_qubits(11)
    b.x(root); b.init_and(root[0], root[1], ladder[0])
    for i in range(1, 11): b.init_and(ladder[i-1], root[i+1], ladder[i])
    yield ladder[-1]
    for i in reversed(range(1, 11)): b.del_and(ladder[i-1], root[i+1], ladder[i])
    b.del_and(root[0], root[1], ladder[0]); b.x(root); b.free(ladder)


@contextmanager
def active_mask(b, mask):
    active = b.alloc_bits(1)[0]
    b.bit_store(active, True, control=mask); b.push_condition(active)
    yield
    b.pop_condition(); b.bit_store(active, False); b.free(km.array([active]))


def f_phase(b, root, value, mask):
    with active_mask(b, mask):
        parity = b.alloc_qubits(1)[0]
        for q, control in zip(value, mask):
            b.push_condition(control); b.cx(q, parity); b.pop_condition()
        with root_flag(b, root) as flag: b.cz(flag, parity)
        for q, control in zip(value, mask):
            b.push_condition(control); b.cx(q, parity); b.pop_condition()
        b.free(km.array([parity]))


def bilinear_phase(b, root, omega, psi, mask):
    terms = product_terms(km.GF2Field(12)); u = b.alloc_qubits(1)[0]
    with active_mask(b, mask), root_flag(b, root) as flag:
        for j in range(12):
            def parity():
                for k in range(12):
                    src = [omega[i] for i, jj in terms[k] if jj == j]
                    if src:
                        b.push_condition(mask[k]); b.cx(km.array(src), u); b.pop_condition()
            parity(); b.ccz(flag, u, psi[j]); parity()
    b.free(km.array([u]))


def division_phase(b, root, numerator, denominator, mask):
    class Body(RecursiveForney): gf = km.GF2Field(12)
    with active_mask(b, mask): Body()._forney_kernel(b, root, denominator, numerator, mask)


def build(kind):
    b = FastBuilder(km.CircuitBuilder())
    root = b.create_quantum_register(12, 'root'); a = b.create_quantum_register(12, 'a')
    c = None if kind == 'F' else b.create_quantum_register(12, 'c')
    mask = b.create_classical_register(12, 'mask')
    fn = dict(F=f_phase, Psi=bilinear_phase, division=division_phase,
              incumbent=incumbent, joint=joint, tower=old_tower)[kind]
    fn(b, root, a, mask) if c is None else fn(b, root, a, c, mask)
    return b.finish_circuit()


def test():
    reports = []; values = np.arange(4096); rng = np.random.default_rng(519032)
    bits = lambda x: ((x[None, :] >> np.arange(12)[:, None]) & 1).astype(bool)
    for kind in ['incumbent', 'division', 'tower', 'joint', 'Psi', 'F']:
        circ = build(kind); regs = circ.register_data; public = set(regs['root']) | set(regs['a'])
        if 'c' in regs: public |= set(regs['c'])
        scratch = km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
        cases = [(f'basis_mask_{j}', values, np.ones(4096, dtype=int), np.zeros(4096, dtype=int), np.full(4096, 1<<j)) for j in range(12)]
        cases += [('all_masks', rng.integers(0,4096,4096), rng.integers(0,4096,4096), np.zeros(4096,dtype=int), values),
                  ('nonroots', values, rng.integers(0,4096,4096), rng.integers(1,4096,4096), values),
                  ('zero_denominators', values, np.zeros(4096,dtype=int), np.zeros(4096,dtype=int), values)]
        means = []
        for seed, (name, av, cv, rv, mv) in enumerate(cases):
            sim = km.Simulator(batch_size=4096, seed=9501+seed); sim.use_same_registers_as(circ)
            for reg, v in [('root', rv), ('a', av), ('mask', mv), *([('c',cv)] if 'c' in regs else [])]: sim.write_across_shots(regs[reg], bits(v))
            sim.do(circ)
            if kind in ['incumbent','tower','joint']: vv = [inv(mul(int(a),int(c))) for a,c in zip(av,cv)]
            elif kind == 'division': vv = [mul(int(a),inv(int(c))) for a,c in zip(av,cv)]
            elif kind == 'Psi': vv = [mul(int(a),int(c)) for a,c in zip(av,cv)]
            else: vv = av
            expected = np.array([r==0 and (int(v)&int(mask)).bit_count()&1 for v,r,mask in zip(vv,rv,mv)],dtype=bool)
            assert np.array_equal(sim.read_phase_flipped_across_shots(),expected),(kind,name,'phase')
            for reg, v in [('root',rv),('a',av),('mask',mv), *([('c',cv)] if 'c' in regs else [])]:
                assert np.array_equal(sim.read_across_shots(regs[reg]),bits(v)),(kind,name,reg)
            assert not sim.read_across_shots(scratch).any(),(kind,name,'scratch')
        path = Path('candidates') / (kind+'_point.kmb'); circ.to_file(str(path),format='kmb')
        reports.append(dict(kind=kind,status='PASS',shots=len(cases)*4096,
            Q=circ.num_qubits,public_Q=len(public),scratch_Q=circ.num_qubits-len(public),
            static_T=circ.operation_statistics().num_toffoli_like,
            scope='point body only; arbitrary runtime masks, zero inputs, nonroots, source/mask restoration, scratch and phase'))
        print(json.dumps(reports[-1]),flush=True)
    Path('results/point_kernels.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__ == '__main__': test()
