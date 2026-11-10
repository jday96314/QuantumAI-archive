"""Measured half-degree reciprocal-derivative lifetime on compact storage.

Adapted from andigu/uncompute (Apache-2.0), revision
64b5a13aa249f1360846657392e6c66817b26d0f, disclosed circuit.py stages 5/7.
The locator recurrence is neither called nor changed. A retained 35-word
compact prior-cofactor bank is explicitly LIVE throughout this schedule.
Unlike Gu's BM layout, it is not a clean cofactor workspace.

This is a component, NOT an accepted full decoder. The upstream Euclid pass
has a limited zero-pivot repair; no new exception repair is attempted here.
"""
import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from published_dqi import isomorphism
from oracle import (EXP, add, product, derivative, evaluate, divide, inverse_mod,
                    known_error_polynomials, mul, power)

ROOT = Path(__file__).resolve().parent
OLD = ROOT / '../../4/codex'


def serialize(m, b, layout, name):
    rows = np.asarray(b._items._rows, dtype=np.uint32).reshape(-1, 9)
    tape = ROOT / 'parameters' / (name+'.npy'); np.save(tape, rows)
    pb = km.CircuitBuilder()
    for prefix, words in [('sigma', layout['sigma']), ('prior', layout['prior'])]:
        for j, word in enumerate(words):
            got = pb.create_quantum_register(12, f'{prefix}{j}')
            assert [q.id for q in got] == word
    for reg in ['flags', 'counter']:
        got = pb.create_quantum_register(len(layout[reg]), reg)
        assert [q.id for q in got] == layout[reg]
    pre = ROOT / 'parameters' / (name+'.pre.kmx'); pb.finish_circuit().to_file(str(pre), format='kmx')
    empty = ROOT / 'parameters' / 'empty'; empty.write_text('')
    path = ROOT / 'candidates' / (name+'.kmx')
    subprocess.run([str(OLD/'translate_prefix'), str(tape), str(len(rows)), '0', str(pre), str(empty), str(empty),
                    str(path), str(b.next_out), str(empty), str(len(rows))], check=True, stdout=subprocess.DEVNULL)
    circ = km.Circuit.from_file(str(path), format='kmx'); circ.to_file(str(path.with_suffix('.kmb')),format='kmb')
    # KMB plus source/layout are sufficient; retain the compact numeric tape
    # but remove the redundant large text spelling.
    path.unlink()
    return circ


def build():
    m = load(); oldswap = m.Em.swap_fn
    m.Em.swap_fn = lambda self, control, pairs: oldswap(self, control, tuple(pairs))
    b = m.Builder()
    sigma = [b.qubits(12) for _ in range(36)]
    prior = [b.qubits(12) for _ in range(35)]
    flags = b.qubits(35); counter = b.qubits(6)
    # 893 Q = 70 compact field words (including pinv), materialized lead, 35 flags
    # + six-bit unpadding counter. All those wires remain live.
    layout = dict(sigma=sigma, prior=prior, flags=flags, counter=counter)
    A = m.Alloc(b); E = m.Em(b, *b.qubits(3), A)
    EPS = A.open('zprec', 6); PRE = A.open('pre', 1)[0]
    JSET, JPAD, JPRE, JF, JSEED = (m.Journal() for _ in range(5)); PSI_IDS = set()
    with b.phase('setup'), E.using(J=JSET), A.logging(PSI_IDS):
        for word in sigma:
            JSET.step(lambda w=word: m.linear_in_place(b,w,m.SQRT_ROWS),
                      lambda w=word: m.linear_in_place(b,w,m.SQUARE_ROWS))
    EW = [sigma[2*i] for i in range(18)]; OW = [sigma[2*i+1] for i in range(18)]
    # This component's counter directly holds floor((35-L)/2). In the
    # compact decoder those bits can alias completion_counter[1:6].
    MBITS = counter[:5]
    with b.phase('pad'), E.using(J=JPAD): m.pad_barrel(E,MBITS,EW,OW)
    with b.phase('pre'), E.using(J=JPRE), A.logging(PSI_IDS): m.psi_pre(E,EW,OW,PRE)
    with b.phase('F'), E.using(J=JF), A.logging(PSI_IDS):
        FF, GG, RECS = m.psi_forward(E,EW,OW,EPS[0],EPS[1:])
    # This is the new storage that a direct transplant must actually pay.
    # It is allocated late, and is 36 full words including Gu's EXT word.
    # Allocating this bank straight from the free list aliases F's journal
    # scratch. Hold those exact IDs while allocating the persistent bank.
    with b.phase('cofactor_allocation'), A.replay(PSI_IDS,'cofactor_guard'):
        AB = [A.open('bezout_'+str(i),12) for i in range(36)]
    AW, BW = AB[:18], AB[18:]
    with b.phase('seed'), E.using(J=JSEED), A.logging(PSI_IDS): m.psi_seed(E,FF[17],MBITS,AW)
    with b.phase('R'):
        with A.replay(PSI_IDS,'psipool'):
            JR, AFIN, BFIN = m.build_R(E,JF,RECS,AW,BW,PRE,m.PSI_WIN,(EW,OW)); JR.forward()
    layout.update(A=AFIN,B=BFIN,eo=EW+OW,forward_T=b._stream_totals()[0],psi_pool=sorted(PSI_IDS))
    forward = serialize(m,b,layout,'psi_forward')
    pool = sorted(PSI_IDS)
    with b.phase('inverse_R'):
        A.reopen('psipool',pool); JR.reverse()
    with b.phase('unseed'): JSEED.reverse()
    for i, word in reversed(list(enumerate(AB))): A.close('bezout_'+str(i),word)
    with b.phase('inverse_F'): JF.reverse()
    with b.phase('inverse_pre'): JPRE.reverse()
    with b.phase('inverse_pad'): JPAD.reverse()
    with b.phase('unsetup'):
        JSET.reverse(); A.close('psipool',pool); A.close('pre',[PRE]); A.close('zprec',EPS)
    both = serialize(m,b,layout,'psi_roundtrip')
    layout.update(roundtrip_T=b._stream_totals()[0],native_Q=both.num_qubits,
                  live_peak=b.peak_qubits,peak_phase=b.peak_phase,phase_summary=b.phase_summary(),
                  public_Q=893,cofactor_words=36,cofactor_Q=432,
                  scope='half-degree Psi lifetime with live compact prior bank; no omega/layout conversion or scan')
    Path('parameters/psi_layout.json').write_text(json.dumps(layout,indent=2)+'\n')
    # Actual allocator events, including hard-ID replay holds, are exported.
    Path('results/psi_liveness.json').write_text(json.dumps(dict(peak=b.peak_qubits,peak_phase=b.peak_phase,
        events=b.wire_log,phases=b.phase_summary()),indent=2)+'\n')
    return forward,both,layout


def test(forward,both,layout,shots=1024):
    rng = np.random.default_rng(530190); m=load(); fw,bk,transform=isomorphism(m)
    sigmas=[];priors=[];shifts=[];known=[];weights=[]
    for shot in range(shots):
        L = shot % 36
        pos = rng.choice(4095,L,replace=False).tolist(); vals=rng.integers(1,4096,L).tolist()
        sig,om,syn = known_error_polynomials(pos,vals)
        # sigma may be scaled, as long as no zero scale is introduced.
        scale=int(rng.integers(1,4096));sig=[mul(x,scale) for x in sig]
        sigmas.append(sig+[0]*(36-len(sig)));priors.append(rng.integers(0,4096,35).tolist())
        shifts.append((35-L)//2);known.append((pos,vals));weights.append(L)
    sigmas=np.array(sigmas).T;priors=np.array(priors).T;shifts=np.array(shifts)
    bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    reports=[]
    for name,circ in [('forward',forward),('roundtrip',both)]:
        sim=km.Simulator(batch_size=shots,seed=70193);sim.use_same_registers_as(circ)
        for prefix,values in [('sigma',sigmas),('prior',priors)]:
            for j,word in enumerate(values):sim.write_across_shots(circ.register_data[f'{prefix}{j}'],bits(np.array([transform(int(v),fw) for v in word]),12))
        fv=rng.integers(0,2,(35,shots)).astype(bool)
        sim.write_across_shots(circ.register_data['flags'],fv)
        sim.write_across_shots(circ.register_data['counter'],bits(shifts,6));sim.do(circ)
        def read(words):
            return [[transform(int(v),bk) for v in sum(sim.read_across_shots(km.array([km.q(q) for q in word])).astype(np.int64)<<np.arange(12)[:,None])]
                    for word in words]
        assert np.array_equal(read(layout['prior']),priors),(name,'prior')
        assert np.array_equal(sim.read_across_shots(circ.register_data['flags']),fv),(name,'flags')
        assert np.array_equal(sim.read_across_shots(circ.register_data['counter']),bits(shifts,6)),(name,'counter')
        failures=[]
        if name=='forward':
            aa=np.array(read(layout['A']));bb=np.array(read(layout['B']))
            for s in range(shots):
                # Psi = B^2 + z*A^2; a coefficient square is binary linear.
                psi=[0]*36
                for i in range(18):psi[2*i]=mul(int(bb[i,s]),int(bb[i,s]));psi[2*i+1]=mul(int(aa[i,s]),int(aa[i,s]))
                pos,_=known[s]
                for j in pos:
                    x=power(2,-j)
                    if mul(evaluate(psi,x),evaluate(derivative(sigmas[:,s]),x))!=1:
                        failures.append(dict(shot=s,degree=weights[s],position=j));break
            # A nonzero phase in a compute-only circuit is also an error.
            phase=np.flatnonzero(sim.read_phase_flipped_across_shots()).tolist()
            reports.append(dict(direction=name,shots=shots,wrong_reciprocals=failures,phase_errors=phase))
        else:
            values_bad=np.flatnonzero((np.array(read(layout['sigma']))!=sigmas).any(axis=0)).tolist()
            public=set(q for key in ['sigma','prior'] for word in layout[key] for q in word)|set(layout['flags'])|set(layout['counter'])
            sc=km.array([km.q(i) for i in range(circ.num_qubits) if i not in public])
            junk=np.flatnonzero(sim.read_across_shots(sc).any(axis=0)).tolist();phase=np.flatnonzero(sim.read_phase_flipped_across_shots()).tolist()
            reports.append(dict(direction=name,shots=shots,value_errors=values_bad,scratch_errors=junk,phase_errors=phase))
    out=dict(status='PASS' if not any(r.get(k) for r in reports for k in ['wrong_reciprocals','phase_errors','value_errors','scratch_errors']) else 'FAIL',
        scope=layout['scope'],seed=530190,Q=both.num_qubits,live_peak=layout['live_peak'],peak_phase=layout['peak_phase'],
        forward_T=layout['forward_T'],roundtrip_T=layout['roundtrip_T'],public_Q=893,cofactor_Q=432,psi_pool_Q=len(layout['psi_pool']),reports=reports)
    Path('results/psi_schedule.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=1024);p.add_argument('--reuse',action='store_true');a=p.parse_args()
    if a.reuse:
        f=km.Circuit.from_file('candidates/psi_forward.kmb',format='kmb');b=km.Circuit.from_file('candidates/psi_roundtrip.kmb',format='kmb');l=json.loads(Path('parameters/psi_layout.json').read_text())
    else:f,b,l=build()
    test(f,b,l,a.shots)
