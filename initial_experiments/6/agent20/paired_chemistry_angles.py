"""Pair two bilinear lookup terms with one Toffoli.

(a0+d2)(a1+d1) = a0*d1 + a1*d2 + d1*d2 + a0*a1.
The d1*d2 correction is compiled into the classical base table; a0*a1 is
shared. This is the paired-combine idea credited to Phuc Thai's public
fouraddr2 compiler, adapted to the existing compressed eight-bit bank.
Prior chemistry sources are read-only. This file modifies the generated
lookup compiler, while preserving the exact physical dictionary and words.
"""
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import shutil
import sys
import time

HERE = Path(__file__).resolve().parent
PRIOR = HERE/'../../5/chem-qrom'
sys.path.insert(0, str(PRIOR.resolve()))
import build as linear
import quadratic
import kickmix as km
import numpy as np
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas


def pair_mask(columns):
    available = set(columns)
    return sum(1 << i for m,i in columns if m == 1 and (2,i) in available)


def combine(b, address, public, extra, columns):
    positions = {tuple(key):q for key,q in zip(columns,extra)}
    paired = sorted(i for m,i in positions if m == 1 and (2,i) in positions)
    both = b.alloc_qubits(1)[0]
    b.init_and(address[0], address[1], both)
    for i in paired:
        d1,d2 = positions[1,i],positions[2,i]
        b.cx(address[0],d2)
        b.cx(address[1],d1)
        b.ccx(d2,d1,public[i])
        b.cx(address[1],d1)
        b.cx(address[0],d2)
        b.cx(both,public[i])
    controls = {1:address[0],2:address[1],3:both}
    for q,(m,i) in zip(extra,columns):
        if m in (1,2) and i in paired:
            continue
        b.ccx(controls[m],q,public[i])
    b.del_and(address[0],address[1],both)
    b.free(km.array([both]))


def patched_build(module, nonlinear):
    source = inspect.getsource(module.build)
    column = 'kept' if nonlinear else 'active'
    old = 'table=[row[0]|(v<<width) for row,v in zip(coeff,bank)]'
    new = f'table=[(row[0]^((row[1]&row[2])&pair_mask({column})))|(v<<width) for row,v in zip(coeff,bank)]'
    assert source.count(old) == 1
    source = source.replace(old,new)
    if nonlinear:
        start = source.index('    both=b.alloc_qubits(1)[0]')
        end = source.index('    decoded=b.alloc_qubits(drop)',start)
        source = source[:start]+'    combine(b,a,public,extra,kept)\n'+source[end:]
    else:
        source = source.replace('    A,J=delta.shape;', '    assert select==2 and not compress\n    A,J=delta.shape;',1)
        start = source.index('    controls={}\n')
        end = source.index('    info=erase(',start)
        source = source[:start]+'    combine(b,address,public,extra,active)\n'+source[end:]
    namespace = dict(module.__dict__)
    namespace.update(pair_mask=pair_mask,combine=combine)
    exec(source,namespace)
    name = 'nonlinear' if nonlinear else 'linear'
    (HERE/'parameters'/f'paired_chemistry_{name}_build_source.py').write_text(source)
    return namespace['build']


def main():
    prior = (HERE/'../../5/chem-census/candidates/structured').resolve()
    original = json.loads((HERE/'../../5/chem-census/results/structured_build.json').read_text())
    transform = json.loads((HERE/'../../5/chem-census/results/structured_transform.json').read_text())
    output = HERE/'candidates/paired_chemistry'
    output.mkdir(exist_ok=True)
    shutil.copyfile(prior/'submission.npz',output/'submission.npz')
    words = np.load(prior/'submission.npz',allow_pickle=False)['angle_words']
    lb,qb = patched_build(linear,False),patched_build(quadratic,True)
    rows = []
    for k,delta in enumerate(batch_deltas(words,4)):
        start = time.time()
        if k == 0:
            circuit,info = lb(delta,select=2,lowbits=7,mix=transform['matrix'])
        else:
            circuit,info = qb(delta,k,model=original['rows'][k]['info']['model'])
        assert circuit.num_qubits <= 128
        path = output/f'angle_k_{k}.kmb'
        circuit.to_file(str(path),format='kmb')
        row = dict(k=k,Q=circuit.num_qubits,T=circuit.operation_statistics().num_toffoli_like,
                   previous_T=original['rows'][k]['T'],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.time()-start)
        rows.append(row)
        print(json.dumps(row),flush=True)
    total = 19900+2*sum(row['T'] for row in rows)
    result = dict(status='BUILT_REQUIRES_NATIVE_CHECKS',rows=rows,total_Toffolis=total,
                  total_qubits=454,conditional_score=48.567337474336654*total*454,
                  physical_npz_unchanged=hashlib.sha256((prior/'submission.npz').read_bytes()).hexdigest()==hashlib.sha256((output/'submission.npz').read_bytes()).hexdigest(),
                  source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),Path(linear.__file__),Path(quadratic.__file__),HERE/'parameters/paired_chemistry_linear_build_source.py',HERE/'parameters/paired_chemistry_nonlinear_build_source.py']})
    (HERE/'results/paired_chemistry_build.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    main()
