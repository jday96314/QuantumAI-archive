"""Check submission files with the host metric's actual circuit/NPZ readers.

The two readers and their constants are loaded from its Python AST, avoiding
its unrelated pandas/orchestrator dependency. No host source is modified.
"""
import ast
from collections.abc import Mapping
import hashlib
import json
import math
from pathlib import Path
import stat
import sys
from typing import Any
import zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'from_host/kickmix'),str(ROOT/'from_host/verifiers')]


def sha256(path):
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(16*1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def main():
    tree=ast.parse((ROOT/'from_host/metric.py').read_text())
    constants={'ALLOWED_INSTRUCTIONS','MAX_QUBIT_OR_BIT_ID','MAX_NPZ_BYTES'}
    selected=[]
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in constants for t in node.targets):selected.append(node)
        if isinstance(node,ast.FunctionDef) and node.name in {'read_circuit','read_npz'}:selected.append(node)
    assert len(selected)==5
    import pathlib
    namespace={'Mapping':Mapping,'Any':Any,'np':np,'math':math,'pathlib':pathlib,'zipfile':zipfile}
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(ROOT/'from_host/metric.py'),'exec'),namespace)
    root=HERE/'submission'
    checks={}
    for path in sorted(root.rglob('*')):
        if path.is_dir():continue
        mode=path.lstat()
        assert stat.S_ISREG(mode.st_mode) and mode.st_nlink==1
        assert path.resolve().is_relative_to(root.resolve())
        if path.suffix=='.kmb':
            limit=1 if path.name=='invsqrt.kmb' else 4165 if path.name=='dqi.kmb' else 128
            c=namespace['read_circuit'](path,limit)
            if path.parent.name=='dfthcsos':assert c.num_qubits<=128
            checks[str(path.relative_to(root))]={'bytes':mode.st_size,'sha256':sha256(path),'qubits':c.num_qubits}
        else:assert path.name=='submission.npz'
    spec={'U':('f',('R','B',76)),'W':('f',('R','B','C')),'H2_I':('f',('R','C')),'H1_I':('f',()),'Hsym_lower':('f',(2926,)),'tree_edges':('iu',(75,2))}
    npz=root/'dfthcsos/submission.npz';arrays=namespace['read_npz'](npz,spec)
    checks['dfthcsos/submission.npz']={'bytes':npz.stat().st_size,'sha256':sha256(npz),'shapes':{k:list(v.shape) for k,v in arrays.items()}}
    files=set(checks)
    assert files=={'invsqrt.kmb','dqi.kmb','dfthcsos/submission.npz'}|{f'dfthcsos/angle_k_{i}.kmb' for i in range(9)}
    document={'status':'PASS','host_metric':'from_host/metric.py','files':checks}
    (HERE/'results/metric-contract.json').write_text(json.dumps(document,indent=2)+'\n')
    print('PASS: host circuit and NPZ readers, allowed instructions, file structure, width caps, and regular-file checks',flush=True)

if __name__=='__main__':main()
