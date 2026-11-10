"""Read-only host and completed-reference source hashes."""
import hashlib,json
from core import ROOT,HERE,INC,SP
paths=[ROOT/'from_host/verifiers/dfthcsos/verify.py']
paths += [ROOT/f'from_host/verifiers/dfthcsos/scoreinfra/{s}.py' for s in ['reference','dfthc_math','circuits','costs_prepare']]
paths += [ROOT/'from_host/verifiers/dfthcsos/FeMoCo76o113e.fcidump',ROOT/'from_host/verifiers/dfthcsos/FeMoCo76o113e.fixed_shifts.npz',SP,INC/'submission.npz']
paths += [ROOT/f'initial_experiments/4/codex/{s}' for s in ['angle_gauge_probe.py','jacobi_l1','jacobi_l1.py','lattice_qubo.py','qubo_anneal.so','monomial_unlookup.py','build_selectswap_chem.py']]
paths += [ROOT/f'initial_experiments/5/chem-qfit/{s}' for s in ['RESULTS.md','core.py','physical.py','discrete.py','exchange.py','legacy_block_exact.py','anneal_exact.py']]
records=[]
for p in paths:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    records.append(dict(path=str(p.relative_to(ROOT)),sha256=h.hexdigest(),bytes=p.stat().st_size))
(HERE/'results/reference_provenance.json').write_text(json.dumps(records,indent=2));print(json.dumps(dict(read_only_reference_files=len(records))))
