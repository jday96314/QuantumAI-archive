"""Assemble submission/ from verified per-task candidates, zip it, and score it.

Usage: package.py INV.kmb INV.json DQI.kmb DQI.json CHEM_DIR CHEM.json
Each JSON is the unchanged host verifier's output for that candidate.
"""
import hashlib, json, math, shutil, sys, zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
inv, inv_json, dqi, dqi_json, chem, chem_json = map(Path, sys.argv[1:7])
results = {k: json.loads(p.read_text()) for k, p in
           (('invsqrt', inv_json), ('dqi', dqi_json), ('dfthcsos', chem_json))}
assert all(r['status'] == 'PASS' for r in results.values())
sub = HERE / 'submission'
if sub.exists():
    shutil.rmtree(sub)
(sub / 'dfthcsos').mkdir(parents=True)
shutil.copyfile(inv, sub / 'invsqrt.kmb')
shutil.copyfile(dqi, sub / 'dqi.kmb')
for f in sorted(chem.glob('angle_k_*.kmb')) + [chem / 'submission.npz']:
    shutil.copyfile(f, sub / 'dfthcsos' / f.name)
base = {'invsqrt': 2076774.233715, 'dqi': 24439463583672.929688, 'dfthcsos': 3355599588.315537}
scores = {k: results[k]['score'] for k in base}
imp = {k: base[k] / scores[k] for k in base}
overall = math.prod(imp.values()) ** (1 / 3)
arc = HERE / 'submission.zip'
with zipfile.ZipFile(arc, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in sorted(sub.rglob('*')):
        if p.is_file():
            z.write(p, p.relative_to(sub).as_posix())
with zipfile.ZipFile(arc) as z:
    assert z.testzip() is None
out = dict(scores=scores, improvements=imp, overall_score=overall,
           sources={k: str(v) for k, v in (('invsqrt', inv_json), ('dqi', dqi_json), ('dfthcsos', chem_json))},
           sha256=hashlib.sha256(arc.read_bytes()).hexdigest())
(HERE / 'results/summary.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
