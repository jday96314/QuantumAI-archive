"""Bind native evidence to delivered artifacts and prove the incumbent untouched."""
import json,hashlib,shutil
from core import *
from oracles import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
cal=json.loads((HERE/'results/calibration.json').read_text())
for name,digest in cal['source_hashes'].items():assert sha(INC/name)==digest
assert sha(ROOT/'initial_experiments/4/codex/submission.zip')=='9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68'
rebuilt=HERE/'candidates/incumbent_rebuilt';rebuilt.mkdir(exist_ok=True)
shutil.copyfile(INC/'submission.npz',rebuilt/'submission.npz')
words=np.load(INC/'submission.npz')['angle_words'];rebuild=[]
for k,delta in enumerate(batch_deltas(words,7)):
    p=rebuilt/f'angle_k_{k}.kmb';build(delta).to_file(str(p),format='kmb')
    rebuild.append(dict(k=k,identical=sha(p)==sha(INC/p.name),sha256=sha(p)))
assert all(r['identical'] for r in rebuild)
(HERE/'results/incumbent_rebuild.json').write_text(json.dumps(dict(status='PASS',circuits=rebuild,npz_identical=sha(rebuilt/'submission.npz')==sha(INC/'submission.npz')),indent=2)+'\n')
rows=[]
for dirname,report in [(str(INC),'incumbent_native'),('greedy_B711','greedy711_native'),('legacy_exact_B700_C360','b700c360_native'),('constrained_B711','constrained711_native')]:
    directory=Path(dirname) if dirname.startswith('/') else HERE/'candidates'/dirname
    result=json.loads((HERE/f'results/{report}.json').read_text())
    z=np.load(directory/'submission.npz');R,B,C=z['W'].shape
    assert submission_sha256(load_circuits(directory),directory/'submission.npz',z['tree_edges'])==result['submission_sha256']
    assert result['status']=='PASS' and result['verified_shots']==64*(76+R*B)
    assert result['resources']['angles_within_cap']
    assert result['metrics']['error_1']<=.02 and result['metrics']['error_2']<=.02
    assert result['metrics']['Lambda']>49 and result['metrics']['E_gap']>3.914898
    assert abs(result['score']-result['lambda_eff']*result['toffoli_count']*result['num_qubits'])<1e-6
    files={p.name:sha(p) for p in directory.iterdir() if p.suffix in ['.npz','.kmb']}
    rows.append(dict(directory=str(directory),native_report=report+'.json',native_submission_sha256=result['submission_sha256'],score=result['score'],files=files))
manifest=dict(status='PASS',incumbent_unchanged=True,archive_unchanged=True,submissions=rows,
              code={p.name:sha(p) for p in HERE.iterdir() if p.suffix in ['.py','.sh']})
(HERE/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
best=min(rows[1:],key=lambda r:r['score'])
(HERE/'selection.json').write_text(json.dumps(dict(replacement_selected=False,selected_source=str(INC),
    selected_score=rows[0]['score'],combined_reference=9.367722153134604,
    best_new_verified_score=best['score'],best_new_directory=best['directory'],
    reason='All new native passes cost more than the preserved incumbent; neither 567M nor 506M reached.'),indent=2)+'\n')
print(json.dumps(dict(status='PASS',incumbent_unchanged=True,rebuilt_byte_identical=True,best_new_score=best['score'])),flush=True)
