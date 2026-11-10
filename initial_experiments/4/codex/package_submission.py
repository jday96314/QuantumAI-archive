"""Bind verified circuit bytes to the best local archive and its score."""
import argparse,hashlib,io,json,math,zipfile
from pathlib import Path
import numpy as np
import kickmix as km
from dfthcsos.scoreinfra.circuits import submission_sha256
BASE={'invsqrt':2076774.233715,'dqi':24439463583672.929688,'dfthcsos':3355599588.315537}
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--inverse',default='candidates/zero_head_affine.kmb');p.add_argument('--inverse-result',default='results/zero_head_affine_verify.json')
 p.add_argument('--dqi',default='candidates/cubic_trace_product_n70_m4095.kmb');p.add_argument('--dqi-result',default='results/cubic_trace_product_n70_m4095_verify.json')
 p.add_argument('--chemistry',default='candidates/norm_cycle2_signed');p.add_argument('--chemistry-result',default='results/norm_cycle2_signed_verify.json');a=p.parse_args()
 inv=Path(a.inverse);chem=Path(a.chemistry)
 paths={'invsqrt':a.inverse_result,'dqi':a.dqi_result,'dfthcsos':a.chemistry_result}
 results={k:json.loads(Path(v).read_text()) for k,v in paths.items()};assert all(r['status']=='PASS' for r in results.values())
 before=json.loads(Path('baseline/submission.summary.json').read_text());assert sha(Path('baseline/submission.zip').read_bytes())==before['sha256']
 with zipfile.ZipFile('baseline/submission.zip') as original:
  files={'dqi.kmb':Path(a.dqi).read_bytes() if a.dqi else original.read('dqi.kmb'),'invsqrt.kmb':inv.read_bytes()}
  assert sha(files['dqi.kmb'])==results['dqi']['circuit_sha256']
  # Factorizations may improve independently of the angle circuits. Their
  # exact arrays are bound to the native verifier's canonical digest below.
 assert sha(files['invsqrt.kmb'])==results['invsqrt']['circuit_sha256']
 for f in sorted(chem.glob('angle_k_*.kmb')):files[f'dfthcsos/{f.name}']=f.read_bytes()
 files['dfthcsos/submission.npz']=(chem/'submission.npz').read_bytes()
 assert len(files)==results['dfthcsos']['num_circuits']+3
 # The host digest includes canonical circuit and NPY bytes in a prescribed order.
 circuits=[km.Circuit.from_file(str(chem/f'angle_k_{k}.kmb'),format='kmb') for k in range(results['dfthcsos']['num_circuits'])]
 with np.load(chem/'submission.npz') as zz:arrays={k:zz[k] for k in ['U','W','H2_I','H1_I','Hsym_lower']};tree=zz['tree_edges']
 digest=submission_sha256(circuits,arrays,tree)
 assert digest==results['dfthcsos']['submission_sha256'],(digest,results['dfthcsos']['submission_sha256'])
 with zipfile.ZipFile('submission.zip.tmp','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as out:
  for name,data in files.items():out.writestr(name,data)
 with zipfile.ZipFile('submission.zip.tmp') as out:
  assert out.testzip() is None
  assert set(out.namelist())==set(files)
  members={name:sha(out.read(name)) for name in out.namelist()}
  assert all(members[name]==sha(data) for name,data in files.items())
 Path('submission.zip.tmp').replace('submission.zip')
 improvements={k:BASE[k]/results[k]['score'] for k in BASE};overall=math.prod(improvements.values())**(1/3)
 summary=dict(overall_score=overall,target_score=14.0,target_reached=overall>=14.0,task_scores={k:results[k]['score'] for k in BASE},improvements=improvements,sources=paths,archive='submission.zip',sha256=sha(Path('submission.zip').read_bytes()))
 Path('submission.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 Path('results/archive_audit.json').write_text(json.dumps(dict(status='PASS',archive_sha256=summary['sha256'],members=members,chemistry_digest=digest,default_shots={k:r['verified_shots'] for k,r in results.items()}),indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
