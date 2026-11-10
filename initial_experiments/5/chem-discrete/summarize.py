"""Recompute final exact-word checkpoint metrics and bind native artifacts."""
import json,hashlib,argparse,io
from scipy.optimize import linear_sum_assignment
from core import *
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x

constants=np.load(HERE/'candidates/native_constants.npz');h0=constants['h0'];constant=float(constants['constant'])
parser=argparse.ArgumentParser();parser.add_argument('--trajectories',action='store_true');parser.add_argument('--repair',action='store_true');parser.add_argument('--records-only',action='store_true');args=parser.parse_args()
rows=[];m=None
existing={}
if args.repair and (HERE/'results/numerical_summary.json').exists():
    existing={r['checkpoint']:r for r in json.loads((HERE/'results/numerical_summary.json').read_text())}
paths=list((HERE/'candidates').glob('*.npz'))
if args.trajectories:paths+=list((HERE/'trajectories').glob('*.npz'))
if args.records_only:
    rows=json.loads((HERE/'results/numerical_summary.json').read_text());paths=[]
for path in sorted(paths):
    payload=path.read_bytes();sha=hashlib.sha256(payload).hexdigest();rel=str(path.relative_to(HERE))
    if path.parent.name=='trajectories' and rel in existing and existing[rel]['sha256']==sha:
        rows.append(existing[rel]);continue
    with np.load(io.BytesIO(payload)) as data:
        if 'U' not in data or 'angle_words' not in data:continue
    z,un,wn,tree,words=load(io.BytesIO(payload))
    if m is None:m=Metric(tree)
    assert np.array_equal(tree,m.tree)
    R,B,C=wn.shape;t=m.tensor(words.reshape(R,B,75));u=m.decode(t);w=m.tensor(wn)
    err=float(m.loss(u,w).sqrt());nativeu=reconstruct_u_vectors(words/128,tree).reshape(R,B,76)
    stored=un/np.linalg.norm(un,axis=-1,keepdims=True);sgn=np.where((nativeu*stored).sum(-1)<0,-1.,1.)
    du=float(abs(nativeu-sgn[...,None]*stored).max());assert du<1e-10,(path,du)
    visited={0}
    for p,c in tree:assert p in visited and c not in visited;visited.add(int(c))
    assert words.dtype.kind in 'iu' and words.min()>=0 and words.max()<256
    row=dict(checkpoint=rel,shape=[R,B,C],brot=8,error_2=err,decoded_vector_max_error_up_to_sign=du,sha256=sha)
    if 'W' in z and z['W'].shape==(R,C):
        off=z['W'];h=h0-np.einsum('rc,rcpq->pq',off,reconstruct_x(nativeu,wn));ev,V=np.linalg.eigh((h+h.T)/2)
        lam=float(abs(ev).sum()+.25*np.square(abs(wn).sum(1)+abs(off+wn.sum(1))).sum());gap=float(constant+abs(ev).sum()-ev.sum()+.5*np.square(off).sum())
        row.update(Lambda=lam,E_gap=gap,lambda_eff=float(np.sqrt(gap*(2*lam-gap))) if 0<gap<2*lam else None)
        if 'onebody_angle_words' in z:
            one=reconstruct_u_vectors(z['onebody_angle_words']/128,tree);one/=np.linalg.norm(one,axis=-1,keepdims=True)
            rr,cc=linear_sum_assignment(-abs(one@V));matched=np.empty(76);matched[rr]=ev[cc]
            row['error_1']=float(np.linalg.norm((one.T*matched)@one-h))
            row['physical_pass']=bool(err<=.020 and row['error_1']<=.020 and lam>49 and 3.914898<gap<2*lam)
    rows.append(row)
    if len(rows)%20==0:print(json.dumps(dict(checkpoints=len(rows),last=row['checkpoint'])),flush=True)
(HERE/'results/numerical_summary.json').write_text(json.dumps(rows,indent=2))
physical=[]
for path in sorted((HERE/'results').glob('*_physical.json')):
    r=json.loads(path.read_text())
    if 'metrics' not in r and 'error_1' in r:
        native=json.loads((HERE/'results/incumbent_native.json').read_text())
        r=dict(fit=str(INC/'submission.npz'),shape=[1,719,350],metrics=r,physical_pass=True,native_status='PASS',resources_static=native['resources'],estimated_score=native['score'],native_record='results/incumbent_native.json')
    if 'resources_static' in r:r['brot']=r['resources_static']['brot']
    tag=Path(r['fit']).stem;nativepath=HERE/f'results/{tag}_native.json'
    if nativepath.exists() and 'circuit_directory' in r:
        native=json.loads(nativepath.read_text());r['native_status']=native['status'];r['native_score']=native['score'];r['native_record']=str(nativepath.relative_to(HERE))
    r['record_file']=str(path.relative_to(HERE));physical.append(r)
(HERE/'results/physical_summary.json').write_text(json.dumps(physical,indent=2))
native=[]
for path in sorted((HERE/'results').glob('*_native.json')):
    r=json.loads(path.read_text());r['record_file']=str(path.relative_to(HERE));native.append(r)
selection=dict(status='no_replacement',incumbent_cost=708547388.9553173,best_new_native_cost=min([r['score'] for r in native if r['record_file']!='results/incumbent_native.json'],default=None),target_110_cost=708547388.9553173/1.1,target_120_cost=708547388.9553173/1.2,target_135_cost=708547388.9553173/1.35,native_results=native)
assert not any(r['score']<selection['incumbent_cost']-1 for r in native),'A winning native candidate needs promotion.'
(HERE/'selection.json').write_text(json.dumps(selection,indent=2));print(json.dumps(dict(numerical_checkpoints=len(rows),physical_records=len(physical),native_records=len(native),selection=selection['status'])),flush=True)
