"""Unchanged full-default native evaluation retaining rejected research metrics.

The host verify() calls these same native functions, then raises on a physical
rejection. This recorder retains their complete result including default shots;
it never overrides thresholds, scoring, seeding, or native oracle semantics.
Passing candidates are additionally run with the host verify.py CLI.
"""
import argparse
from census_common import *
from dfthcsos.verify import evaluate_circuits,rejection_reason,load_tree_edges
from dfthcsos.scoreinfra.circuits import load_circuits,SHOTS_PER_ADDRESS
p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--salt',required=True);p.add_argument('--out',required=True);a=p.parse_args()
path=Path(a.candidate);inp=load_inputs(path/'submission.npz');edges=load_tree_edges(None,path)
cs=load_circuits(path);r=evaluate_circuits(inp,cs,edges,salt=a.salt)
words=r.pop('angle_words');expected=np.load(path/'submission.npz')['angle_words']
assert np.array_equal(words,expected)
r['recovered_words_match_saved']=True
reason=rejection_reason(r['metrics'],r['resources'])
assert r['verified_shots']==inp.num_vectors*SHOTS_PER_ADDRESS
status='PASS' if reason is None else 'REJECT'
if not np.isfinite(r['score']):r['score']=None
r.update(status=status,rejection_reason=reason,salt=a.salt,default_shots_per_address=SHOTS_PER_ADDRESS,
    conditional_objective=r['metrics']['lambda_eff']*r['resources']['toffoli']*r['resources']['qubits'])
write(a.out,r);print(json.dumps(r,allow_nan=False),flush=True)
