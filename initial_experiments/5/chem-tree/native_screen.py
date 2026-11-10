"""Full native ANGLE recovery/resources for frozen research submissions.

Uses the unchanged full native path even when physics is known to fail. Such
records are explicitly REJECT, with score infinity encoded as null in JSON.
"""
import argparse
from core import *
from dfthcsos.verify import evaluate_circuits,rejection_reason

def main():
    p=argparse.ArgumentParser();p.add_argument('submissions',nargs='+');a=p.parse_args()
    for path in a.submissions:
        src=Path(path);z=np.load(src/'submission.npz');inp=load_inputs(src/'submission.npz')
        salt='chem-tree-frozen-native-20261024-'+src.name
        r=evaluate_circuits(inp,src,z['tree_edges'],salt=salt)
        recovered=r.pop('angle_words');assert np.array_equal(recovered,z['angle_words'])
        np.savez(HERE/f'results/{src.name}_native_words.npz',angle_words=recovered)
        reason=rejection_reason(r['metrics'],r['resources'])
        r['status']='PASS' if reason is None else 'REJECT';r['reason']=reason;r['salt']=salt
        if not np.isfinite(r['score']):r['score']=None
        report(HERE/f'results/{src.name}_native.json',r)
if __name__=='__main__':main()
