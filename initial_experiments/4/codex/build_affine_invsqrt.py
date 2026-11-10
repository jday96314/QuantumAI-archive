"""Rebuild the fully verified affine inverse-square-root candidate."""
import argparse,hashlib,json
from pathlib import Path
from affine_invsqrt import build,compact

root=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--output',default='candidates/affine_best.kmb');args=parser.parse_args()
path=root/args.output
c=compact(build(limit=90,keep=4,delay=2,product=17,deferred=3),path)
c.to_file(str(path),format='kmb')
print(json.dumps(dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(path),
                     sha256=hashlib.sha256(path.read_bytes()).hexdigest())),flush=True)
