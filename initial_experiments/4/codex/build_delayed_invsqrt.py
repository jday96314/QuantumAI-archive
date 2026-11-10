"""Rebuild the delayed-coefficient inverse-square-root candidate."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
from optimize_workspace import compact
root=Path(__file__).resolve().parent
a=argparse.ArgumentParser();a.add_argument('--output',default='candidates/delayed_rebuilt.kmb');args=a.parse_args()
p.table_for=lambda *unused:np.load(root/'candidates/margin_curve8_t5.9e-08.npy')
c=p.build(7,(27,18,8),lookup_frac=8,accum_bits=[18,9,None],carry_bits=1,rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,defer_workspace=True,qubit_limit=90,curve_keep_bits=4,defer_curve_phase=True,product_limit=17,borrow_integration=True,monomial_delete=True,delay_high_bits=4)
path=root/args.output;c=compact(c,path);c.to_file(str(path),format='kmb')
print(json.dumps(dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())),flush=True)
