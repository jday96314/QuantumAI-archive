"""Build the current exact workspace revision of the inherited approximation."""
from pathlib import Path
import argparse,json
import numpy as np
import binary_centered_invsqrt as p
from optimize_workspace import compact
root=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--output',default='candidates/sign_anf.kmb');a=parser.parse_args()
p.table_for=lambda *args:np.load(root/'candidates/margin_curve8_t5.9e-08.npy')
c=p.build(7,(27,18,8),lookup_frac=8,accum_bits=[18,9,None],carry_bits=1,
 rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,
 defer_workspace=True,qubit_limit=94,curve_keep_bits=4,defer_curve_phase=True,
 product_limit=17,borrow_integration=True,monomial_delete=True)
c=compact(c,root/a.output);c.to_file(str(root/a.output),format='kmb')
print(json.dumps(dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,output=a.output)),flush=True)
