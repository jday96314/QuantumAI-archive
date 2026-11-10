"""Preserve the selected approximation, share adjacent Booth products.

The frozen round-4 source is loaded read-only. All new artifacts stay local.
"""
import argparse,json,inspect
from pathlib import Path
import kickmix as km
import numpy as np
import binary_centered_invsqrt as poly
import affine_invsqrt as affine
import fixed_top_affine_invsqrt as fixed
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact
ROOT=Path.cwd().parents[2]


def build(limit=90,product=17):
 source=(ROOT/'initial_experiments/4/codex/parameters/affine_build_source.py').read_text().replace('from booth_adder import accumulate','from booth_pair_borrow_adder import accumulate')
 ns=dict(vars(poly));ns.update(table_for=lambda *args:affine.table(),restore_affine=fixed.restore_affine,cshift=zero.cshift);exec(source,ns)
 return ns['build'](7,(27,18,8),lookup_frac=8,accum_bits=[18,9,None],carry_bits=1,rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,defer_workspace=True,qubit_limit=limit,curve_keep_bits=4,defer_curve_phase=True,product_limit=product,borrow_integration=True,monomial_delete=True,delay_high_bits=2)


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=90);p.add_argument('--product',type=int,default=17);a=p.parse_args();path=Path(f'candidates/pair_booth_borrow_q{a.limit}_p{a.product}.kmb');c=compact(build(a.limit,a.product),path);c.to_file(str(path),format='kmb');r=dict(status='BUILT',scope='unchanged approximation; requires native inverse-square-root verification',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,output=str(path));Path(f'results/pair_booth_borrow_q{a.limit}_p{a.product}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
