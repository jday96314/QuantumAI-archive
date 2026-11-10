"""Complete compact decoder using the checked cubic trace phase kernel."""
import argparse,json,hashlib
from pathlib import Path
import compact_product_decoder as decoder
from joint_cubic_inverse_product_kernel import apply

p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--m',type=int,default=4095);p.add_argument('--output',default='candidates/joint_cubic_trace_product_n70_m4095.kmb');a=p.parse_args();decoder.phase_kernel=apply
c=decoder.build(a.n,3,22,a.m,Path(a.output),linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
paths=['build_joint_cubic_compact_decoder.py','cubic_inverse_product_kernel.py','cubic_trace_field.py','boolean_lookup.py','compact_product_decoder.py','joint_cubic_inverse_product_kernel.py','parameters/cubic_joint_quadratic_16.txt']
manifest=dict(scope='complete compact decoder with GF16 cubic trace inverse-product phase kernel; requires full native verification',parameters=vars(a),source_sha256={s:hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in paths},circuit_sha256=hashlib.sha256(Path(a.output).read_bytes()).hexdigest())
# The underlying compact decoder's pinned manifest supplies unchanged core and
# read-only upstream dependencies; include it by digest as well.
f=Path('parameters/compact_product_slim_source_manifest.json');manifest['underlying_manifest_sha256']=hashlib.sha256(f.read_bytes()).hexdigest();Path('parameters/'+Path(a.output).stem+'_source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
