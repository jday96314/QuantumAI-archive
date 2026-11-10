"""Reproduce the inverse-square-root circuit with delayed phase repair.

Run from the repository root after sourcing this experiment's env.sh.
The coefficient table is fixed; all arithmetic and phase-oracle semantics
are unchanged. Private wires are recycled after resets or their final use.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import kickmix as km
import binary_centered_invsqrt as polynomial

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=Path(__file__).resolve().parent
table=np.load(root/'margin_curve8_t5.9e-08.npy')
polynomial.table_for=lambda *args:table
c=polynomial.build(7,(27,18,8),lookup_frac=8,accum_bits=[18,9,None],
                  carry_bits=1,rounded=True,booth=True,reuse_normalized_top=True,
                  output_extra_bits=10,defer_workspace=True,qubit_limit=90,
                  curve_keep_bits=4,defer_curve_phase=True,product_limit=17,
                  borrow_integration=True)
a.output.parent.mkdir(parents=True,exist_ok=True)
raw=a.output.with_suffix('.raw.kmx');compact=a.output.with_suffix('.kmx')
c.to_file(str(raw),format='kmx')
subprocess.run([str(Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/compact_wires')),str(raw),str(compact),'22'],check=True)
c=km.Circuit.from_file(str(compact));c.to_file(str(a.output),format='kmb')
print(json.dumps(dict(qubits=c.num_qubits,toffolis=c.operation_statistics().num_toffoli_like,output=str(a.output))),flush=True)
