"""Write angle_k_*.kmb + submission.npz for a finalized factorization."""
import sys,shutil
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from chemistry_qrom import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
src,out=Path(sys.argv[1]),Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
for old in out.glob('angle_k_*.kmb'):old.unlink()
z=dict(np.load(src))
np.savez(out/'submission.npz',**{k:z[k] for k in ('U','W','H2_I','H1_I','Hsym_lower','tree_edges','angle_words')})
for i,delta in enumerate(batch_deltas(z['angle_words'],9)):
    c=build(delta,12);assert c.num_qubits<=128,c.num_qubits
    c.to_file(out/f'angle_k_{i}.kmb',format='kmb');print(i,c.num_qubits,c.operation_statistics().num_toffoli_like,flush=True)
