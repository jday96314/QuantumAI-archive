"""Find where the 530 workspace qubits are allocated in stream_capture."""
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3")
sys.path[:0] = [str(ROOT), str(ROOT / "vendor/uncompute")]

import import_published as p
m = p.load()

# Let's inspect the allocations in BM
print("B:", m.B)
print("CHAIN_WINS:", m.CHAIN_WINS)
