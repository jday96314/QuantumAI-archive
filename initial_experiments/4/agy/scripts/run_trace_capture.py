import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3")
sys.path[:0] = [str(ROOT), str(ROOT / "vendor/uncompute")]

import import_published as p
m = p.load()

from stream_capture import capture

# Let's inspect the names allocated during BM
# We can run capture and inspect the allocator
capture('test_stream_alloc')
