import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km
import inspect

b = km.CircuitBuilder()
print("init_lookup:", inspect.getdoc(b.init_lookup))
print("del_lookup:", inspect.getdoc(b.del_lookup))
