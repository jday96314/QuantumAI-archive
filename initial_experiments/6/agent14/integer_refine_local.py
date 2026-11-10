"""Run the attributed integer proposal solver with strictly local outputs."""
import sys
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
source=ROOT/'initial_experiments/5/chem-audit/integer_refine.py'
exec(compile(source.read_text(),str(source),'exec'))
