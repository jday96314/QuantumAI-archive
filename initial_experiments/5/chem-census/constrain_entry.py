"""Read-only audit augmented-Lagrangian fit with local checkpoints/results."""
from census_common import ROOT
source=ROOT/'initial_experiments/5/chem-audit/refit.py'
text=source.read_text().replace('from audit_core import *',
    'from census_common import *\nfrom audit_core import checkpoint,write_json')
exec(compile(text,str(source),'exec'))
