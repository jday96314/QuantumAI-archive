"""Reproduce a fixed circuit's default sample and save failed shots for diagnosis.

This does not alter the score or select hashes/salts. Saved cases are a common
regression sample for comparing actual arithmetic changes.
"""
import sys,json,inspect
from pathlib import Path
import numpy as np
import kickmix as km
sys.path.insert(0,str(Path("../../../from_host/verifiers/dqi").resolve()))
import verify_submission_dqi as v
path=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
source=inspect.getsource(v.verify_circuit) if hasattr(v,'verify_circuit') else ''
print([n for n in vars(v) if n.startswith('verify')],flush=True)
source=inspect.getsource(v.verify)
needle='        stats = sim.operation_statistics()'
extra='''        failed = val_err | phase_err
        if failed.any():
            np.savez(CAPTURE_PATH / f'{start:06d}.npz', syndromes=syndromes[failed], mbuc=mbuc[failed], phases=phases[failed], value_errors=val_err[failed], phase_errors=phase_err[failed], shot_indices=np.flatnonzero(failed)+start, sim_seed=np.array(str(sim_seed)))
'''
assert source.count(needle)==1
source=source.replace(needle,extra+needle);ns=dict(vars(v));ns['CAPTURE_PATH']=out;exec(source,ns);c=km.Circuit.from_file(str(path),format='kmb');report=ns['verify'](c,v.ProblemSpecDQI.benchmark_instance(),max_workers=3);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
