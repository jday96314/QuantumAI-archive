"""Read-only adoption of dqi-core's checked components, not its full builder."""
import importlib.util
from pathlib import Path

CORE=Path(__file__).resolve().parent/'../dqi-core'
def _read_module(name,filename):
    spec=importlib.util.spec_from_file_location(name,CORE/'components'/filename)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

TrimmedArithmetic=_read_module('dispatch_trimmed_middle','trimmed_middle.py').TrimmedArithmetic
budget_inverse=_read_module('dispatch_budget_inverse','budget_inverse.py').bind

def configure(m,limit=961):
    m.inv_word=budget_inverse(m,limit=limit)
    return TrimmedArithmetic

if __name__=='__main__':
    import sys,json
    sys.path.insert(0,str(CORE))
    spec=importlib.util.spec_from_file_location('readonly_core_attribution',CORE/'attribute_core.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    report=mod.profile(trimmed=True,retain=True)
    assert report['forward_T']==312103 and report['reverse_T']==312103
    assert report['peak_qubits']==961
    Path('results/core_components_adoption.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['forward_T','reverse_T','peak_qubits','retained_decisions']}))
