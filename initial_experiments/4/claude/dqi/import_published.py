"""Import the Apache-2.0 decoder at andigu/uncompute, preserving attribution.

Upstream revision 64b5a13aa249f1360846657392e6c66817b26d0f.
The upstream fixed-mask sweep is replaced by a runtime-mask sweep, and
native/public field bases are converted by reversible binary linear maps.
"""
from pathlib import Path
import sys, importlib.util, json, time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/vendor/uncompute')
sys.path.insert(0,str(ROOT))
SOURCE=ROOT/'tasks/dqi-reed-solomon/circuits/discovered/circuit.py'

def load():
    spec=importlib.util.spec_from_file_location('published_decoder',SOURCE)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def capture():
    mod=load();mod.ZMASKI=[0]*70
    hook={}
    def sweep(b,E,A,omw,ebw,abw):
        hook.update(index=len(b._items),omega=[list(x) for x in omw],
                    eo=[list(x) for x in ebw],ab=[list(x) for x in abw],
                    free=list(b.free_qubits))
    mod.emit_sweep=sweep
    mod.assert_metrics=lambda b,p:b._stream_totals()[0]
    start=time.time();b,t=mod.build()
    assert not b._items._objs, [(k,type(v)) for k,v in b._items._objs.items()][:5]
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    np.save(HERE/'parameters/published_tape.npy',rows)
    hook.update(peak=b.peak_qubits,toffoli=t,ops=len(b._items),
                maxq=b.next_qubit,outcomes=b.next_out,bits=b.next_bit)
    (HERE/'parameters/published_hook.json').write_text(json.dumps(hook,indent=2))
    print({k:v for k,v in hook.items() if k not in ['omega','eo','ab','free']},
          'free wires',len(hook.get('free',[])), 'seconds',time.time()-start,flush=True)

if __name__=='__main__':capture()
