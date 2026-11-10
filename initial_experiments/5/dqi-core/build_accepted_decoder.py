"""Build the selected schedule with trimmed middle products and retained inverses.

All generated files are confined to this directory. The copied selected builder
changes only its arithmetic class, inverse binding, and read-only translator
path; the original cubic scan and all schedule parameters are preserved.
Field arithmetic: andigu/uncompute, Apache-2.0, revision
64b5a13aa249f1360846657392e6c66817b26d0f; see the read-only vendor LICENSE.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import kickmix as km
from components import trimmed_decoder as decoder
from cubic_inverse_product_kernel import apply

HERE=Path(__file__).resolve().parent
R4=HERE.parents[1]/'4/codex'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk:=f.read(8<<20):h.update(chunk)
    return h.hexdigest()

def main():
    begin=time.time();os.chdir(HERE)
    for path in ['candidates','results','parameters']:Path(path).mkdir(exist_ok=True)
    decoder.phase_kernel=apply
    raw=HERE/'candidates/trimmed_decoder_raw.kmb'
    c=decoder.build(70,3,22,4095,raw,linear_padding=True,virtual_zero=True,
                    virtual_lead=True,low_inverse=True,lazy_flags=True)
    assert c.operation_statistics().num_toffoli_like==994358
    del c
    binary=HERE/'components/compact_wires'
    subprocess.run(['g++','-O3',str(HERE/'components/compact_wires.cpp'),'-o',str(binary)],check=True)
    compact=HERE/'candidates/trimmed_decoder.kmx'
    with (HERE/'results/trimmed_compaction.log').open('w') as log:
        subprocess.run([str(binary),str(raw.with_suffix('.kmx')),str(compact),'840'],check=True,stderr=log)
    c=km.Circuit.from_file(str(compact),format='kmx')
    assert c.num_qubits==961,c.num_qubits
    assert c.operation_statistics().num_toffoli_like==994358
    target=HERE/'candidates/trimmed_decoder.kmb';c.to_file(str(target),format='kmb')
    report=dict(status='BUILT_REQUIRES_FULL_VERIFICATION',Q=c.num_qubits,
        static_T=c.operation_statistics().num_toffoli_like,core_T=624206,scan_T=370152,
        N=70,M=4095,block=3,serial_after=22,circuit_sha256=sha(target),
        incumbent_circuit_sha256='4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0',
        sources={str(p.relative_to(HERE)):sha(p) for p in [Path(__file__),HERE/'components/trimmed_decoder.py',HERE/'components/trimmed_middle.py',HERE/'components/budget_inverse.py',HERE/'components/compact_wires.cpp']},
        inherited_sources={str(p.relative_to(R4)):sha(p) for p in [R4/'compact_product_decoder.py',R4/'cubic_inverse_product_kernel.py',R4/'parameters/compact_product_slim_source_manifest.json']},
        seconds=time.time()-begin)
    (HERE/'parameters/trimmed_decoder_source_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'results/trimmed_decoder_build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)

if __name__=='__main__':main()
