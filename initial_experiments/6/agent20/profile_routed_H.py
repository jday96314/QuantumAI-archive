"""Exact accepted-core arithmetic count with local routed scale controls."""
from pathlib import Path
import json
import routed_quantum_control as routing
import H_decoder as decoder
routing.install()
HERE=Path(__file__).resolve().parent
build=decoder.builder(profile=True)
report=build(70,3,100,4095,HERE/'candidates/routed_H_profile_only.kmb',linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
report={k:v for k,v in report.items() if k!='rows'};report['scope']='exact routed-scale arithmetic profile; no native circuit or submission acceptance'
(HERE/'results/routed_H_core_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
