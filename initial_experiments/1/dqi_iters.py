"""Build the host decoder with a configurable synchronized iteration budget."""
import argparse,json,time
from pathlib import Path
import kickmix as km
from dqi.dqi_reference_circuit import ProosAndZalkaDecoder
from dqi.verify_submission_dqi import verify,ProblemSpecDQI
p=argparse.ArgumentParser();p.add_argument('--compact',action='store_true');p.add_argument('--fast',action='store_true');p.add_argument('--tower',action='store_true');p.add_argument('--full-crt',action='store_true');p.add_argument('--crt',action='store_true');p.add_argument('--iters',type=int,required=True);p.add_argument('--shots',type=int,default=4096);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();b=km.CircuitBuilder()
if a.fast:
 from fast_dqi import FastBuilder
 b=FastBuilder(b)
elif a.tower:
 from tower_inverse import TowerBuilder
 b=TowerBuilder(b)
elif a.full_crt:
 from crt_multiplier import FullCRTBuilder
 b=FullCRTBuilder(b)
elif a.crt:
 from crt_multiplier import CRTBuilder
 b=CRTBuilder(b)
s=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
e=[b.create_classical_register(12,f'E{i}') for i in range(4095)]
if a.compact:
 from compact_dqi import CompactDecoder
 Decoder=CompactDecoder
else:Decoder=ProosAndZalkaDecoder
start=time.monotonic();Decoder(70,km.GF2Field(12)).build_syndrome_decoding_mbuc(b,syndromes=s,error_mbuc=e,n_iters=a.iters)
c=b.finish_circuit();c.to_file(a.output,format='kmb');print('BUILT',a.iters,c.num_qubits,c.operation_statistics().num_toffoli_like,'seconds',time.monotonic()-start,flush=True)
try:r=verify(c,ProblemSpecDQI.benchmark_instance(num_shots=a.shots),max_workers=4)
except Exception as ex:r={'status':'FAIL','error':str(ex)}
print(json.dumps(r,indent=2));a.output.with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n')
