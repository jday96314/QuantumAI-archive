"""Compose raw extraction and unnormalized CF synthesis into the phase oracle."""
import argparse,json,time,subprocess,hashlib
from pathlib import Path
import kickmix as km
from raw_cf_synthesize_destructive import layout
from raw_cf_middle_lowtemp import write_middle


def build(N=70,M=4095):
 start=time.time();n=N//2;stem=Path(f'candidates/raw_cf_destructive_n{N}_m{M}');b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 for j in range(M):b.create_classical_register(12,f'E{j}')
 head=stem.with_suffix('.headers');head.write_text(str(b.finish_circuit())+'\n');empty=stem.with_suffix('.empty');empty.write_text('');L=layout(N);private=12*M+24*n*(n-1)//2+24*n;scan=stem.with_suffix('.scan.kmx');live={q for w in L['N'][:n]+L['B'][:n] for q in w};scratch=[q for q in range(12*N+240) if q not in live]
 with scan.open('w') as f:
  f.write(f"X q{L['B'][n][0]}\n");middle=write_middle(f,N,M,private,L['N'][:n],L['B'][:n],scratch);f.write(f"X q{L['B'][n][0]}\n")
 full=stem.with_suffix('.kmx');subprocess.run(['./schur_compose',str(full),str(head),str(empty),f'candidates/schur_extract_raw_compact_n{N}_forward.kmx',str(12*M),f'candidates/raw_cf_synthesize_destructive_n{N}_forward.kmx',str(12*M+24*n*(n-1)//2),str(scan),f'candidates/raw_cf_synthesize_destructive_n{N}_inverse.kmx',str(12*M+24*n*(n-1)//2),f'candidates/schur_extract_raw_compact_n{N}_inverse.kmx',str(12*M),str(empty)],check=True);c=km.Circuit.from_file(str(full),format='kmx');out=stem.with_suffix('.kmb');c.to_file(str(out),format='kmb');r=dict(status='BUILT',scope='requires complete native verification',N=N,M=M,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,sha256=hashlib.file_digest(out.open('rb'),'sha256').hexdigest(),middle=middle,seconds=time.time()-start);Path(f'results/raw_cf_destructive_n{N}_m{M}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--m',type=int,default=4095);a=p.parse_args();build(a.n,a.m)
