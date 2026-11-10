"""Compose tested Schur components into the benchmark phase oracle.

This initial architecture leaves non-completion zero Schur pivots unsupported.
Full verifier error and repetition cost must be measured before promotion.
"""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
import kickmix as km
from schur_scan import write_scan
from extension_inverse import setup
from schur_affine_synthesize_filled import layout


def build(N=70,M=4095,repair=True):
 start=time.time();tag=f'schur_affine_filled_n{N}_m{M}';stem=Path('candidates')/tag;m,k,fw,bk,tr=setup();b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 for j in range(M):b.create_classical_register(12,f'E{j}')
 head=stem.with_suffix('.headers');head.write_text(str(b.finish_circuit())+'\n');begin=stem.with_suffix('.begin');end=stem.with_suffix('.end');g=tr(1,fw);gg=''.join(f'X q{12*N+i}\n' for i in range(12) if g>>i&1);begin.write_text(gg);end.write_text(gg)
 private=12*M+24*(N//2)*(N//2-1)//2;scan=stem.with_suffix('.scan.kmx')
 with scan.open('w') as f:L=layout(N);kernel=write_scan(f,N,M,private+12,repair,banks=(L['N'],L['B']),scratch_base=max(L['nodework'])+1)
 print(json.dumps(dict(stage='scan_text',seconds=time.time()-start,**kernel)),flush=True)
 extract_tag='schur_extract_index';synth_tag='schur_affine_synthesize_filled'
 ef=Path(f'candidates/{extract_tag}_n{N}_forward.kmx');ei=Path(f'candidates/{extract_tag}_n{N}_inverse.kmx');sf=Path(f'candidates/{synth_tag}_n{N}_forward.kmx');si=Path(f'candidates/{synth_tag}_n{N}_inverse.kmx');full=stem.with_suffix('.kmx')
 subprocess.run(['./schur_compose',str(full),str(head),str(begin),str(ef),str(12*M),str(sf),str(private),str(scan),str(si),str(private),str(ei),str(12*M),str(end)],check=True)
 print(json.dumps(dict(stage='composed_text',seconds=time.time()-start)),flush=True);c=km.Circuit.from_file(str(full),format='kmx');out=stem.with_suffix('.kmb');c.to_file(str(out),format='kmb');digest=hashlib.file_digest(out.open('rb'),'sha256').hexdigest();r=dict(status='BUILT',scope='requires complete native verification; coherent node repair' if repair else 'requires complete native verification; no zero-pivot repair',N=N,M=M,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,sha256=digest,seconds=time.time()-start);Path(f'results/{tag}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--m',type=int,default=4095);p.add_argument('--repair',action='store_true');a=p.parse_args();build(a.n,a.m,True)
