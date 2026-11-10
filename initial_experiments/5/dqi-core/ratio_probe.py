"""Count a running local-ratio recurrence, including early stages."""
import json
from pathlib import Path
from import_published import load
from series_euclid_block import apply
from components.trimmed_middle import TrimmedArithmetic
from components.budget_inverse import bind

HERE=Path(__file__).resolve().parent

def main():
    rows=[]
    for k in [1,10,19,31]:
        m=load();b=m.Builder();W=[b.qubits(12) for _ in range(70)];F=b.qubits(k+3)
        class Count(m.Em):
            T=0
            def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
            def emit_free(self,rt):pass
        E=Count(b,*b.qubits(3),m.Alloc(b));m.inv_word=bind(m)
        math=TrimmedArithmetic(m,E);math.virtual_zero=True
        apply(m,E,math,W[:k-1],W[k-1:2*k-1],W[2*k-1:-1],W[-1],3,flags=F[k:])
        row=dict(k=k,block=3,forward_T=E.T+b._cur_stats['toffoli'],Q=b.peak_qubits)
        rows.append(row);print(json.dumps(row),flush=True)
    report=dict(scope='exact forward native-emitter counts, optimized full-ratio recurrence; no new full native verification',rows=rows)
    (HERE/'results/ratio_probe.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
