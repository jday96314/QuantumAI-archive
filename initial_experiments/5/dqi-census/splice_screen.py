"""Conditional splice of inherited exact liveness/count records; no bridge assumed."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    inherited=ROOT/'../dqi-hybrid/results'
    bm=json.loads((inherited/'round3_liveness.json').read_text())['phase_statistics']['phases']
    ee=json.loads((inherited/'round4_liveness.json').read_text())['phase_statistics']['phases']
    rows=[]
    for k in [1,4,7,10,13,16,19,22]:
        end=2*k-1
        prefix=sum(p['toffoli'] for p in bm if (p['name'].startswith('bm/r') or p['name'].startswith('reverse/r')) and int(p['name'].split('/r')[1])<=end)
        q=max(p['peak'] for p in bm if p['name'].startswith('bm/r') and int(p['name'].split('/r')[1])<=end)
        suffix=2*sum(p['toffoli'] for p in ee if p['name'].startswith('forward/serial'))
        if k<22:
            suffix+=2*sum(p['toffoli'] for p in ee if p['name'].startswith('forward/block') and int(p['name'].split('block')[1])>=(k-1)//3)
        rows.append(dict(k=k,last_BM_round=end,prefix_T=prefix,prefix_peak=q,compact_suffix_T=suffix,total_T_zero_cost_handoff=prefix+suffix+370152,Q_lower_bound=max(961,q),scope='inherited exact phase counts; conditional compatible-state splice; handoff, alignment, final scale and new gap repair omitted'))
    # Explicit old scan=370152. No bridge or reconstituted final scale is free
    # in a proposed implementation; zero is only a favorable rejection screen.
    (ROOT/'results/lazy_handoff_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows))
if __name__=='__main__':main()
