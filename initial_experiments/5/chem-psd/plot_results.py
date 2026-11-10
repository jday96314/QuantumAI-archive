"""Standalone rank/error and geometry convergence artifacts (mHa)."""
import os,json,csv
from pathlib import Path
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(name):return json.loads((HERE/f'results/{name}.json').read_text())
def main():
    records=[];fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,name,title in zip(axes.flat,['prior_r4b180','prior_r4b224','prior_r2b351','claude_r2b351'],
           ['Prior R4/B180','Prior R4/B224','Prior R2/B351','Stronger R2/B351 basin']):
        if name.startswith('prior'):
            xs=[];ys=[]
            for f in (HERE/'results').glob(name+'_c*.json'):
                z=json.loads(f.read_text());xs.append(z['rank']);ys.append(1000*z['diagnostics']['direct_error_2'])
            full=load(name+'_psd')['diagnostics'];pairs=sorted(zip(xs,ys));ax.plot(*zip(*pairs),'o-',label='Rank-C feasible fits')
            for c,e in pairs:records.append(dict(dictionary=name,quantized=False,C=c,error_mHa=e))
        else:
            best={}
            for nm in ['claude_r2b351_cont','claude_r2b351_warm_cont','claude_r2b351_threshold']:
                if not (HERE/f'results/{nm}_capacity.json').exists():continue
                for r in load(nm+'_capacity')['ranks']:best[r['C']]=min(best.get(r['C'],float('inf')),1000*r['direct_error_2'])
            pairs=sorted(best.items());ax.plot(*zip(*pairs),'o-',label='Best continuous rank fit')
            for c,e in pairs:records.append(dict(dictionary=name,quantized=False,C=c,error_mHa=e))
            full=load('claude_r2b351_psd')['diagnostics']
            quant=load('claude_r2b351_q8_capacity');pairs=[(r['C'],1000*r['direct_error_2']) for r in quant['ranks']]
            ax.plot(*zip(*pairs),'s-',label='Decoded eight-bit rank fit')
            ax.axhspan(1000*quant['full_PSD']['dual_error_lower_bound'],1000*quant['full_PSD']['direct_error_2'],color='tab:orange',alpha=.25,label='Decoded full-PSD interval')
            for c,e in pairs:records.append(dict(dictionary=name,quantized=True,C=c,error_mHa=e))
        ax.axhspan(1000*full['dual_error_lower_bound'],1000*full['direct_error_2'],color='tab:blue',alpha=.2,label='Continuous full-PSD interval')
        ax.axhline(20,color='black',ls='--',lw=1,label='Physical two-body threshold')
        ax.set(title=title,xlabel='Columns per block C',ylabel='Two-body error (mHa)',yscale='log');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Fixed-dictionary rank capacity: finite fits and convex primal/dual intervals')
    fig.savefig(HERE/'results/rank_curves.png',dpi=180);fig.savefig(HERE/'results/rank_curves.pdf');plt.close(fig)
    with (HERE/'results/rank_curves.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for name,label in [('trust_r4b224_cpu','R4/B224 continuation'),('trust_prior_r2b351_cpu','R2/B351 continuation'),
                        ('trust_eigen_r2b351','Individual eigenatoms R2/B351'),('trust_graph_r4b224','Shared covariance graph R4/B224')]:
        z=load(name);hist=z['history'];axes[0].plot([0]+[r['step']+1 for r in hist],[1000*z['initial']['direct_error_2']]+[1000*r['error_2'] for r in hist],label=label)
    for name,label in [('integer_claude_r2b351','Gradient, 145 mHa start'),('integer_newton_r2b351','Newton, 68 mHa start'),('integer_projected_r2b351','Projected Newton, 68 mHa start')]:
        z=load(name);hist=z['history'];axes[1].plot([0]+[r['step']+1 for r in hist],[1000*z['initial']['direct_error_2']]+[1000*r['error_2'] for r in hist],label=label)
    for ax,title in zip(axes,['Continuous full-PSD geometry searches','Exact decoded integer-vector searches']):
        ax.axhline(20,color='black',ls='--',lw=1);ax.set(title=title,xlabel='Accepted geometry iteration within run',ylabel='Two-body error (mHa)',yscale='log');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(HERE/'results/geometry_curves.png',dpi=180);fig.savefig(HERE/'results/geometry_curves.pdf')
    (HERE/'results/plot_captions.json').write_text(json.dumps(dict(rank_curves='Finite rank fits are upper bounds; full-PSD bands are fixed-dictionary numerical primal/dual intervals. All errors use native signed target.',geometry_curves='Each run has its own explicitly labelled initialization and finite budget. These are unscored physical residuals, not accepted task costs.'),indent=2)+'\n')
if __name__=='__main__':main()
