"""Standalone research figures using read-only completed plotting dependencies."""
from pathlib import Path
import os,sys,json,csv
HERE=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(HERE/'.matplotlib')
import numpy as np
sys.path.insert(0,str(HERE.parent/'chem-psd/.plotdeps'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def write(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(HERE/f'results/{name}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def main():
    controls=['incumbent','b639_float','b351_float','b479_exact'];rows=[];maps=[]
    for c in controls:
        data=json.loads((HERE/f'results/screen_{c}.json').read_text());rows.extend(data);maps.append({r['family']:r for r in data})
    keys=list(maps[0]);values=np.array([[m[k]['error_refit']*1000 for m in maps] for k in keys])
    fig,ax=plt.subplots(figsize=(8,9));im=ax.imshow(values,norm=LogNorm(vmin=19,vmax=values.max()),cmap='magma_r',aspect='auto')
    ax.set_xticks(range(4),['R1/B719/C350\naccepted','R1/B639/C350\ncontinuous','R2/B351/C189\ncontinuous','R2/B479/C175\nexact']);ax.set_yticks(range(len(keys)),keys)
    for i in range(len(keys)):
        for j in range(4):ax.text(j,i,f'{values[i,j]:.1f}',ha='center',va='center',fontsize=8,color='white' if values[i,j]>150 else 'black')
    ax.set_title('Exact decoded tensor error after matched 12-cycle refits (mHa)')
    fig.colorbar(im,ax=ax,label='Two-body error (mHa, logarithmic colors)',shrink=.6);write(fig,'tree_screens')
    with (HERE/'results/tree_screens.csv').open('w') as f:
        cols=['control','family','error_refit','frozen_error','vector_rms','linear_tensor_delta','expected_isotropic_tensor_rms','seconds']
        writer=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for tag,label in [('topology_b351','Adaptive topology'),('random_b351','Random balanced trees')]:
        data=json.loads((HERE/f'results/{tag}.json').read_text());axs[0].plot([r['step']+1 for r in data],[r['best']*1000 for r in data],label=label)
    axs[0].set(xlabel='Shared-tree proposals (24 block updates each)',ylabel='Best two-body error (mHa)',title='B351/C189: matched topology-search work')
    for tag,label in [('repair_topology_b351','Selected new tree + integer moves'),('repair_incumbent_b351','Incumbent tree + integer moves')]:
        path=HERE/f'results/{tag}.json'
        if path.exists():
            data=json.loads(path.read_text());axs[1].plot([r['block_eigensolves'] for r in data],[r['error_2']*1000 for r in data],label=label)
    for tag,label in [('fixed_topology_b351','New tree: weights only'),('fixed_incumbent_b351','Incumbent tree: weights only')]:
        path=HERE/f'results/{tag}.json'
        if path.exists():
            data=json.loads(path.read_text());axs[1].axhline(data['final']*1000,ls=':',label=label,color='tab:blue' if 'topology' in tag else 'tab:orange')
    axs[1].set(xlabel='Block updates within integer-repair stage',ylabel='Two-body error (mHa)',title='Conditional repair comparison')
    for ax in axs:
        ax.axhline(20,color='black',ls='--',lw=1,label='Physical threshold');ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.text(.5,-.04,'The new tree required a separate 5,760-update topology search; repair curves do not include that work.',ha='center',fontsize=9)
    fig.tight_layout();write(fig,'search_trajectories')
if __name__=='__main__':main()
