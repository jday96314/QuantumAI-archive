"""Standalone scientific progress plots from preserved trajectory JSONs."""
import os,json,sys
from pathlib import Path
here=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(here/'results/mplconfig')
if (here/'results/plot_dependencies').is_dir():sys.path.insert(0,str(here/'results/plot_dependencies'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(path):return json.loads((here/path).read_text())
fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
ax=axs[0,0]
prefix='native_' if (here/'candidates/native_matched_thermal_B479.npz.json').exists() and (here/'results/native_matched_fixed_B479.json').exists() else ''
for mode in ['fixed','single','coupled','thermal']:
    p=f'results/{prefix}matched_fixed_B479.json' if mode=='fixed' else f'candidates/{prefix}matched_{mode}_B479.npz.json'
    rows=read(p);ax.plot([r['PSD_block_sweeps'] for r in rows],[r.get('best',r['error_2'])*1000 for r in rows],label=mode)
ax.set(title='B479: matched exact block eigensolves',xlabel='PSD block eigensolves',ylabel='Best two-body residual (mHa)');ax.legend()
ax=axs[0,1];r=read('candidates/coupled_B639.npz.json');s=read('candidates/long_norm_B639.npz.json');shift=r[-1]['seconds']
ax.plot([0]+[x['seconds']/60 for x in r],[102.0062267]+[x['best']*1000 for x in r],label='coupled integer')
ax.plot([(shift+x['seconds'])/60 for x in s],[x['best']*1000 for x in s],label='normalization-aware continuation')
if (here/'candidates/native_norm_B639.npz.json').exists():
    ns=read('candidates/native_norm_B639.npz.json');shift+=s[-1]['seconds']
    ax.plot([(shift+x['seconds'])/60 for x in ns],[x['best']*1000 for x in ns],label='direct native-target continuation')
ax.axhline(20,color='black',ls='--',lw=1,label='physical threshold');ax.set(title='B639: direct integer continuation',xlabel='Elapsed compute time (minutes)',ylabel='Best two-body residual (mHa)');ax.legend(fontsize=8)
ax=axs[1,0];r=read('candidates/psd_geometry_B479_full.npz.json');s=read('candidates/psd_geometry_B479_full.npz_rank_audit.json')
ax.plot([x['seconds']/60 for x in r],[x['best']*1000 for x in r],label='full PSD geometry')
ax.plot([x['seconds']/60 for x in s],[x['error_2']*1000 for x in s],'o-',label='actual C175 refits')
bound=read('results/admm_fullPSD_B479_polish_dual.json')['dual_lower_bound']*1000
ax.axhline(20,color='black',ls='--',lw=1);ax.axhline(bound,color='gray',ls=':',lw=1,label='old dictionary PSD certificate')
ax.set(title='B479: changed words and rank-C audits',xlabel='Elapsed compute time (minutes)',ylabel='Two-body residual (mHa)');ax.legend(fontsize=8)
ax=axs[1,1];r=read('candidates/feasible_B711.npz.json');s=read('candidates/final_feasible_B711.npz.json')
for rows,label in [(r,'feasible weight stage'),(s,'final exact-assignment stage')]:
    ax.plot([x['step'] for x in rows],[x['lambda_eff']*30408*478/1e6 for x in rows],label=label)
ax.axhline(708547388.9553173/1e6,color='black',ls='--',lw=1,label='selected incumbent')
ax.set(title='B711: normalization with physical feasibility',xlabel='Optimizer steps',ylabel='Estimated task cost (millions)');ax.legend(fontsize=8)
for ax in axs.flat:ax.grid(alpha=.2)
fig.suptitle('Direct integer chemistry: exact-word numerical controls',fontsize=14)
fig.savefig(here/'results/progress.png',dpi=180);fig.savefig(here/'results/progress.pdf')
print('results/progress.png and results/progress.pdf')
