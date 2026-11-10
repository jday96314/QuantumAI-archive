"""Standalone export of measured traces and rounding/refit checkpoints.

Plot dependencies are isolated in plot_deps; the verifier/training environment
is not changed. All values come from saved float64 native-metric checkpoints.
"""
from pathlib import Path
import sys,json,re,os
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'cache/matplotlib'))
sys.path.insert(0,str(HERE/'plot_deps'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter

data=[]
for p in sorted((HERE/'results').glob('scratch_*.json')):
    z=json.loads(p.read_text())
    if z.get('rows'):
        shape=re.search(r'r\db\d+',p.stem).group();reg=int(p.stem.rsplit('_l',1)[1]);data.append((shape,reg,z))
shapes=sorted({s for s,_,_ in data});colors={s:plt.get_cmap('tab10')(i) for i,s in enumerate(shapes)}
styles={30:':',300:'-',3000:'--'}
fig,axes=plt.subplots(1,2,figsize=(14,6))
left,right=axes
for shape,reg,z in data:
    rs=[r['metrics']['float'] for r in z['rows'] if r['step']>0]
    left.plot([r['error_2']*1000 for r in rs],[r['Lambda'] for r in rs],
        linestyle=styles[reg],marker='.',ms=4,lw=1.2,color=colors[shape],alpha=.8)
for shape in shapes:
    points=[r for s,_,z in data if s==shape for r in z['rows'] if r['metrics']['float']['Lambda']<=320]
    if not points:continue
    cp=min(points,key=lambda r:r['metrics']['float']['error_2'])['metrics'];stage=[]
    for kind,marker in [('float','o'),('rounded','^'),('refit','x')]:
        r=cp[kind];xy=(r['error_2']*1000,r['Lambda']);stage.append(xy)
        right.plot(*xy,marker=marker,color=colors[shape],ms=7,linestyle='none')
    for start,end in zip(stage,stage[1:]):
        right.annotate('',xy=end,xytext=start,arrowprops=dict(arrowstyle='->',color=colors[shape],alpha=.65,lw=1.1))

selected=json.loads((HERE/'selection.json').read_text())['metrics']
for ax in axes:
    ax.scatter([selected['error_2']*1000],[selected['Lambda']],marker='*',s=160,color='black',zorder=10)
    ax.axvline(20,color='black',lw=.7);ax.axhline(303.5283147431962,color='black',lw=.7)
    ax.set_xscale('log');ax.set_xlabel('Two-body Frobenius error (mHa)');ax.set_ylabel('Lambda (Ha)')
    ax.grid(alpha=.2)
left.set_title('Observed float optimization traces');right.set_title('Rounding/refit of float models with Lambda ≤ 320')
left.set_xlim(15,1000);right.set_xlim(15,260)
for ax,ticks in [(left,[20,30,50,100,200,500,1000]),(right,[20,30,50,100,200])]:
    ax.set_xticks(ticks,labels=[str(x) for x in ticks]);ax.xaxis.set_minor_formatter(NullFormatter())
handles=[Line2D([0],[0],color=colors[s],lw=2,label=s.upper().replace('B',' / B')) for s in shapes]
handles += [Line2D([0],[0],color='gray',ls=styles[r],label=f'λ_reg = {r}') for r in [30,300,3000]]
handles += [Line2D([0],[0],color='black',marker='*',linestyle='none',ms=10,label='native banked selection')]
left.legend(handles=handles,fontsize=8,loc='upper right',ncol=2)
right.legend(handles=[Line2D([0],[0],color='gray',marker=m,ls='none',label=k) for k,m in [('float','o'),('rounded','^'),('refit','x')]],loc='upper right')
fig.suptitle('Round 5 chemistry audit: measured checkpoints, fixed native shift',fontsize=13)
fig.tight_layout();fig.savefig(HERE/'results/frontier.png',dpi=180);fig.savefig(HERE/'results/frontier.svg')
print(json.dumps(dict(matplotlib=matplotlib.__version__,plots=['results/frontier.png','results/frontier.svg'])))
