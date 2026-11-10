"""Execute attributed read-only exact-word solver with strictly local outputs.

Only its star import is redirected: all proposals, acceptance logic, CLI
parameters and native physical constraints are preserved. The source digest
is part of source_manifest.json. No imported source file is modified.
"""
from census_common import ROOT
source=ROOT/'initial_experiments/5/chem-audit/integer_refine.py'
text=source.read_text().replace('from audit_core import *',
    'from census_common import *\nfrom audit_core import checkpoint,write_json')
# In bounded mode, retain the largest feasible part of a projected weight
# update. Rejecting every full least-squares projection would never test
# bounded continuation. The unchanged --unbounded path supplies a control.
helper='''
def capped_blend(uu,warm,proposal):
    if a.unbounded:return proposal,float(m.loss(uu,proposal).sqrt())
    lo,hi=0.,1.
    if not within_caps(uu,warm):return warm,float(m.loss(uu,warm).sqrt())
    if within_caps(uu,proposal):return proposal,float(m.loss(uu,proposal).sqrt())
    for _ in range(25):
        alpha=(lo+hi)/2
        if within_caps(uu,warm+alpha*(proposal-warm)):lo=alpha
        else:hi=alpha
    trial=warm+lo*(proposal-warm)
    return trial,float(m.loss(uu,trial).sqrt())
'''
text=text.replace('start=time.time();accepted=0;stale=0',helper+'\nstart=time.time();accepted=0;stale=0')
text=text.replace('if within_caps(u,wp):wt=wp','wp,pe=capped_blend(u,wt,wp)\n    if within_caps(u,wp):wt=wp')
text=text.replace('feasible=within_caps(uu,ww)','ww,ne=capped_blend(uu,wt,ww)\n            feasible=within_caps(uu,ww)')
exec(compile(text,str(source),'exec'))
