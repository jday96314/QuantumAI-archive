"""Locate the first reverse boundary whose basis values differ."""
import sys,json
from pathlib import Path
import kickmix as km
import numpy as np
path=Path(sys.argv[1]);prefix=sys.argv[2] if len(sys.argv)>2 else 'stream3'
hook=json.load(open(path.parent.parent/'parameters'/f'{prefix}_hook.json'))
sim=km.Simulator(batch_size=1,seed=8);reg=km.array([km.q(i) for i in range(hook['maxq'])])
snap={};lines=[];phase='pre';previous=None
def process():
    global previous
    if phase in ['bm/r68','bm/r69','reverse/r69']:
        print('before',phase,'lines',len(lines),'syn',np.flatnonzero(sim.read_across_shots(reg)[:840,0]).tolist(),flush=True)
    if phase.startswith('bm/r'):snap[phase[3:]]=sim.read_across_shots(reg)[:,0].copy()
    if lines:sim.do(km.Circuit(''.join(lines)))
    a=sim.read_across_shots(reg)[:,0]
    if phase.startswith('reverse/'):
        r=phase.split('/')[-1]
        expected=snap.get(r)
        if expected is not None:
            diff=np.flatnonzero(a^expected)
            print(phase,'different',len(diff),diff[:30].tolist(),flush=True)
            if len(diff):
                print('actual ones',np.flatnonzero(a).tolist(),flush=True)
                print('expected ones',np.flatnonzero(expected).tolist(),flush=True)
                closest=sorted((int(np.count_nonzero(a^v)),k) for k,v in snap.items())[:5]
                print('closest',closest,flush=True)
                raise SystemExit(1)
    previous=a
with path.open() as f:
    for line in f:
        if line.startswith('INFO_TAG "DBG:'):
            process();lines=[];phase=line.split('DBG:',1)[1].split('"')[0]
        else:lines.append(line)
process()
