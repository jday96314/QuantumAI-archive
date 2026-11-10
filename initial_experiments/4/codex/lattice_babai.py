"""Integer quadratic proposals using LLL reduction and Babai nearest planes.

The metric is real; an integer-scaled Cholesky basis approximates it. Every
proposal is evaluated again using the original metric, and chemistry acceptance
uses the actual nonlinear projected loss. This routine does not claim exact CVP.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor/lattice_python'))
from fpylll import IntegerMatrix,LLL,GSO


def proposals(H,g,max_step=32,delta=.99):
    H=.5*(H+H.T);g=np.asarray(g)
    active=np.flatnonzero(H.diagonal()>max(1e-30,float(H.diagonal().max())*1e-14))
    if not len(active):return []
    h=H[np.ix_(active,active)];gg=g[active]
    e,v=np.linalg.eigh(h);e=np.maximum(e,max(float(e.max())*1e-10,1e-30));h=(v*e)@v.T
    L=np.linalg.cholesky(h);scale=1e11/abs(L).max();B0=np.rint(L*scale).astype(np.int64)
    B=IntegerMatrix.from_matrix(B0.tolist());U=IntegerMatrix.identity(len(active));LLL.reduction(B,U,delta=delta)
    M=GSO.Mat(B,update=True);transform=np.array([[int(U[i,j]) for j in range(len(active))] for i in range(len(active))],dtype=np.int64)
    d=np.linalg.solve(h,-gg);target=tuple(float(x) for x in d@B0)
    coeff=np.array(M.babai(target),dtype=np.int64);nearest=coeff@transform
    candidates=[nearest,np.rint(d).astype(np.int64)]
    # Short reduced-basis directions can cross a barrier even if every original
    # one-coordinate step increases the quadratic loss.
    for direction in transform:
        curvature=float(direction@h@direction)
        k=int(np.rint(-float(direction@gg)/max(curvature,1e-30)))
        if k:candidates.append(k*direction)
    out=[];seen=set()
    for candidate in candidates:
        largest=max(1,int(abs(candidate).max()));trust=min(1.,max_step/largest)
        for shrink in sorted(set([1.,.5,.25,.125,trust,.5*trust,.25*trust]),reverse=True):
            move=np.rint(candidate*shrink).astype(np.int64)
            if not np.any(move) or abs(move).max()>max_step:continue
            full=np.zeros(len(g),dtype=np.int64);full[active]=move;key=tuple(full)
            if key in seen:continue
            seen.add(key);value=float(g@full+.5*full@H@full)
            if value<0:out.append((value,full))
    return sorted(out,key=lambda t:t[0])


if __name__=='__main__':
    import itertools,json
    H=np.array([[1.,.99],[.99,1.]]);g=np.array([.1,-.1])
    rows=proposals(H,g);assert rows
    brute=min((float(g@v+.5*v@H@v),v.tolist()) for v in map(np.array,itertools.product(range(-12,13),repeat=2)))
    assert np.isclose(rows[0][0],brute[0])
    assert all(float(g@v+.5*v@H@v)>0 for v in map(np.array,[(1,0),(-1,0),(0,1),(0,-1)]))
    out=dict(status='PASS',scope='correlated quadratic barrier: every unit coordinate step is uphill; lattice proposal reaches exhaustive optimum',
             proposal=rows[0][1].tolist(),objective=rows[0][0],exhaustive=brute)
    Path('results/lattice_babai_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
