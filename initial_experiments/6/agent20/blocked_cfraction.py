"""Two-pair C-fraction extraction with one long normalized reciprocal.

For head C-fraction coefficients c0..c3, let d=c1+c2, e=d+c3,
h=c1*c3, u=c0*c1, p=u*c2*c3, v=u*d, w=u*d^2+p.
A=1+e*x+h*x^2, D=1+d*x, P=(e*w+h*v)+(h*w)*x.
The remaining tail R maps to A*T/(p+x*D*T), T=A*R+P.
Zero p is armed to one to define a reversible total map. Supported generic
and completed charts are checked independently; other charts are not claimed.
"""
from contextlib import contextmanager
import importlib
import json
import sys
from pathlib import Path

import dqi_budget_probe as probe
from cube_arithmetic import CubeArithmetic
from import_published import load
from low_workspace_inverse import bind

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(probe.ROUND5/'dqi-core-phase'))
cf=importlib.import_module('cfraction_decoder')


@contextmanager
def parameters(m,E,math,H):
    c0,c1,c2,c3=H
    with m.opened(E,'block_d',12) as d,m.opened(E,'block_e',12) as e,m.opened(E,'block_h',12) as h,m.opened(E,'block_p',12) as p,m.opened(E,'block_v',12) as v,m.opened(E,'block_w',12) as w:
        def linear(src,dst):
            return [('L',(a,),0,(b,)) for a,b in zip(src,dst)]
        def prod(a,b,out):return m.prod_rt(out,a,b)
        dr=linear(c1,d)+linear(c2,d);er=linear(d,e)+linear(c3,e);hr=prod(c1,c3,h)
        E.mac(dr);E.mac(er);E.mac(hr)
        with m.opened(E,'block_u',12) as u:
            ur=prod(c0,c1,u);E.mac(ur)
            with m.opened(E,'block_q',12) as q:
                qr=prod(c2,c3,q);pr=prod(u,q,p)
                E.mac(qr);E.mac(pr);E.free(qr)
            vr=prod(u,d,v);E.mac(vr)
            wr=[];m.bil_terms(wr,m.KT,m.spec_plain(u),math.spec([(1,d,2)]),tgt=w)
            wr+=linear(p,w);E.mac(wr);E.free(ur)
        with m.opened(E,'block_P0',12) as p0,m.opened(E,'block_P1',12) as p1:
            r0=prod(e,w,p0)+prod(h,v,p0);r1=prod(h,w,p1)
            E.mac(r0);E.mac(r1)
            yield [(),e,h],[(),d],p,(p0,p1)
            E.free(r1);E.free(r0)
        with m.opened(E,'block_u',12) as u:
            ur=prod(c0,c1,u);E.mac(ur)
            wr=[];m.bil_terms(wr,m.KT,m.spec_plain(u),math.spec([(1,d,2)]),tgt=w);wr+=linear(p,w)
            E.free(wr);E.free(prod(u,d,v))
            with m.opened(E,'block_q',12) as q:
                qr=prod(c2,c3,q);E.mac(qr);E.free(prod(u,q,p));E.free(qr)
            E.free(ur)
        E.free(hr);E.free(er);E.free(dr)


def block(m,E,math,W,inverse=False):
    assert len(W)>=4 and len(W)%2==0
    H=W[:4];R=W[4:]
    def head(back=False):
        if back:
            cf.peel_pair(m,E,math,H[2:],True);cf.peel_pair(m,E,math,H,True)
        else:
            cf.peel_pair(m,E,math,H);cf.peel_pair(m,E,math,H[2:])
    if not inverse:head()
    if R:
        with parameters(m,E,math,H) as (A,D,p,P):
            def shift():
                for src,dst in zip(P,R):
                    for a,b in zip(src,dst):E.cx(a,b)
            if inverse:
                math.triangular(A,1,R,True)
                math.triangular(D,1,R)
                math.inverse([(),*R])
                cf.scale(m,E,R,p)
                math.triangular(D,1,R,True)
                shift();math.triangular(A,1,R,True)
            else:
                math.triangular(A,1,R);shift()
                math.triangular(D,1,R)
                cf.scale(m,E,R,p,True)
                math.inverse([(),*R])
                math.triangular(D,1,R,True)
                math.triangular(A,1,R)
    if inverse:head(True)


def extract(m,E,math,W,inverse=False):
    steps=list(range(0,len(W),4))
    for k in (reversed(steps) if inverse else steps):
        if len(W)-k==2:cf.peel_pair(m,E,math,W[k:],inverse)
        else:block(m,E,math,W[k:],inverse)


def profile(N=70):
    m=load();m.inv_word=bind(m);b=m.Builder();b._items=probe.count_module().Builder()._items
    W=[b.qubits(12) for _ in range(N)]
    class Counter(m.Em):
        extra=0
        def emit_mac(self,rt):self.extra+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=Counter(b,*b.qubits(3),m.Alloc(b));math=CubeArithmetic(m,E)
    reports=[]
    for name,fn in [('extract',lambda:extract(m,E,math,W)),('synthesize',lambda:cf.synthesize(m,E,W)),
                    ('synthesize_inverse',lambda:cf.synthesize(m,E,W,True)),('extract_inverse',lambda:extract(m,E,math,W,True))]:
        before=b._cur_stats['toffoli']+E.extra;fn()
        reports.append(dict(stage=name,T=b._cur_stats['toffoli']+E.extra-before,Q=b.peak_qubits))
    return dict(N=N,stages=reports,Q=b.peak_qubits,T=sum(r['T'] for r in reports),
                whole_oracle_static_Toffolis=sum(r['T'] for r in reports)+80*4095)


if __name__=='__main__':
    r=profile();Path('results/blocked_cfraction_profile.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
