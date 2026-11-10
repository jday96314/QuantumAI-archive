"""Complete small Jacobi oracle, preserving the standard native scan contract."""
from pathlib import Path
import json
import jacobi_decoder as jacobi
from cube_arithmetic import CubeArithmetic
from import_published import load
from low_workspace_inverse import bind

HERE=Path(__file__).resolve().parent
cf=jacobi.cf
EMIT_EXTRACT=jacobi.blocked_extract
EMIT_SYNTHESIZE=jacobi.synthesize
EMIT_ARITHMETIC=CubeArithmetic


def model(F,s):
    W=list(s)
    import check_jacobi
    class Adapter:
        gmul=staticmethod(F.mul);gpow=staticmethod(F.pow)
        ginv=staticmethod(lambda a:F.inv(a) if a else 0)
    for i in range(0,len(W),4):
        if len(W)-i==2:W[i:]=check_jacobi.pair(Adapter,W[i:])
        else:W[i:]=check_jacobi.model(Adapter,W[i:])
    def fold(W):
        if len(W)==2:return [W[0]],[W[1]]
        a,b=W[:2];P,Q=fold(W[2:]);t=len(P)
        numerator=[a,*[F.mul(a or 1,q) for q in Q]]
        denominator=[b^Q[0],*[P[i]^F.mul(b,Q[i])^(Q[i+1] if i+1<t else 0) for i in range(t)]]
        return numerator,denominator
    P,Q=fold(W)
    from algebra import trim
    return trim([1,*Q]),trim(P)


def emit(N,mode='both',count_only=False):
    m=load();m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(N)]
    if count_only:b._items=jacobi.cfblock.probe.count_module().Builder()._items
    class Counter(m.Em):
        extra_T=0
        def emit_mac(self,rt):self.extra_T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=(Counter if count_only else m.Em)(b,*b.qubits(3),m.Alloc(b));math=EMIT_ARITHMETIC(m,E)
    stages=[];total=lambda:b._cur_stats['toffoli']+getattr(E,'extra_T',0)
    for name,fn in [('extract',lambda:EMIT_EXTRACT(m,E,math,W)),('synthesize',lambda:EMIT_SYNTHESIZE(m,E,W)),
                    ('synthesize_inverse',lambda:EMIT_SYNTHESIZE(m,E,W,True)),('extract_inverse',lambda:EMIT_EXTRACT(m,E,math,W,True))]:
        before=total();fn();stages.append(dict(name=name,T=total()-before,peak=b.peak_qubits))
        if name=='synthesize':hook=len(b._items)
    report=dict(N=N,Q=b.peak_qubits,T=total(),stages=stages,public_Q=12*N,scratch_Q=b.peak_qubits-12*N,
                final_live=b.active_qubits,scope='Jacobi blocked arithmetic; generic/completed charts; needs native task verification')
    return m,b,W,hook,report


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=14);p.add_argument('--m',type=int,default=127);p.add_argument('--shots',type=int,default=4096);args=p.parse_args()
    cf.HERE=HERE;cf.TRANSLATOR=HERE.parent.parent/'4/codex/translate_prefix'
    cf.emit=emit;cf.model=model;cf.synthesize=jacobi.synthesize
    cf.test(args.n,args.m,args.shots,'jacobi-oracle-a')
