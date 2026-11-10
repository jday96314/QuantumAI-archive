"""Complete small leading-zero Jacobi oracle against its literal model."""
from pathlib import Path
import argparse,json
import check_jacobi_oracle as wrapper
import check_lookahead_jacobi as component
import lookahead_jacobi as lookahead

HERE=Path(__file__).resolve().parent


def model(F,s):
    class Adapter:
        gmul=staticmethod(F.mul);gpow=staticmethod(F.pow)
        ginv=staticmethod(lambda a:F.inv(a) if a else 0)
    W=list(s)
    for i in range(0,len(W),4):
        if len(W)-i==2:W[i:]=component.check_jacobi.pair(Adapter,W[i:])
        else:W[i:]=component.model(Adapter,W[i:])
    def one(a,b,P,Q):
        return [a,*[F.mul(a or 1,q) for q in Q]],[b^(Q[0] if Q else 0),*[p^F.mul(b,q)^(Q[i+1] if i+1<len(Q) else 0) for i,(p,q) in enumerate(zip(P,Q))]]
    def fold(W):
        if len(W)==2:return [W[0]],[W[1]]
        a,b,h,e=W[:4];P,Q=fold(W[4:]) if len(W)>4 else ([],[])
        if a:
            d=e^b;c=h^F.mul(b,d);mid=one(c,d,P,Q);return one(a,b,*mid)
        get=lambda i:Q[i] if i<len(Q) else 0
        return [0,b,*[F.mul(b or 1,q) for q in Q]],[e^get(0),h^F.mul(e,get(0))^get(1),*[p^F.mul(h,get(i))^F.mul(e,get(i+1))^get(i+2) for i,p in enumerate(P)]]
    P,Q=fold(W)
    from algebra import trim
    return trim([1,*Q]),trim(P)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=14);p.add_argument('--m',type=int,default=127);p.add_argument('--shots',type=int,default=4096);p.add_argument('--profile',action='store_true');p.add_argument('--shared',action='store_true');args=p.parse_args()
    wrapper.EMIT_EXTRACT=lookahead.extract;wrapper.EMIT_SYNTHESIZE=lookahead.synthesize;wrapper.EMIT_ARITHMETIC=lookahead.Arithmetic
    if args.shared:
        import shared_jacobi_synthesis
        wrapper.EMIT_SYNTHESIZE=shared_jacobi_synthesis.synthesize
    tag='lookahead_shared' if args.shared else 'lookahead_jacobi'
    if args.profile:
        *_,report=wrapper.emit(args.n,count_only=True)
        report['scope']='leading-zero Jacobi chart; zero second pivots remain unsupported; exact emitted resource profile'
        Path(f'results/{tag}_full_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    else:
        workspace=HERE/('lookahead_shared' if args.shared else 'lookahead')
        for sub in ('candidates','parameters','results'):(workspace/sub).mkdir(parents=True,exist_ok=True)
        cf=wrapper.cf;cf.HERE=workspace;cf.TRANSLATOR=HERE.parent.parent/'4/codex/translate_prefix'
        cf.emit=wrapper.emit;cf.model=model
        cf.test(args.n,args.m,args.shots,'lookahead-shared-oracle-a' if args.shared else 'lookahead-jacobi-oracle-a')
