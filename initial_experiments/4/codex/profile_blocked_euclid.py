"""Exact arithmetic-gadget gate count without materializing Clifford tapes.

Emitter products are counted after the same term normalization as real builds.
This resource model omits complete-decoder handling, layout and phase scan.
"""
import argparse,json,time
from pathlib import Path
from import_published import load
from mobius_series import BoundedArithmetic
from blocked_monic_euclid import apply as blocked
from monic_euclid_step import apply as serial
MODULE=load()


def profile(N,k,block=None):
    m=MODULE;b=m.Builder();words=[b.qubits(12) for _ in range(N)]
    Bp=words[:k-1];B=words[k-1:2*k-1];C=words[2*k-1:-1];pinv=words[-1]
    class Emitter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(term[0]=='P' for term in rt)
        def emit_free(self,rt):pass
    sA,sB,sT=b.qubits(3);E=Emitter(b,sA,sB,sT,m.Alloc(b));math=BoundedArithmetic(m,E)
    if block is None:serial(m,E,math,Bp,B,C,pinv)
    else:blocked(m,E,math,Bp,B,C,pinv,block)
    return dict(N=N,k=k,block=block,Q=b.next_qubit,T=E.T+b._cur_stats['toffoli'])


def whole(N,width,padding=False,deferred=False,direct=False,thin=False,suffix=False):
    from blocked_euclid_sequence import frames
    from monic_euclid_local import first_step
    m=MODULE;b=m.Builder();words=[b.qubits(12) for _ in range(N)];flags=b.qubits(N//2)
    class Emitter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(term[0]=='P' for term in rt)
        def emit_free(self,rt):pass
    sA,sB,sT=b.qubits(3);E=Emitter(b,sA,sB,sT,m.Alloc(b));math=BoundedArithmetic(m,E)
    with m.eq_test(E,words[-1],0) as top:E.cx(top,flags[0])
    E.cx(flags[0],words[-1][0]);first_step(m,E,math,words)
    stages,_=frames(words,width)
    for k,count,bp,bb,c,p in stages:blocked(m,E,math,bp,bb,c,p,count,flags=flags[k:k+count],padding=padding,monomial=deferred,thin=thin,suffix=suffix)
    if deferred:
        from deferred_euclid_padding import apply,direct_repad
        postprocess=direct_repad if direct else apply
        postprocess(m,E,words,flags)
    return dict(N=N,block=width,padding=padding,deferred=deferred,direct=direct,thin=thin,suffix=suffix,Q=b.next_qubit,T=E.T+b._cur_stats['toffoli'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--all-stages',action='store_true');p.add_argument('--resume-log');p.add_argument('--whole',action='store_true');p.add_argument('--deferred',action='store_true');p.add_argument('--direct',action='store_true');a=p.parse_args();began=time.time();rows=[]
    if a.deferred:
        for width in [2,4,6,8]:
            row=whole(a.n,width,deferred=True,direct=a.direct);rows.append(row);print(json.dumps(row),flush=True)
        tag='deferred_padding'+('_direct' if a.direct else '')+f'_whole_profiles_n{a.n}'
        Path(f'results/{tag}.json').write_text(json.dumps(dict(scope='full forward compact transform with monomial core completion and deferred common-factor squarefree padding; excludes H map, H division, phase scan and repetition',rows=rows,seconds=time.time()-began),indent=2)+'\n')
        return
    if a.whole:
        for width in [2,4,6,8]:
            for padding in [False,True]:
                row=whole(a.n,width,padding);rows.append(row);print(json.dumps(row),flush=True)
        Path(f'results/blocked_whole_profiles_n{a.n}.json').write_text(json.dumps(dict(scope='full forward compact transform resource counts; excludes Wronskian, unpadding, phase scan and repetition',rows=rows,seconds=time.time()-began),indent=2)+'\n')
        return
    if a.resume_log:
        for line in Path(a.resume_log).read_text().splitlines():
            try:row=json.loads(line)
            except ValueError:continue
            if row.get('N')==a.n and 'k' in row and 'block' in row:rows.append(row);print(json.dumps(row),flush=True)
    done={(r['k'],r['block']) for r in rows}
    for k in (range(1,a.n//2) if a.all_stages else [1,5,9,13,17,21,25,29,33]):
        if k>=a.n//2:continue
        for block in [None,2,3,4,6,8]:
            if block and k+block>a.n//2:continue
            if (k,block) in done:continue
            row=profile(a.n,k,block);rows.append(row);print(json.dumps(row),flush=True)
    out=dict(scope='generic Euclid arithmetic only, exact emitter product counts; excludes degree handling, first step, scan, and repetition',N=a.n,rows=rows,seconds=time.time()-began)
    Path(f'results/blocked_euclid_profile_n{a.n}'+('_all' if a.all_stages else '')+'.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
