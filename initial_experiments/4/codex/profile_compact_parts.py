"""Attribute exact compact-transform counts to its arithmetic components."""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import json
from pathlib import Path

import blocked_monic_euclid as block
from blocked_euclid_sequence import frames
from deferred_euclid_padding import apply as postprocess
from import_published import load
from mobius_series import BoundedArithmetic
from monic_euclid_local import first_step


def profile(N,width):
    m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)];flags=b.qubits(N//2)
    class Counter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(term[0]=='P' for term in rt)
        def emit_free(self,rt):pass
    E=Counter(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E)
    counts=defaultdict(int);stack=[]
    total=lambda:E.T+b._cur_stats['toffoli']
    @contextmanager
    def region(name):
        before=total();stack.append(name)
        try:yield
        finally:
            stack.pop();counts[name]+=total()-before
    def wrapped(fn,name,outer_only=False):
        def call(*args,**kwargs):
            if outer_only and stack:return fn(*args,**kwargs)
            with region(name):return fn(*args,**kwargs)
        return call
    originals={key:getattr(block,key) for key in ['cofactor','cf','mobius']}
    for key,name in [('cofactor','cofactor_three_passes'),('cf','head_continued_fraction'),('mobius','head_mobius')]:
        setattr(block,key,wrapped(originals[key],name))
    math.triangular=wrapped(math.triangular,'tail_triangular',True)
    rows=[]
    try:
        with region('first_step'):
            with m.eq_test(E,words[-1],0) as top:E.cx(top,flags[0])
            E.cx(flags[0],words[-1][0]);first_step(m,E,math,words)
        stages,_=frames(words,width)
        for k,count,bp,bb,c,p in stages:
            before=total();prior=dict(counts)
            block.apply(m,E,math,bp,bb,c,p,count,flags=flags[k:k+count],monomial=True)
            row=dict(k=k,block=count,T=total()-before)
            row['parts']={name:value-prior.get(name,0) for name,value in counts.items() if value-prior.get(name,0)}
            row['parts']['core_padding_and_other']=row['T']-sum(row['parts'].values())
            rows.append(row)
        with region('final_common_padding'):postprocess(m,E,words,flags)
    finally:
        for key,fn in originals.items():setattr(block,key,fn)
    parts=dict(counts);parts['core_padding_and_other']=total()-sum(counts.values())
    return dict(N=N,block=width,Q=b.next_qubit,T=total(),parts=parts,rows=rows,
                scope='exact forward component count; excludes Wronskian, H division, scan and repetition')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--block',type=int,default=4);a=p.parse_args()
    out=profile(a.n,a.block)
    Path(f'results/compact_parts_n{a.n}_b{a.block}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='rows'}),flush=True)
