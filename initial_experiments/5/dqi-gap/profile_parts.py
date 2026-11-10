"""Exact arithmetic attribution for incumbent and unified block implementations."""
import json,argparse
from pathlib import Path
from collections import defaultdict
from contextlib import contextmanager
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from profile_hybrid import DiscardTape
import fused_gap_block as gap
import blocked_monic_euclid as old


def profile(k,new):
    m=load();m.inv_word=bind(m);b=m.Builder();b._items=DiscardTape();W=[b.qubits(12) for _ in range(70)];F=b.qubits(2*(k+3)+1);Z=[(),b.qubits(12)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.virtual_zero=True;counts=defaultdict(int);stack=[];fn=gap if new else old
    def wrap(func,name):
        def call(*args,**kwargs):
            if stack:return func(*args,**kwargs)
            before=b._cur_stats['toffoli'];stack.append(name)
            try:return func(*args,**kwargs)
            finally:stack.pop();counts[name]+=b._cur_stats['toffoli']-before
        return call
    names=['cofactor','selected_cofactor','mobius'] if new else ['cofactor','cf','mobius']
    saved={name:getattr(fn,name) for name in names}
    for name in names:setattr(fn,name,wrap(saved[name],name))
    math.triangular=wrap(math.triangular,'triangular');math.triangular_suffix=wrap(math.triangular_suffix,'suffix')
    original_inv=m.inv_word
    @contextmanager
    def inv(E,word):
        begin=b._cur_stats['toffoli']
        with original_inv(E,word) as out:
            counts['inverse_gates']+=b._cur_stats['toffoli']-begin
            yield out
            begin=b._cur_stats['toffoli']
        counts['inverse_gates']+=b._cur_stats['toffoli']-begin
    m.inv_word=inv
    original_scale=m.mon_slot
    def scale(*args):
        before=b._cur_stats['toffoli'];r=original_scale(*args);counts['field_scaling']+=b._cur_stats['toffoli']-before;return r
    m.mon_slot=scale
    try:
        bp=W[:k-1];bb=W[k-1:2*k-1];c=W[2*k-1:-1];p=W[-1]
        if new:gap.apply(m,E,math,bp,bb,c,p,Z,k,3,F[:k+3],F[k+3:-1],F[-1])
        else:old.apply(m,E,math,bp,bb,c,p,3,flags=F[:3],monomial=True,suffix=True,linear_padding=True)
    finally:
        for name,func in saved.items():setattr(fn,name,func)
    total=b._cur_stats['toffoli'];exclusive={name:counts[name] for name in names+['triangular','suffix']};exclusive['other']=total-sum(exclusive.values())
    return dict(k=k,new=new,T=total,Q=b.peak_qubits,exclusive=exclusive,overlapping=dict(counts))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--k',type=int,default=10);a=p.parse_args();r=[profile(a.k,False),profile(a.k,True)];Path(f'results/profile_parts_k{a.k}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
