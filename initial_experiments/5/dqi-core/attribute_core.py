"""Exact emitted-Toffoli accounting without materializing Clifford tapes.

Reuses round-5 dqi-hybrid replay and round-4 profile_compact_parts methodology.
The arithmetic is the Apache-2.0 andigu/uncompute implementation, revision
64b5a13aa249f1360846657392e6c66817b26d0f, via read-only import_published.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import json
from pathlib import Path
import time
import sys
from import_published import load
import blocked_monic_euclid as block
import monomial_serial_euclid as serial
from blocked_euclid_sequence import frames
from monic_euclid_local import first_step
from mobius_series import BoundedArithmetic

HERE=Path(__file__).resolve().parent

def profile(width=3, serial_after=22, n=70, both=True, trimmed=False, retain=False):
    start=time.time();m=load();b=m.Builder()
    words=[b.qubits(12) for _ in range(n)];alloc=m.Alloc(b)
    flags=list(alloc.open('first_flag',1));stack=[];counts=Counter();inclusive=Counter();calls=[]
    arithmetic_calls=[];arithmetic_depth=0;current_name='init'
    class CountEmitter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=CountEmitter(b,*b.qubits(3),alloc)
    from low_workspace_inverse import bind
    lowinv=bind(m)
    retained_decisions=Counter()
    @contextmanager
    def budget_inverse(E,source):
        # chain_live holds 12 output + 18 intermediate bits. Compare the
        # actual live allocation, not the historical physical frontier.
        use=retain and b.active_qubits+42<=961
        retained_decisions['retained' if use else 'low_workspace']+=1
        if use:
            with m.opened(E,'retained_inverse_output',12) as out:
                with m.chain_live(E,source,out):yield out
        else:
            with lowinv(E,source) as out:yield out
    m.inv_word=budget_inverse
    if trimmed:
        from components.trimmed_middle import TrimmedArithmetic
        math=TrimmedArithmetic(m,E)
    else:math=BoundedArithmetic(m,E)
    math.virtual_zero=True
    total=lambda:E.T+b._cur_stats['toffoli']
    @contextmanager
    def region(name):
        before=total();stack.append([name,0])
        try:yield
        finally:
            name,children=stack.pop();delta=total()-before
            counts[name]+=delta-children
            if not any(frame[0]==name for frame in stack):inclusive[name]+=delta
            if stack:stack[-1][1]+=delta
    def wrap(fn,name):
        def call(*args,**kwargs):
            with region(name):return fn(*args,**kwargs)
        return call
    def ctxwrap(fn,name):
        @contextmanager
        def call(*args,**kwargs):
            cm=fn(*args,**kwargs)
            with region(name):value=cm.__enter__()
            try:yield value
            except BaseException:
                with region(name):suppress=cm.__exit__(*sys.exc_info())
                if not suppress:raise
            else:
                with region(name):cm.__exit__(None,None,None)
        return call
    m.inv_word=ctxwrap(m.inv_word,'inverses')
    m.eq_test=ctxwrap(m.eq_test,'flags_completion')
    m.mon_slot=wrap(m.mon_slot,'scaling')
    def arithmetic_wrap(fn,suffix=False):
        def call(*args,**kwargs):
            nonlocal arithmetic_depth
            before=total();top=arithmetic_depth==0
            parent=[frame[0] for frame in stack]
            arithmetic_depth+=1
            try:
                with region('series_division_middle'):return fn(*args,**kwargs)
            finally:
                arithmetic_depth-=1
                if top:
                    ii=4 if suffix else 3
                    inverse=args[ii] if len(args)>ii else kwargs.get('inverse',False)
                    arithmetic_calls.append(dict(region=current_name,T=total()-before,
                        operation=('suffix_' if suffix else '')+('division' if inverse else 'multiply'),
                        coefficients=len(args[0])-1,q=args[1],target_words=len(args[2]),
                        area='head' if any(name.startswith('six_moment_') for name in parent) else 'tail'))
        return call
    math.triangular=arithmetic_wrap(math.triangular)
    math.triangular_suffix=arithmetic_wrap(math.triangular_suffix,True)
    originals={key:getattr(block,key) for key in ['cofactor','cf','mobius']}
    for key,cat in [('cofactor','cofactor_updates'),('cf','six_moment_CF'),('mobius','six_moment_Mobius')]:
        setattr(block,key,wrap(originals[key],cat))
    # Padding matrices are entirely controlled linear circuits.
    import controlled_linear
    origlinear=controlled_linear.apply
    controlled_linear.apply=wrap(origlinear,'padding')
    import monic_euclid_step,monic_euclid_local
    inverse_originals=[]
    for module in [monic_euclid_step,monic_euclid_local]:
        for key in ['inverse_keep','replace_inverse']:
            fn=getattr(module,key);inverse_originals.append((module,key,fn))
            setattr(module,key,wrap(fn,'inverses'))
    origgeneric=serial.generic;serial.generic=wrap(origgeneric,'cofactor_updates')
    stages,final=frames(words,width);bp,bb,_,p=final
    def run(name,fn):
        nonlocal current_name
        current_name=name
        before=total();prior=counts.copy();peak_before=b.peak_qubits
        with region('other'):fn()
        calls.append(dict(name=name,T=total()-before,
                          parts=dict(counts-prior),peak_so_far=b.peak_qubits))
    def first(inverse=False):
        if not inverse:
            with m.eq_test(E,words[-1],0) as z:E.cx(z,flags[0])
            E.cx(flags[0],words[-1][0])
        first_step(m,E,math,words,inverse)
        if inverse:
            E.cx(flags[0],words[-1][0])
            with m.eq_test(E,words[-1],0) as z:E.cx(z,flags[0])
    def scale(inverse=False):
        with m.inv_word(E,p) as pi,m.opened(E,'scan_scale_product',12) as temp:
            for word in bp:m.mon_slot(E,word,pi if inverse else p,p if inverse else pi,temp)
    def stage(k,count,oldp,old,c,pin,inverse=False):
        if k>=serial_after:
            if inverse:
                ff,_=block.schedule(oldp,old,list(reversed(c))[:2*count])
                for j in reversed(range(count)):
                    prev,cur,_,_=ff[j]
                    serial.apply(m,E,math,prev,cur,c[:len(c)-2*j],pin,flags[k+j],True,linear_padding=True)
            else:
                for j in range(count):oldp,old,c,pin=serial.apply(m,E,math,oldp,old,c,pin,flags[k+j],linear_padding=True)
        else:block.apply(m,E,math,oldp,old,c,pin,count,inverse,flags=flags[k:k+count],monomial=True,suffix=True,linear_padding=True)
    try:
        run('forward/first',first)
        for args in stages:
            k,count,*_=args;flags.extend(alloc.open(f'flags_{k}',count))
            run(f'forward/k{k}/b{count}',lambda args=args:stage(*args))
        run('forward/scale',scale);forward=total();forward_parts=dict(counts);forward_inclusive=dict(inclusive)
        if both:
            run('reverse/scale',lambda:scale(True))
            for args in reversed(stages):
                k,count,*_=args
                run(f'reverse/k{k}/b{count}',lambda args=args:stage(*args,True))
                alloc.close(f'flags_{k}',flags[k:k+count])
            run('reverse/first',lambda:first(True));alloc.close('first_flag',flags[:1])
    finally:
        for key,fn in originals.items():setattr(block,key,fn)
        controlled_linear.apply=origlinear
        serial.generic=origgeneric
        for module,key,fn in inverse_originals:setattr(module,key,fn)
    reverse_parts={key:counts[key]-forward_parts.get(key,0) for key in counts}
    assert sum(counts.values())==total()
    return dict(n=n,width=width,serial_after=serial_after,forward_T=forward,reverse_T=total()-forward,
                forward_parts=forward_parts,reverse_parts=reverse_parts,peak_qubits=b.peak_qubits,
                wire_frontier=b.next_qubit,calls=calls,seconds=time.time()-start,
                forward_inclusive=forward_inclusive,retained_decisions=dict(retained_decisions),
                arithmetic_calls=arithmetic_calls,
                trimmed=trimmed,retain=retain,
                method='exact product-term and direct-builder CCX accounting; Clifford emission suppressed')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--block',type=int,default=3);ap.add_argument('--serial-after',type=int,default=22)
    ap.add_argument('--forward-only',action='store_true');ap.add_argument('--trimmed',action='store_true');ap.add_argument('--retain',action='store_true');a=ap.parse_args()
    report=profile(a.block,a.serial_after,both=not a.forward_only,trimmed=a.trimmed,retain=a.retain)
    if a.block==3 and a.serial_after==22 and not (a.trimmed or a.retain):
        assert report['forward_T']==324901,report['forward_T']
        if not a.forward_only:assert report['reverse_T']==324901
    suffix=('_trimmed' if a.trimmed else '')+('_retain' if a.retain else '')
    path=HERE/f'results/attribution_b{a.block}_s{a.serial_after}{suffix}.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='calls'}),flush=True)
