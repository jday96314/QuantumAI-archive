"""Bounded-storage, measurement-correct arithmetic DAG evaluation.

Cache eviction first measures a derived word, then reconstructs its direct
parents to repay the phase. Thus nonadjacent repeated polynomial subexpressions
can share a small cache without retaining every nonlinear intermediate.
"""
from contextlib import contextmanager
from collections import Counter


class Evaluator:
    def __init__(self,m,E,math,compiler,variables,capacity=10):
        self.m,self.E,self.math=m,E,math
        self.nodes=compiler.nodes;self.operations=compiler.operations;self.variables=variables
        self.capacity=capacity;self.cache={};self.pins=Counter();self.tick=0;self.age={}
        self.products=0;self.evictions=0;self.peak=0;self.busy=set();self.building=set();self.serial=0
        self.ancestors={}
        def ancestry(i):
            if i not in self.ancestors:self.ancestors[i]=frozenset({i}).union(*(ancestry(j) for j in self.deps(i)))
            return self.ancestors[i]
        for i in range(len(self.nodes)):ancestry(i)

    def deps(self,node):
        kind,*args=self.nodes[node]
        if kind=='mul' or kind=='add':return args
        if kind=='square':return args[:1]
        return []

    def touch(self,node):self.tick+=1;self.age[node]=self.tick

    @contextmanager
    def value(self,node):
        kind,*args=self.nodes[node]
        if kind=='var':yield self.m.spec_plain(self.variables[args[0]]);return
        if kind=='zero':yield [frozenset() for _ in range(12)];return
        if kind=='one':raise ValueError('constant one must simplify before emission')
        if kind=='square':
            with self.value(args[0]) as source:yield list(map(frozenset,self.math.squared_spec(source)))
            return
        self.ensure(node);self.pins[node]+=1;self.touch(node)
        try:yield self.m.spec_plain(self.cache[node][1])
        finally:self.pins[node]-=1

    @contextmanager
    def direct(self,node):
        from contextlib import ExitStack
        kind,*args=self.nodes[node]
        with ExitStack() as stack:
            values=[stack.enter_context(self.value(i)) for i in args]
            if kind=='add':
                spec=[frozenset() for _ in range(12)]
                for value in values:spec=[a^b for a,b in zip(spec,value)]
                yield kind,spec
            elif kind=='mul':yield kind,values
            else:raise ValueError(kind)

    def terms(self,kind,value,target):
        if kind=='add':return [('L',tuple(sorted(s)),0,(t,)) for s,t in zip(value,target) if s]
        out=[];self.m.bil_terms(out,self.m.KT,*value,tgt=target);return out

    def phase(self,terms,masks):
        b=self.E.b
        for term in terms:
            for target in term[-1]:
                cond=masks[target]
                if term[0]=='L':
                    for source in term[1]:b.z(source,cond=cond)
                    if term[2]:b.neg(cond=cond)
                else:self.E.pivot_pair(term[1],term[2],term[3],term[4],lambda a,c,cond=cond:b.cz(a,c,cond=cond))

    def remove(self,node):
        assert not self.pins[node]
        name,word=self.cache.pop(node);self.busy.add(node)
        masks={q:self.E.b.hmr(q) for q in word}
        self.E.A.close(name,word)
        # Parent reconstruction may itself evict words, but ancestors cannot
        # depend on this measured descendant.
        with self.direct(node) as (kind,value):self.phase(self.terms(kind,value,word),masks)
        self.busy.remove(node);self.evictions+=1

    def room(self):
        while len(self.cache)>=self.capacity:
            choices=[i for i in self.cache if not self.pins[i] and i not in self.busy and not(self.ancestors[i]&self.building)]
            if not choices:
                # A wide XOR has several simultaneous parents. The actual
                # high-water mark is reported, never hidden by a hard claim.
                return
            # Prefer cheap-to-recreate expressions, with recency as tie break.
            victim=min(choices,key=lambda i:(len(self.operations[i]),self.age[i]))
            self.remove(victim)

    def ensure(self,node):
        if node in self.cache:return
        self.building.add(node)
        with self.direct(node) as (kind,value):
            self.room();self.serial+=1;name=f'field_dag_{self.serial}'
            word=self.E.A.open(name,12);self.cache[node]=(name,word);self.touch(node)
            terms=self.terms(kind,value,word);self.E.mac(terms)
            self.products+=sum(t[0]=='P' for t in terms)
            self.peak=max(self.peak,len(self.cache))
        self.building.remove(node)

    def xor(self,node,target):
        with self.value(node) as value:self.E.xor([('L',tuple(sorted(s)),0,(t,)) for s,t in zip(value,target) if s])

    def clear(self):
        while self.cache:
            # Descendants first avoid recreating ancestors that are still live.
            node=max(self.cache,key=lambda i:len(self.operations[i]));self.remove(node)


class TreeEvaluator(Evaluator):
    """Recompute nonlinear XOR children for the final measurement repayment.

    Pure binary linear combinations stay virtual. A sum with several nonlinear
    children is accumulated one child at a time, so its full expression tree
    need not remain live at once.
    """
    @contextmanager
    def value(self,node):
        from contextlib import ExitStack
        kind,*args=self.nodes[node]
        if kind in ('var','zero','square','one'):
            with super().value(node) as value:yield value
            return
        if kind=='add' and sum(bool(self.operations[j]) for j in args)<=1:
            with ExitStack() as stack:
                values=[stack.enter_context(self.value(i)) for i in args]
                out=[frozenset() for _ in range(12)]
                for value in values:out=[a^b for a,b in zip(out,value)]
                yield out
            return
        self.serial+=1;name=f'field_tree_{self.serial}'
        with self.m.opened(self.E,name,12) as word:
            if kind=='mul':
                with self.value(args[0]) as left,self.value(args[1]) as right:
                    terms=self.terms('mul',(left,right),word);self.E.mac(terms)
                    self.products+=sum(t[0]=='P' for t in terms)
                    yield self.m.spec_plain(word)
                    self.E.free(terms)
            elif kind=='add':
                for j in args:
                    with self.value(j) as value:self.E.xor(self.terms('add',value,word))
                yield self.m.spec_plain(word)
                masks={q:self.E.b.hmr(q) for q in word}
                for j in args:
                    with self.value(j) as value:self.phase(self.terms('add',value,word),masks)
            else:raise ValueError(kind)
    def clear(self):pass


def profile(pairs,capacity=10,tree=False):
    import dqi_budget_probe as probe
    from import_published import load
    from cube_arithmetic import CubeArithmetic
    from continuant_polynomials import derive,Compiler
    m=load();b=m.Builder();b._items=probe.count_module().Builder()._items
    inputs=[b.qubits(12) for _ in range(2*pairs)];outputs=[b.qubits(12) for _ in range(pairs+1)]
    class CounterEmitter(m.Em):
        extra=0
        def emit_mac(self,rt):self.extra+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=CounterEmitter(b,*b.qubits(3),m.Alloc(b));math=CubeArithmetic(m,E)
    polys=derive(pairs);compiler=Compiler(2*pairs);roots=[compiler.compile(polys['p']),*[compiler.compile(p) for p in polys['P']]]
    evaluator=(TreeEvaluator if tree else Evaluator)(m,E,math,compiler,inputs,capacity)
    for node,out in zip(roots,outputs):evaluator.xor(node,out)
    evaluator.clear()
    return dict(pairs=pairs,capacity=capacity,tree=tree,Q=b.peak_qubits,T=b._cur_stats['toffoli']+E.extra,peak_cache=evaluator.peak,products=evaluator.products,evictions=evaluator.evictions)


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,nargs='+',default=[3,4]);p.add_argument('--capacity',type=int,nargs='+',default=[6,8,10]);p.add_argument('--tree',action='store_true');args=p.parse_args();rows=[]
    for pairs in args.pairs:
        for capacity in args.capacity:
            row=profile(pairs,capacity,args.tree);rows.append(row);print(json.dumps(row),flush=True)
            Path('results/field_'+('tree' if args.tree else 'dag')+'_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
