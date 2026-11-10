"""One-gap H/Sylvester repair with a retained field word and source-only tags.

At a missing degree-one pivot, r has degree n-2. Put k=floor(padding/2),
c=gamma*x^k+x*c1, retaining the promised k leading payload zeros. Replace
r by x*r, recurse with c1, and add gamma*x^k*r. One retained gamma word
supports the first such gap. Every step has an explicit inverse, including
unsupported source/flag states. No full-decoder acceptance is claimed here.
"""
from contextlib import contextmanager
import fast_sylvester as fast
from compact_monomial_unpad import degree_count,rotate


class Repair:
    def __init__(self,m,E,completion):
        self.m=m;self.E=E;self.completion=completion;self.full=len(completion)
        self.tags=E.A.open('H_gap_tags',max(2,(self.full//2).bit_length()))
        self.gamma=E.A.open('H_gap_gamma',12)
        self.count=E.A.open('H_gap_padding_count',self.full.bit_length())
        degree_count(E.b,completion,self.count)

    def close(self):
        E=self.E;degree_count(E.b,self.completion,self.count)
        E.A.close('H_gap_padding_count',self.count)
        E.A.close('H_gap_gamma',self.gamma);E.A.close('H_gap_tags',self.tags)

    @contextmanager
    def selected(self,depth):
        m,E=self.m,self.E
        with m.opened(E,'H_gap_selected',1) as word:
            flag=word[0]
            with m.eq_test(E,self.tags,depth+1) as bit:E.cx(bit,flag)
            yield flag
            with m.eq_test(E,self.tags,depth+1) as bit:E.cx(bit,flag)

    def toggle_tag(self,r,depth,inverse=False):
        m,E=self.m,self.E
        with m.eq_test(E,r[-1],0) as zero,m.eq_test(E,r[-2],0) as nextzero:
            E.x(nextzero)
            with m.opened(E,'H_gap_eligible',1) as elig,m.opened(E,'H_gap_tag_write',1) as select:
                rt=m.one_rt((zero,),0,(nextzero,),0,elig[0]);E.mac(rt)
                def eligibility():
                    with m.eq_test(E,self.tags,0) as empty:E.xor(m.one_rt((elig[0],),0,(empty,),0,select[0]))
                def identity():
                    with m.eq_test(E,self.tags,depth+1) as equal:E.cx(equal,select[0])
                if inverse:identity()
                else:eligibility()
                for j,tag in enumerate(self.tags):
                    if (depth+1)>>j&1:E.cx(select[0],tag)
                if inverse:eligibility()
                else:identity()
                E.free(rt)
            E.x(nextzero)

    @contextmanager
    def k_equal(self,k):
        m,E=self.m,self.E
        with m.opened(E,'H_gap_k_equal',1) as word:
            bit=word[0]
            indices=(self.full-2*k,self.full-2*k-2)
            for index in indices:
                if index==self.full:E.x(bit)
                elif index>=0:E.cx(self.completion[index],bit)
            yield bit
            for index in reversed(indices):
                if index==self.full:E.x(bit)
                elif index>=0:E.cx(self.completion[index],bit)

    def rotation(self,bank,control,amount,inverse=False):
        E=self.E
        class WordBuilder:
            def cswap(self,ctl,a,b):E.swaps(ctl,list(zip(a,b)))
        rotate(WordBuilder(),bank,control,amount,inverse)

    def payload_shift(self,c,depth,inverse=False,flag=None):
        m,E=self.m,self.E;n=len(c)
        assert flag is not None
        if inverse:E.swaps(flag,list(zip(c[-1],self.gamma)))
        if not inverse:self.rotation(c,flag,1)
        for k in (reversed(range(1,n-1)) if inverse else range(1,n-1)):
            with self.k_equal(k) as equal,m.opened(E,'H_gap_selected_k',1) as select:
                rt=m.one_rt((flag,),0,(equal,),0,select[0]);E.mac(rt)
                E.swaps(select[0],list(zip(c[k-1],c[-1])));E.free(rt)
        if inverse:self.rotation(c,flag,1,True)
        if not inverse:E.swaps(flag,list(zip(c[-1],self.gamma)))

    def correction(self,r,low,depth,inverse=False,flag=None):
        m,E=self.m,self.E;n=len(r)
        assert flag is not None
        # Rotation is only needed for k<=n-2 on supported repaired charts.
        controls=self.count[1:1+(n-2).bit_length()]
        def orient(back=False):
            for j in (reversed(range(len(controls))) if back else range(len(controls))):
                self.rotation(r,controls[j],1<<j,not back)
        with m.opened(E,'H_gap_masked_gamma',12) as gamma,m.opened(E,'H_gap_correction_product',12) as product:
            rt=[term for x,y in zip(self.gamma,gamma) for term in m.one_rt((flag,),0,(x,),0,y)]
            E.mac(rt);orient()
            for i in (reversed(range(n)) if inverse else range(n)):
                pr=m.prod_rt(product,gamma,r[i]);E.mac(pr)
                for x,y in zip(product,low[i]):E.cx(x,y)
                index=self.full-2*i-2
                if i+n<len(low) and index>=0:
                    wrap=self.completion[index]
                    for p,a,b in zip(product,low[i],low[i+n]):
                        E.cx(a,b);E.b.ccx(wrap,p,a);E.cx(a,b)
                E.free(pr)
            orient(True);E.free(rt)

    def apply(self,m,E,f,g,c,d,inverse=False,g_one=False,depth=0):
        mm,n=len(f)-1,len(g)-1
        assert mm>=n>=0 and len(c)==n and len(d)==mm
        with m.opened(E,'repaired_sylv_zero',1) as zw:
            zero=zw[0]
            with m.eq_test(E,g[-1],0) as z:E.cx(z,zero)
            E.cx(zero,g[-1][0])
            if not n:fast.syl.scale(m,E,d,g[0],inverse,g_one)
            else:
                k=n-1;q=f[n:];r=f[:n]
                low=fast.output_layout(n,k,d[:k],c);high=d[k:]
                def shear():E.xor(m.polyprod_rt(q,c,[[w] for w in d]))
                def cross():
                    targets=[[low[i+k]] if i<n else [] for i in range(len(g)+len(high)-1)]
                    E.xor(m.polyprod_rt(g,high,targets))
                fast.divide(m,E,f,g,one=g_one)
                if inverse:
                    fast.syl.triangular(m,E,g,high,True,g_one);cross()
                    if n>=2:
                        with self.selected(depth) as flag:
                            self.correction(r,low,depth,True,flag)
                            self.rotation(r,flag,1,True)
                    self.apply(m,E,g,r,d[:k],c,True,depth=depth+1)
                    if n>=2:
                        with self.selected(depth) as flag:
                            self.rotation(r,flag,1);self.payload_shift(c,depth,True,flag)
                        self.toggle_tag(r,depth,True)
                    shear()
                else:
                    shear()
                    if n>=2:
                        self.toggle_tag(r,depth)
                        with self.selected(depth) as flag:
                            self.payload_shift(c,depth,flag=flag);self.rotation(r,flag,1,True)
                    self.apply(m,E,g,r,d[:k],c,depth=depth+1)
                    if n>=2:
                        with self.selected(depth) as flag:
                            self.rotation(r,flag,1);self.correction(r,low,depth,flag=flag)
                    cross();fast.syl.triangular(m,E,g,high,False,g_one)
                fast.divide(m,E,f,g,True,g_one)
            E.cx(zero,g[-1][0])
            with m.eq_test(E,g[-1],0) as z:E.cx(z,zero)
        return fast.output_layout(mm,n,c,d)
