"""Low-workspace Möbius map T -> A*T/(p+z^2*D*T) modulo z^n.

A,D have constant one, p is nonzero. Conjugating a single unit-series inverse
avoids inversions of T[0], which can be zero. This research circuit uses the
attributed upstream field/gate generator; it is not a submission decoder.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from series_inverse import SeriesArithmetic
from block_euclid_ratio import series_divide,series_product

ROOT=Path(__file__).resolve().parent


class BoundedArithmetic(SeriesArithmetic):
    def __init__(self,*args):
        super().__init__(*args);self.known_zero=set()
    def spec(self,terms):
        return super().spec([(c,w,q) for c,w,q in terms if tuple(w) not in self.known_zero])
    def sparse_middle(self,A,q,low,high):
        m=self.m;n=len(low)+len(high)-1;pts,I=m.vinv(n)
        pairs=[(i,j,(len(low)+i-j)//q) for i in range(len(high)) for j in range(len(low))
               if len(low)+i-j>=q and (len(low)+i-j)%q==0 and
               (len(low)+i-j)//q<len(A) and tuple(A[(len(low)+i-j)//q]) not in self.known_zero]
        if len(pairs)<=n:
            for i,j,h in pairs:
                self.product(self.spec([(1,A[h],q)]),self.spec([(1,low[j],1)]),lambda tb,i=i:[high[i][b] for b in tb])
            return
        for a,point in enumerate(pts):
            left=self.spec([(I[j-1][a],A[j//q],q) for j in range(q,n+1,q) if j//q<len(A)])
            if not any(left):continue
            right=self.spec([(m.gpow(point,j),word,1) for j,word in enumerate(reversed(low))])
            powers=[m.gpow(point,i) for i in range(len(high))]
            def targets(tb):
                return [high[i][k] for i,power in enumerate(powers) for k in m.bitlist(m.gmul(power,m.tmask(tb)),12)]
            self.product(left,right,targets)

    def triangular_suffix(self,A,q,target,start,inverse=False):
        """Multiply only the suffix by a unit series, preserving its prefix.

        This is the lower triangular block [[I,0],[M,A_tail]]. Its inverse
        also leaves the prefix fixed; the recursion uses the original low
        coefficients before each forward cross update and after its inverse.
        """
        if len(target)<=1 or start>=len(target):return
        if start<=0:
            self.triangular(A,q,target,inverse);return
        cut=len(target)//2;lo,hi=target[:cut],target[cut:]
        def middle():
            skip=max(0,start-cut)
            if not skip:self.sparse_middle(A,q,lo,hi);return
            m=self.m;n=len(lo)+len(hi)-1
            pairs=[(i,j,(len(lo)+i-j)//q)
                   for i in range(skip,len(hi)) for j in range(len(lo))
                   if len(lo)+i-j>=q and (len(lo)+i-j)%q==0
                   and (len(lo)+i-j)//q<len(A)
                   and tuple(A[(len(lo)+i-j)//q]) not in self.known_zero]
            if len(pairs)<=n:
                for i,j,h in pairs:
                    self.product(self.spec([(1,A[h],q)]),self.spec([(1,lo[j],1)]),lambda tb,i=i:[hi[i][b] for b in tb])
                return
            pts,I=m.vinv(n)
            for a,point in enumerate(pts):
                left=self.spec([(I[j-1][a],A[j//q],q) for j in range(q,n+1,q) if j//q<len(A)])
                if not any(left):continue
                right=self.spec([(m.gpow(point,j),word,1) for j,word in enumerate(reversed(lo))])
                powers=[m.gpow(point,i) for i in range(skip,len(hi))]
                def targets(tb):
                    return [hi[i][k] for i,power in zip(range(skip,len(hi)),powers)
                            for k in m.bitlist(m.gmul(power,m.tmask(tb)),12)]
                self.product(left,right,targets)
        if inverse:
            self.triangular_suffix(A,q,lo,start,True);middle()
            self.triangular_suffix(A,q,hi,max(0,start-cut),True)
        else:
            self.triangular_suffix(A,q,hi,max(0,start-cut));middle()
            self.triangular_suffix(A,q,lo,start)


def apply(m,E,math,A,D,p,T,inverse=False,p_is_inverse=False):
    def scale(inv):
        with m.inv_word(E,p) as iv,m.opened(E,'mobius_product',12) as product:
            scale_inverse=bool(inv)^bool(p_is_inverse)
            for word in T:m.mon_slot(E,word,p if scale_inverse else iv,iv if scale_inverse else p,product)
    def reciprocal_shift():
        # Coefficient0 is implicitly one and never accessed by the reciprocal
        # emitter. Coefficient1 remains zero throughout that computation.
        if getattr(math,'virtual_zero',False):
            zero=()
            math.known_zero.add(tuple(zero))
            math.inverse([A[0],zero,*T])
            math.known_zero.remove(tuple(zero))
        else:
            with m.opened(E,'mobius_zero',12) as zero:
                math.known_zero.add(tuple(zero))
                math.inverse([A[0],zero,*T])
                math.known_zero.remove(tuple(zero))
    if inverse:
        math.triangular(A,1,T,True)
        math.triangular(D,1,T)
        reciprocal_shift()
        scale(True)
        math.triangular(D,1,T,True)
    else:
        math.triangular(D,1,T)
        scale(False)
        reciprocal_shift()
        math.triangular(D,1,T,True)
        math.triangular(A,1,T)


def build(n,k,which):
    m=load();b=m.Builder()
    A=[b.qubits(12) for _ in range(k+1)];D=[b.qubits(12) for _ in range(k)]
    p=b.qubits(12);T=[b.qubits(12) for _ in range(n)]
    sA,sB,sT=b.qubits(3);alloc=m.Alloc(b);E=m.Em(b,sA,sB,sT,alloc);math=BoundedArithmetic(m,E)
    if which in ['forward','both']:apply(m,E,math,A,D,p,T)
    if which in ['inverse','both']:apply(m,E,math,A,D,p,T,True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'mobius_n{n}_k{k}_{which}'
    tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for name,words in [('A',A),('D',D),('p',[p]),('T',T)]:
        for i in range(len(words)):kb.create_quantum_register(12,f'{name}{i}')
    pre=ROOT/f'parameters/{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters/mobius.empty';empty.write_text('');out=ROOT/f'candidates/{tag}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb')
    return c


def test(n,k,shots=512):
    rng=np.random.default_rng(948521);aa=rng.integers(0,4096,(k+1,shots));dd=rng.integers(0,4096,(k,shots))
    aa[0]=1;dd[0]=1;pp=rng.integers(1,4096,(1,shots));tt=rng.integers(0,4096,(n,shots));tt[:,0]=0;tt[0,1::7]=0
    expected=[]
    for shot in range(shots):
        t=tt[:,shot].tolist();a=aa[:,shot].tolist();d=dd[:,shot].tolist();p=int(pp[0,shot])
        expected.append(series_divide(series_product(t,a,n),[p,0]+series_product(t,d,n-2),n))
    expected=np.array(expected).T;results=[];began=time.time()
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for which in ['forward','inverse','both']:
        c=build(n,k,which);sim=km.Simulator(batch_size=shots,seed=7192);sim.use_same_registers_as(c)
        inputs={'A':aa,'D':dd,'p':pp,'T':expected if which=='inverse' else tt}
        for name,values in inputs.items():
            for i,row in enumerate(values):sim.write_across_shots(c.register_data[f'{name}{i}'],bits(row))
        sim.do(c)
        for name,values in inputs.items():
            want=(expected if which=='forward' else tt) if name=='T' else values
            for i,row in enumerate(want):
                assert np.array_equal(sim.read_across_shots(c.register_data[f'{name}{i}']),bits(row)),(which,name,i)
        public={q for reg in c.register_data.values() for q in reg}
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        results.append(dict(which=which,Q=c.num_qubits,public_qubits=len(public),extra_qubits=c.num_qubits-len(public),T=c.operation_statistics().num_toffoli_like))
    out=dict(status='PASS',scope='Möbius series map, inverse and composition, including zero T constant; monic A/D and nonzero p',
             n=n,k=k,shots=shots,circuits=results,seconds=time.time()-began)
    (ROOT/f'results/mobius_n{n}_k{k}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=16);p.add_argument('--k',type=int,default=8);p.add_argument('--shots',type=int,default=512)
    a=p.parse_args();test(a.n,a.k,a.shots)
