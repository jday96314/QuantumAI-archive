"""Native census probes, importing Apache-2.0 field emitters read-only.

No public verifier is changed or invoked. Outputs remain in this directory.
C probes the disjoint-register rational update; A probes Frobenius reduction
and an explicitly charged reciprocal-product involution. These are components.
"""
import argparse, hashlib, json, subprocess, time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from mobius_series import BoundedArithmetic
from published_dqi import isomorphism
from oracle import product, divide, derivative, inverse_mod, mul, power, add

ROOT=Path(__file__).resolve().parent
TRANSLATOR=ROOT/'../../4/codex/translate_prefix'

def rational(m,E,math,sigma,tau,alpha,V,shift=1):
    # sigma[0]=1. This is an involution on (sigma,V) at fixed tau,alpha.
    # All coefficient arrays are disjoint. Sharing sigma with V is INVALID.
    math.triangular(sigma,1,V,True)
    for i in range(shift,len(sigma)):
        if i-shift<len(tau):E.emit_mac(m.prod_rt(sigma[i],alpha,tau[i-shift]))
    math.triangular(sigma,1,V)


def modular_mac(m,E,math,a,c,sigma,out,square=False):
    """out ^= a*c mod (x^d+sigma). Quantum modulus; arbitrary target."""
    d=len(sigma)
    with m.opened(E,'census_convolution',12*(2*d-1)) as flat:
        w=[flat[12*i:12*i+12] for i in range(2*d-1)]
        if square:
            rt=[]
            rows=math.rows(1,2)
            for i,word in enumerate(a):m.lin_rows_terms(rt,rows,word,w[2*i])
        else:rt=m.polyprod_rt(a,c,[[word] for word in w])
        high=w[d:]; coeff=[()]+list(reversed(sigma))
        def reduction(inverse=False):
            if inverse:
                if high:E.emit_mac(m.polyprod_rt(sigma,high,[[w[i]] if i<d else [] for i in range(2*d-2)]))
                math.triangular(coeff,1,list(reversed(high)))
            else:
                math.triangular(coeff,1,list(reversed(high)),True)
                if high:E.emit_mac(m.polyprod_rt(sigma,high,[[w[i]] if i<d else [] for i in range(2*d-2)]))
        E.emit_mac(rt);reduction()
        for x,y in zip(w,out):
            for xx,yy in zip(x,y):E.cx(xx,yy)
        reduction(True);E.emit_mac(rt)


def inverse_xor(m,E,math,a,sigma,out):
    """a^4094 in split quotient algebra, zero totalized at EACH root.

    Every intermediate stays live. No assumption that Frobenius modulo a
    quantum modulus is Clifford. Each square calls charged reduction.
    """
    d=len(a);ops=[];held=[];v=a
    for step in range(10):
        flat=E.A.open(f'census_sq_{step}',12*d)
        sq=[flat[12*i:12*i+12] for i in range(d)];held.append((f'census_sq_{step}',flat))
        modular_mac(m,E,math,v,v,sigma,sq,True);ops.append((v,v,sq,True))
        flat=E.A.open(f'census_mul_{step}',12*d)
        nxt=[flat[12*i:12*i+12] for i in range(d)];held.append((f'census_mul_{step}',flat))
        modular_mac(m,E,math,sq,a,sigma,nxt);ops.append((sq,a,nxt,False));v=nxt
    flat=E.A.open('census_sq_final',12*d)
    sq=[flat[12*i:12*i+12] for i in range(d)];held.append(('census_sq_final',flat))
    modular_mac(m,E,math,v,v,sigma,sq,True);ops.append((v,v,sq,True))
    for x,y in zip(sq,out):
        for xx,yy in zip(x,y):E.cx(xx,yy)
    for (x,y,z,square),(name,flat) in zip(reversed(ops),reversed(held)):
        modular_mac(m,E,math,x,y,sigma,z,square);E.A.close(name,flat)


def reciprocal_product_xor(m,E,math,a,sigma,out):
    d=len(a)
    with m.opened(E,'census_derivative',12*d) as flat:
        df=[flat[12*i:12*i+12] for i in range(d)]
        def copy_df():
            for i in range(1,d,2):
                for x,y in zip(sigma[i],df[i-1]):E.cx(x,y)
            if d%2:E.x(df[-1][0])
        copy_df()
        with m.opened(E,'census_abprime',12*d) as cflat:
            c=[cflat[12*i:12*i+12] for i in range(d)]
            modular_mac(m,E,math,a,df,sigma,c)
            inverse_xor(m,E,math,c,sigma,out)
            modular_mac(m,E,math,a,df,sigma,c)
        copy_df()


def fprime(m,E,math,a,sigma,work):
    reciprocal_product_xor(m,E,math,a,sigma,work)
    E.swaps(None,m.word_pairs(a,work))
    reciprocal_product_xor(m,E,math,a,sigma,work)


def quadratic_xor(m,E,math,a,sigma,out):
    """Degree-2 norm formula, using ONE scalar inverse, no modular powers.

    B=x²+b1*x+b0, norm(A)=a0²+a0*a1*b1+a1²*b0.
    f(A)=(a0+a1*b1+a1*x)/(b1*norm(A)). Nonunits are totalized
    to zero here; the retained-work involution still restores arbitrary data.
    """
    assert len(a)==len(sigma)==len(out)==2
    from low_workspace_inverse import bind
    inv_word=bind(m)
    with m.opened(E,'quadratic_u',12) as u,m.opened(E,'quadratic_norm',12) as norm,m.opened(E,'quadratic_den',12) as den:
        urt=m.prod_rt(u,a[1],sigma[1]);E.emit_mac(urt)
        def norm_terms(erase=False):
            with m.opened(E,'quadratic_square',12) as sq:
                rt=[];m.lin_rows_terms(rt,math.rows(1,2),a[1],sq);E.emit_mac(rt)
                nr=[];m.lin_rows_terms(nr,math.rows(1,2),a[0],norm)
                nr+=m.prod_rt(norm,a[0],u)+m.prod_rt(norm,sq,sigma[0])
                (E.emit_free if erase else E.emit_mac)(nr)
                E.emit_mac(rt)
        norm_terms()
        dr=m.prod_rt(den,sigma[1],norm);E.emit_mac(dr)
        with inv_word(E,den) as iv:
            E.emit_mac(m.prod_rt(out[0],a[0],iv)+m.prod_rt(out[0],u,iv)+m.prod_rt(out[1],a[1],iv))
        E.emit_free(dr);norm_terms(True);E.emit_free(urt)


def quadratic(m,E,math,a,sigma,work):
    quadratic_xor(m,E,math,a,sigma,work)
    E.swaps(None,m.word_pairs(a,work))
    quadratic_xor(m,E,math,a,sigma,work)


def quadratic_low_xor(m,E,math,a,sigma,out):
    """Same norm formula, release/recompute norm around scalar inversion.

    Squared a1 enters as a parity specification, with no materialized word.
    The additional retained out bank is still charged as public component data.
    """
    assert len(a)==len(sigma)==len(out)==2
    from low_workspace_inverse import bind
    with m.opened(E,'quad_low_u',12) as u,m.opened(E,'quad_low_den',12) as den:
        ur=m.prod_rt(u,a[1],sigma[1]);E.emit_mac(ur)
        def with_norm(callback):
            with m.opened(E,'quad_low_norm',12) as norm:
                nr=[];m.lin_rows_terms(nr,math.rows(1,2),a[0],norm)
                nr+=m.prod_rt(norm,a[0],u)
                m.bil_terms(nr,m.KT,math.spec([(1,a[1],2)]),m.spec_plain(sigma[0]),norm)
                E.emit_mac(nr);callback(norm);E.emit_free(nr)
        with_norm(lambda norm:E.emit_mac(m.prod_rt(den,sigma[1],norm)))
        with bind(m)(E,den) as iv:
            E.emit_mac(m.prod_rt(out[0],a[0],iv)+m.prod_rt(out[0],u,iv)+m.prod_rt(out[1],a[1],iv))
        with_norm(lambda norm:E.emit_free(m.prod_rt(den,sigma[1],norm)))
        E.emit_free(ur)


def quadratic_low(m,E,math,a,sigma,work):
    quadratic_low_xor(m,E,math,a,sigma,work)
    E.swaps(None,m.word_pairs(a,work))
    quadratic_low_xor(m,E,math,a,sigma,work)


def make(kind,n=8,d=3,repeat=1,count=False):
    m=load();b=m.Builder()
    if kind=='rational':sizes={'sigma':d+1,'tau':d,'alpha':1,'V':n}
    else:sizes={'a':d,'sigma':d,'out':d}
    regs={name:[b.qubits(12) for _ in range(size)] for name,size in sizes.items()}
    class Counter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=(Counter if count else m.Em)(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E)
    math.known_zero.add(())
    for _ in range(repeat):
        if kind=='rational':rational(m,E,math,regs['sigma'],regs['tau'],regs['alpha'][0],regs['V'])
        elif kind=='square':modular_mac(m,E,math,regs['a'],regs['a'],regs['sigma'],regs['out'],True)
        elif kind=='fprime':fprime(m,E,math,regs['a'],regs['sigma'],regs['out'])
        elif kind=='quadratic':quadratic(m,E,math,regs['a'],regs['sigma'],regs['out'])
        elif kind=='quadratic_low':quadratic_low(m,E,math,regs['a'],regs['sigma'],regs['out'])
        else:raise ValueError(kind)
    T=(E.T if count else 0)+b._stream_totals()[0]
    return m,b,regs,dict(kind=kind,n=n,d=d,repeat=repeat,Q=b.peak_qubits,T=T,public_Q=12*sum(sizes.values()),scratch_Q=b.peak_qubits-12*sum(sizes.values()))


def export(b,regs,tag):
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tape=ROOT/'parameters'/f'{tag}.npy';np.save(tape,rows)
    pb=km.CircuitBuilder()
    for name,words in regs.items():
        for i in range(len(words)):pb.create_quantum_register(12,f'{name}{i}')
    pre=ROOT/'parameters'/f'{tag}.pre.kmx';pb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters/empty';empty.write_text('');path=ROOT/'candidates'/f'{tag}.kmx'
    subprocess.run([str(TRANSLATOR),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(path),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    circ=km.Circuit.from_file(str(path),format='kmx');circ.to_file(str(path.with_suffix('.kmb')),format='kmb');path.unlink()
    return circ


def series_div(a,b,n):
    out=[0]*n
    assert b[0]==1
    for i in range(n):
        out[i]=int(a[i]) if i<len(a) else 0
        for j in range(1,min(i+1,len(b))):out[i]^=mul(b[j],out[i-j])
    return out


def padded(a,n):return np.asarray(list(a)+[0]*(n-len(a)),dtype=np.int64)


def check(kind,n,d,shots,repeat=1,seed=2026102426):
    start=time.time();m,b,regs,info=make(kind,n,d,repeat)
    tag=f'{kind}_n{n}_d{d}_r{repeat}';circ=export(b,regs,tag)
    rng=np.random.default_rng(seed+d+n)
    inputs={name:rng.integers(0,4096,(len(words),shots)) for name,words in regs.items()}
    if kind=='rational':
        inputs['sigma'][0]=1
        inputs['alpha'][:,0]=0;inputs['tau'][:,1]=0;inputs['V'][:,2]=0
    if kind in ['fprime','quadratic','quadratic_low'] and repeat%2:
        # Split, squarefree B and unit A: native midpoint map contract.
        for s in range(shots):
            roots=rng.choice(np.arange(1,4096),d,replace=False).tolist();B=[1]
            for root in roots:B=product(B,[root,1])
            inputs['sigma'][:,s]=B[:-1]
            # A constant nonzero guarantees a unit; half shots general units.
            if s%2==0:
                inputs['a'][:,s]=0;inputs['a'][0,s]=int(rng.integers(1,4096))
            else:
                while True:
                    try:inverse_mod(inputs['a'][:,s].tolist(),B);break
                    except ValueError:inputs['a'][:,s]=rng.integers(0,4096,d)
        inputs['out'][:]=0
    if kind in ['fprime','quadratic','quadratic_low'] and repeat%2==0:
        inputs['sigma'][:,0]=0;inputs['a'][:,0]=0
        inputs['sigma'][:,1]=0;inputs['sigma'][0,1]=1
        inputs['a'][:,2]=0;inputs['out'][:,3]=0
        if d==2 and shots>=4096:inputs['sigma'][1,:4096]=np.arange(4096)
    want={k:v.copy() for k,v in inputs.items()}
    for s in range(shots):
        if kind=='rational' and repeat%2:
            sig=inputs['sigma'][:,s].tolist();tau=inputs['tau'][:,s].tolist();alpha=int(inputs['alpha'][0,s])
            S=series_div(inputs['V'][:,s],sig,n);new=sig[:]
            for i in range(1,d+1):new[i]^=mul(alpha,tau[i-1])
            want['sigma'][:,s]=new;want['V'][:,s]=padded(product(S,new)[:n],n)
        elif kind=='square' and repeat%2:
            rem=divide(product(inputs['a'][:,s],inputs['a'][:,s]),[*inputs['sigma'][:,s],1])[1]
            want['out'][:,s]^=padded(rem,d)
        elif kind in ['fprime','quadratic','quadratic_low'] and repeat%2:
            B=[*inputs['sigma'][:,s],1]
            f=inverse_mod(divide(product(inputs['a'][:,s],derivative(B)),B)[1],B)
            want['a'][:,s]=padded(f,d)
    fw,bk,transform=isomorphism(m)
    bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(circ)
    for name,arr in inputs.items():
        for i,v in enumerate(arr):sim.write_across_shots(circ.register_data[f'{name}{i}'],bits(np.array([transform(int(x),fw) for x in v])))
    sim.do(circ)
    mismatches=0
    for name,arr in want.items():
        for i,v in enumerate(arr):
            got=sim.read_across_shots(circ.register_data[f'{name}{i}'])
            mismatches+=int(np.count_nonzero(got^bits(np.array([transform(int(x),fw) for x in v]))))
    public={q for reg in circ.register_data.values() for q in reg}
    scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
    scratch_bad=int(np.count_nonzero(sim.read_across_shots(scratch)))
    phase_bad=int(np.count_nonzero(sim.read_phase_flipped_across_shots()))
    profile=make(kind,n,d,repeat,True)[3]
    assert profile['Q']==circ.num_qubits,(profile,circ.num_qubits)
    assert profile['T']==circ.operation_statistics().num_toffoli_like,(profile,circ.operation_statistics().num_toffoli_like)
    result=dict(**info,status='PASS' if not(mismatches or scratch_bad or phase_bad) else 'FAIL',shots=shots,seed=seed,mismatching_bits=mismatches,scratch_bits=scratch_bad,phase_errors=phase_bad,seconds=time.time()-start,kmb=str((ROOT/'candidates'/f'{tag}.kmb').relative_to(ROOT)),sha256=hashlib.sha256((ROOT/'candidates'/f'{tag}.kmb').read_bytes()).hexdigest(),scope='native component; field converted from independent host-basis reference; source values, targets, scratch and measurement phase checked; not full decoder verification')
    (ROOT/'results'/f'{tag}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    assert result['status']=='PASS',result
    return result


def profiles():
    result=[]
    for kind,n,d in [('rational',8,3),('rational',22,11),('rational',34,17),('rational',70,35),('square',0,3),('square',0,17),('square',0,35),('fprime',0,3),('fprime',0,35)]:
        info=make(kind,n,d,1,True)[3];result.append(info);print(json.dumps(info),flush=True)
    (ROOT/'results/profiles.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('kind',choices=['rational','square','fprime','quadratic','quadratic_low','profiles']);p.add_argument('--n',type=int,default=8);p.add_argument('--d',type=int,default=3);p.add_argument('--shots',type=int,default=512);p.add_argument('--repeat',type=int,default=1)
    a=p.parse_args()
    if a.kind=='profiles':profiles()
    else:check(a.kind,a.n,a.d,a.shots,a.repeat)
