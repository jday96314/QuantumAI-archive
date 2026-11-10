"""Explicit reversible polynomial transforms; no free 'invertible map'.

The monic modulus has quantum low coefficients. Schoolbook convolution and
long division are used as a working control, not claimed optimal. Every XOR
MAC computes a full product, divides, copies a remainder, reverses division,
and erases the original product. Primitive GF12 multiplication is the same
42-Toffoli attributed CRT implementation as the selected round-4 decoder.

The dry emitter traverses exactly the same allocation/gate calls as native
emission; degree-2/3 native cross-checks test arbitrary targets and moduli.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from tower_inverse import product as field_product, linear_xor
from oracle import mul, power, product, divide, derivative, inverse_mod, known_error_polynomials


class Meter:
    def __init__(self):
        self.live=set();self.free_ids=[];self.frontier=0;self.peak=0;self.T=0;self.events=[];self.stage='init'
    def alloc_qubits(self,n):
        out=[]
        for _ in range(n):
            if self.free_ids:q=self.free_ids.pop()
            else:q=self.frontier;self.frontier+=1
            assert q not in self.live
            self.live.add(q);out.append(km.q(q))
        self.peak=max(self.peak,len(self.live));self.events.append([self.stage,'alloc',n,len(self.live),self.T])
        return km.array(out)
    def free(self,qs):
        for q in qs:self.live.remove(q.id);self.free_ids.append(q.id)
        self.events.append([self.stage,'free',len(qs),len(self.live),self.T])
    def cx(self,*args):pass
    def init_and(self,*args):self.T+=1
    def del_and(self,*args):pass
    def swap(self,*args):pass


def words(b,n):return [b.alloc_qubits(12) for _ in range(n)]
def drop(b,ws):
    for w in reversed(ws):b.free(w)


def fmultiply(b,a,c,out):
    if isinstance(b,Meter):
        scratch=b.alloc_qubits(3);b.T+=42;b.free(scratch)
    else:field_product(b,a,c,out)


def modular_mac(b, a, c, sigma, out, square=False):
    """out ^= (a*c mod monic sigma); all inputs and workspace restored.

    sigma contains low d coefficients; its x^d coefficient is implicit one.
    out may contain arbitrary data. a,c,out,sigma are disjoint (a=c allowed).
    """
    d=len(sigma);assert len(a)==len(c)==len(out)==d
    work=words(b,2*d-1)
    columns=[mul(1<<j,1<<j) for j in range(12)]
    def convolution():
        if square:
            for i in range(d):linear_xor(b,a[i],work[2*i],columns)
        else:
            for i in range(d):
                for j in range(d):fmultiply(b,a[i],c[j],work[i+j])
    def reduce(inverse=False):
        for k in (range(d,2*d-1) if inverse else reversed(range(d,2*d-1))):
            # Keep high coefficients as quotient data; each elimination is
            # an explicit reversible multiply-add, including the inverse.
            for j in (reversed(range(d)) if inverse else range(d)):
                fmultiply(b,work[k],sigma[j],work[k-d+j])
    convolution();reduce()
    for i in range(d):b.cx(work[i],out[i])
    reduce(True);convolution();drop(b,work)


def derivative_words(b,sigma):
    d=len(sigma);out=words(b,d)
    for i in range(1,d,2):b.cx(sigma[i],out[i-1])
    if d%2:
        # The implicit leading coefficient one supplies derivative x^(d-1).
        b.x(out[d-1][0])
    return out


def f_overwrite(b,omega,psi,sigma,inverse=False,work=None):
    """F=omega*Psi mod sigma via copy/swap/inverse multiplication.

    Retain work until the paired inverse, even if sigma is nonsquarefree or
    Psi is invalid. Then arbitrary inputs still restore exactly. At the
    midpoint work=omega_old + F*sigma_prime; correct paths have work=0.
    It is NEVER treated as an unconditionally free/reusable clean register.
    """
    own=work is None
    if own:work=words(b,len(omega))
    deriv=derivative_words(b,sigma)
    if inverse:
        modular_mac(b,omega,deriv,sigma,work)
        for a,c in zip(omega,work):b.swap(a,c)
        modular_mac(b,omega,psi,sigma,work)
    else:
        modular_mac(b,omega,psi,sigma,work)
        for a,c in zip(omega,work):b.swap(a,c)
        modular_mac(b,omega,deriv,sigma,work)
    # Derivative is only a Clifford copy; erase it explicitly.
    for i in range(1,len(sigma),2):b.cx(sigma[i],deriv[i-1])
    if len(sigma)%2:b.x(deriv[-1][0])
    drop(b,deriv)
    return work


def inversion_xor(b,a,sigma,out):
    """out ^= a^4094 mod sigma, total reversible even for nonsplit sigma.

    This gives pointwise inverse (zero->zero) only for split squarefree sigma.
    All 21 intermediate polynomials remain live for compute/copy/uncompute.
    This is a deliberately explicit expensive alternative to half-degree EEA.
    """
    live=[];ops=[];v=a;d=len(a)
    for _ in range(10):
        sq=words(b,d);modular_mac(b,v,v,sigma,sq,True);live.append(sq);ops.append((v,v,sq,True))
        nxt=words(b,d);modular_mac(b,sq,a,sigma,nxt);live.append(nxt);ops.append((sq,a,nxt,False));v=nxt
    sq=words(b,d);modular_mac(b,v,v,sigma,sq,True);live.append(sq);ops.append((v,v,sq,True))
    for x,y in zip(sq,out):b.cx(x,y)
    for x,y,z,is_square in reversed(ops):modular_mac(b,x,y,sigma,z,is_square);drop(b,live.pop())


def native_test(d):
    rng=np.random.default_rng(530125+d);reports=[];n=1024
    bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    def run(label, emit, expected, inputs):
        b=FastBuilder(km.CircuitBuilder());regs={name:[b.create_quantum_register(12,f'{name}{i}') for i in range(d)] for name in inputs}
        emit(b,regs);c=b.finish_circuit();sim=km.Simulator(batch_size=n,seed=991+d);sim.use_same_registers_as(c)
        for name,arr in inputs.items():
            for reg,v in zip(regs[name],arr):sim.write_across_shots(reg,bits(v))
        sim.do(c)
        for name,arr in expected.items():
            for reg,v in zip(regs[name],arr):assert np.array_equal(sim.read_across_shots(reg),bits(v)),(label,d,name)
        public=set(q for ws in regs.values() for reg in ws for q in reg)
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(label,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(label,'phase')
        c.to_file(f'candidates/{label}_d{d}.kmb',format='kmb')
        reports.append(dict(kind=label,degree=d,shots=n,status='PASS',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    inputs={k:rng.integers(0,4096,(d,n)) for k in ['a','c','sigma','out']}
    expected={k:v.copy() for k,v in inputs.items()}
    for s in range(n):
        r=divide(product(inputs['a'][:,s],inputs['c'][:,s]),[*inputs['sigma'][:,s],1])[1]
        for i,x in enumerate(r):expected['out'][i,s]^=x
    run('modular_mac',lambda b,r:modular_mac(b,r['a'],r['c'],r['sigma'],r['out']),expected,inputs)
    inputs={k:rng.integers(0,4096,(d,n)) for k in ['omega','psi','sigma']}
    def roundtrip(b,r):
        w=f_overwrite(b,r['omega'],r['psi'],r['sigma']);f_overwrite(b,r['omega'],r['psi'],r['sigma'],True,w);drop(b,w)
    run('F_roundtrip',roundtrip,inputs,inputs)
    inputs={k:np.zeros((d,n),dtype=np.int64) for k in ['omega','psi','sigma','out']}
    expected={k:v.copy() for k,v in inputs.items()}
    for s in range(n):
        pos=rng.choice(4095,d,replace=False).tolist();vals=rng.integers(1,4096,d).tolist()
        sigma,omega,_=known_error_polynomials(pos,vals)
        # Convert the independent oracle's constant-normalized locator to
        # the monic polynomial used by the compact decoder.
        sigma=[mul(x,power(sigma[-1],4094)) for x in sigma]
        psi=inverse_mod(derivative(sigma),sigma)
        ff=divide(product(omega,psi),sigma)[1]
        for name,poly in [('omega',omega),('sigma',sigma[:-1]),('psi',psi)]:
            inputs[name][:len(poly),s]=poly;expected[name][:len(poly),s]=poly
        expected['omega'][:,s]=0;expected['omega'][:len(ff),s]=ff
    run('F_forward',lambda b,r:f_overwrite(b,r['omega'],r['psi'],r['sigma'],work=r['out']),expected,inputs)
    inputs={k:rng.integers(0,4096,(d,n)) for k in ['a','sigma','out']}
    expected={k:v.copy() for k,v in inputs.items()}
    for s in range(n):
        sigma=[*inputs['sigma'][:,s],1];a=list(inputs['a'][:,s]);v=a
        for _ in range(10):v=divide(product(divide(product(v,v),sigma)[1],a),sigma)[1]
        v=divide(product(v,v),sigma)[1]
        for i,x in enumerate(v):expected['out'][i,s]^=x
    run('inverse_xor',lambda b,r:inversion_xor(b,r['a'],r['sigma'],r['out']),expected,inputs)
    return reports


def budgets():
    reports=[]
    for d in [1,2,3,17,35]:
        for kind in ['multiply','square','F_overwrite_lifetime','inverse_xor']:
            b=Meter();sigma=words(b,d);a=words(b,d);c=words(b,d);out=words(b,d)
            baseline=48*d
            if kind=='multiply':modular_mac(b,a,c,sigma,out)
            elif kind=='square':modular_mac(b,a,a,sigma,out,True)
            elif kind=='inverse_xor':inversion_xor(b,a,sigma,out)
            else:
                # Include a full persistent 35-word compact prior bank when
                # translating these component widths to a full decoder.
                # a=omega, c=Psi; out is the retained old-omega work bank.
                # d extra derivative words live only during either pass.
                # The meter accepts X as a free operation below.
                f_overwrite(b,a,c,sigma,work=out);f_overwrite(b,a,c,sigma,True,out)
            assert len(b.live)==baseline
            reports.append(dict(kind=kind,degree=d,T=b.T,Q=b.peak,public_Q=baseline,scratch_Q=b.peak-baseline,
                                scope='exact dry-emitter count, schoolbook control, native cross-check at degrees 2/3'))
    return reports


Meter.x=lambda self,*args:None


if __name__=='__main__':
    checks=native_test(2)+native_test(3)
    dry=budgets()
    for check in checks:
        kind={'modular_mac':'multiply','F_roundtrip':'F_overwrite_lifetime','F_forward':'F_overwrite_lifetime','inverse_xor':'inverse_xor'}[check['kind']]
        row=next(r for r in dry if r['kind']==kind and r['degree']==check['degree'])
        assert check['T']==row['T']//(2 if check['kind']=='F_forward' else 1),(check,row)
        # F test has three public arrays, whereas the dry meter retains its
        # fourth array explicitly as work; total peak is identical.
        assert check['Q']==row['Q']-(12*check['degree'] if check['kind']=='inverse_xor' else 0),(check,row)
    result=dict(status='PASS',native_checks=checks,exact_emitter_budgets=dry,
                exclusions=['no compact-prior -> omega adapter included','no degree selector/padding adapter included','no native whole-decoder cost claim'])
    Path('results/modular_schedule.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
