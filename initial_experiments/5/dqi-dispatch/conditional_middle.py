"""Coherent control of trimmed series arithmetic with one field-word cache.

Each alias-free bilinear product masks its left evaluated operand, emits the
ordinary field product, then erases the mask. Linear Frobenius terms become
controlled linear terms. This avoids adding a quantum control to measurement
instructions or materializing multiple heads simultaneously.
"""
from core_components import TrimmedArithmetic

class _LinearEmitter:
    def __init__(self,E,flag):self.base=E;self.flag=flag
    def __getattr__(self,name):return getattr(self.base,name)
    def emit_mac(self,rt):
        assert all(t[0]=='L' for t in rt)
        self.base.emit_mac([('P',(self.flag,),0,t[1],t[2],t[3]) for t in rt])

class ConditionalArithmetic(TrimmedArithmetic):
    def __init__(self,m,E,flag):
        super().__init__(m,_LinearEmitter(E,flag));self.base=E;self.flag=flag;self.known_zero.add(())
    def product(self,left,right,targets):
        m=self.m;E=self.base
        with m.opened(E,'conditional_evaluation',12) as masked:
            rt=[('P',(self.flag,),0,tuple(bits),0,(masked[i],)) for i,bits in enumerate(left) if bits]
            E.mac(rt);out=[];m.bil_terms(out,m.KT,m.spec_plain(masked),right,tws_of=targets);E.emit_mac(out);E.free(rt)

def mobius(m,E,math,A,D,p,T,flag,inverse=False):
    """Conditional ordinary Möbius head, p is the stored reciprocal scale."""
    def scale(back):
        with m.opened(E,'conditional_diagonal',12) as factor:
            E.x(factor[0]);rt=[('P',(flag,),0,(a,),0,(b,)) for a,b in zip(p,factor)]
            E.mac(rt);E.cx(flag,factor[0])
            with m.inv_word(E,factor) as pi,m.opened(E,'conditional_scale_product',12) as temp:
                for word in T:
                    if word:m.mon_slot(E,word,pi if back else factor,factor if back else pi,temp)
            E.cx(flag,factor[0]);E.x(factor[0]);E.free(rt)
    def reciprocal():math.inverse([A[0],(),*T])
    if inverse:
        math.triangular(A,1,T,True);math.triangular(D,1,T);reciprocal();scale(True);math.triangular(D,1,T,True)
    else:
        math.triangular(D,1,T);scale(False);reciprocal();math.triangular(D,1,T,True);math.triangular(A,1,T)
