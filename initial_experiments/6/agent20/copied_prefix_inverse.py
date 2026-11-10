"""Unit-series inversion by retaining one reciprocal prefix.

Let k=ceil(n/2), A=F mod x^k, U=A^-1 mod x^k. Then F*U^2=F^-1
mod x^n, and its first k coefficients equal U. Copy A, invert the copy,
multiply F in place by U(x)^2, and erase U against the new prefix. The
constant coefficient remains implicit one and its physical word is ignored.
A width budget selects this method at compile time, with the prior exact
parity recurrence as fallback. All input/phase behavior stays reversible.
"""
from contextlib import ExitStack
from parity_series_inverse import ParityArithmetic


class CopiedArithmetic(ParityArithmetic):
    def __init__(self,*args,limit=1004,recursive=False,**kwargs):
        super().__init__(*args,**kwargs);self.limit=limit;self.recursive=recursive
        self.copied_calls=[];self._copy_depth=0

    def inverse(self,F):
        n=len(F);k=(n+1)//2
        needed=sum(12 for w in F[1:k] if tuple(w) not in self.known_zero)
        if n<4 or (self._copy_depth and not self.recursive) or self.E.b.active_qubits+needed+4>self.limit:
            return super().inverse(F)
        m,E=self.m,self.E;self.copied_calls.append((n,needed,E.b.active_qubits))
        with ExitStack() as stack:
            had_empty=() in self.known_zero;self.known_zero.add(())
            U=[()]
            for j,word in enumerate(F[1:k],1):
                if tuple(word) in self.known_zero:U.append(())
                else:
                    q=stack.enter_context(m.opened(E,f'copied_unit_prefix_{n}_{j}',12));U.append(q)
                    for a,b in zip(word,q):E.cx(a,b)
            self._copy_depth+=1
            try:self.inverse(U)
            finally:self._copy_depth-=1
            # Multiply the nonconstant target first. The omitted F0=1
            # contribution is an affine Frobenius-square addition afterward.
            self.triangular(U,2,[(),*F[1:]])
            for j,word in enumerate(U[1:],1):
                if word and 2*j<n:
                    rt=[];m.lin_rows_terms(rt,self.rows(1,2),word,F[2*j]);E.mac(rt)
            if not had_empty:self.known_zero.remove(())
            for old,new in zip(F[1:k],U[1:]):
                for a,b in zip(old,new):E.cx(a,b)
