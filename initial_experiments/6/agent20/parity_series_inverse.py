"""In-place reciprocal via even/odd decomposition in characteristic two.

For f(x)=E(x^2)+x O(x^2), put R=O/E. Then
  G=E^-1/(1+y R^2), H=R G, 1/f(x)=G(x^2)+x H(x^2).
Thus O/=E; E=inverse(E); E/=(1+y O^2); O*=E.
All four steps are reversible in the original words. The constant E0 is
implicitly one, as in the attributed read-only SeriesArithmetic interface.
"""
from cube_arithmetic import CubeArithmetic


class ParityArithmetic(CubeArithmetic):
    def shifted_square_middle(self,A,low,high):
        m=self.m;cut=len(low)
        pairs=[(i,j,(cut+i-j-1)//2) for i in range(len(high)) for j in range(len(low))
               if cut+i-j>=1 and (cut+i-j)%2==1
               and (cut+i-j-1)//2<len(A)
               and tuple(A[(cut+i-j-1)//2]) not in self.known_zero]
        linear=[(i,k) for i,j,k in pairs if tuple(low[j])==self.constant_word]
        pairs=[(i,j,k) for i,j,k in pairs if tuple(low[j])!=self.constant_word]
        for i,k in linear:
            rt=[];m.lin_rows_terms(rt,self.rows(1,2),A[k],high[i]);self.E.emit_mac(rt)
        # Both an ordinary rank-n middle product and direct sparse products
        # are valid. Evaluation/interpolation has free binary linear maps.
        n=len(low)+len(high)-1
        if len(pairs)<=n:
            for i,j,k in pairs:
                self.product(self.spec([(1,A[k],2)]),self.spec([(1,low[j],1)]),
                             lambda tb,i=i:[high[i][b] for b in tb])
        elif pairs:
            points,I=m.vinv(n)
            for a,point in enumerate(points):
                left=self.spec([(I[2*k][a],word,2) for k,word in enumerate(A) if 2*k<n])
                right=self.spec([(m.gpow(point,j),word,1) for j,word in enumerate(reversed(low))
                                 if tuple(word)!=self.constant_word])
                def target(tb):
                    return [high[i][b] for i in range(len(high))
                            for b in m.bitlist(m.gmul(m.gpow(point,i),m.tmask(tb)),12)]
                self.product(left,right,target)

    def shifted_square(self,A,target,inverse=False):
        if len(target)<=1:return
        cut=len(target)//2;lo,hi=target[:cut],target[cut:]
        if inverse:
            self.shifted_square(A,lo,True)
            self.shifted_square_middle(A,lo,hi)
            self.shifted_square(A,hi,True)
        else:
            self.shifted_square(A,hi)
            self.shifted_square_middle(A,lo,hi)
            self.shifted_square(A,lo)

    def inverse(self,F):
        if len(F)<=1:return
        self.constant_word=tuple(F[0])
        even,odd=F[::2],F[1::2]
        self.triangular(even,1,odd,True)
        self.inverse(even)
        self.shifted_square(odd,even,True)
        self.triangular(even,1,odd)

