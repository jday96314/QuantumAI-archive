"""Repair one degree-two Euclid quotient using a persistent boundary offset."""
import kickmix as km
from repair_dqi import RepairDecoder
from toom_dqi import plan,lin_xor


class OffsetDecoder(RepairDecoder):
    def _relation_offset(self,b,X,Y,i,active,offset,bad,q0,q1,q2,which,compute):
        T=2*self.t
        if not compute:
            b.cswap(bad,q1,q0);b.cswap(bad,q2,q1)
        selected=b.alloc_qubits(1)[0];b.cx(offset,selected);b.cx(bad,selected)
        extras=[]
        def reg(array,k):
            if 0<=k<=T:return array[k]
            z=b.alloc_qubits(12);extras.append(z);return z
        if which=='r':
            xr=[reg(X,T+1-i-j) for j in range(4)]
            yr=[reg(Y,T-i-j) for j in range(4)]
            for j in range(3):b.cswap(offset,xr[j],xr[j+1])
            for j in range(3):b.cswap(selected,yr[j],yr[j+1])
            top,nxt,third=xr[:3];lead,ynext,y2=yr[:3]
        else:
            xr=[reg(X,T-i-1+j) for j in range(4)]
            yr=[reg(Y,T-i+j) for j in range(min(4,i+1))]
            for j in reversed(range(3)):b.cswap(selected,xr[j],xr[j+1])
            for j in reversed(range(len(yr)-1)):b.cswap(offset,yr[j],yr[j+1])
            top,nxt,third=xr[1:4]
            lead=yr[1];ynext=yr[2] if len(yr)>2 else None;y2=yr[3] if len(yr)>3 else None
        g=self._gate(b,lead,active)

        def highest():
            if compute:self._mul(b,top,g,q1)
            else:b.del_gf2_mul(top,g,target=q1,field=self.gf)

        def middle():
            tmp=b.alloc_qubits(12);b.cx(nxt,tmp);products=[]
            if ynext is not None:
                self._mul(b,q1,ynext,tmp);products.append((q1,ynext))
            if compute:self._mul(b,tmp,g,q0)
            else:b.del_gf2_mul(tmp,g,target=q0,field=self.gf)
            self._clear_expression(b,tmp,nxt,products)

        def lowest():
            tmp=b.alloc_qubits(12);b.cx(third,tmp);products=[]
            for lhs,rhs in [(q0,ynext),(q1,y2)]:
                if rhs is not None:
                    self._mul(b,lhs,rhs,tmp);products.append((lhs,rhs))
            prod=b.alloc_qubits(12);self._mul(b,tmp,g,prod)
            if compute:b.init_and(bad,prod,q2)
            else:b.del_and(bad,prod,q2)
            b.del_gf2_mul(tmp,g,target=prod,field=self.gf);b.free(prod)
            self._clear_expression(b,tmp,third,products)

        for fn in ([highest,middle,lowest] if compute else [lowest,middle,highest]):fn()
        self._ungate(b,lead,active,g)
        if which=='r':
            for j in reversed(range(3)):b.cswap(selected,yr[j],yr[j+1])
            for j in reversed(range(3)):b.cswap(offset,xr[j],xr[j+1])
        else:
            for j in range(len(yr)-1):b.cswap(offset,yr[j],yr[j+1])
            for j in range(3):b.cswap(selected,xr[j],xr[j+1])
        for r in extras:b.free(r)
        b.cx(bad,selected);b.cx(offset,selected);b.free(selected)
        if compute:
            b.cswap(bad,q2,q1);b.cswap(bad,q1,q0)

    def _poly_bounded(self,b,q0,q1,q2,ys,targets):
        # The two potentially highest terms vanish on both supported layouts.
        D=len(targets)-1
        pts,E,out=plan(D+1)
        p=b.alloc_qubits(12)
        for j,a in enumerate(pts):
            ev=[(E[j][k],ys[k],ys[0]) for k in range(1,len(ys))]
            qev=[(a,q1,q0),(self.gf.mul(a,a),q2,q0)]
            lin_xor(b,ev);lin_xor(b,qev)
            self._mul(b,q0,ys[0],p)
            lin_xor(b,[(out[j][k],p,targets[k]) for k in range(D+1)])
            b.del_gf2_mul(q0,ys[0],target=p,field=self.gf)
            lin_xor(b,qev);lin_xor(b,ev)
        b.free(p)

    def _step_offset(self,b,X,Y,i,active,offset,bad,forward):
        T=2*self.t
        q0=b.alloc_qubits(12);q1=b.alloc_qubits(12);q2=b.alloc_qubits(12)
        self._relation_offset(b,X,Y,i,active,offset,bad,q0,q1,q2,'r' if forward else 's',True)
        extra=b.alloc_qubits(12);b.cswap(offset,Y[T-i],extra)
        self._poly_bounded(b,q0,q1,q2,[Y[k] for k in range(T-i+1)],
                           [X[k] for k in range(T-i+2)])
        self._poly_bounded(b,q0,q1,q2,[Y[T-d] for d in range(i)]+[extra],
                           [X[T-d] for d in range(i+2)])
        b.cswap(offset,Y[T-i],extra);b.free(extra)
        self._relation_offset(b,X,Y,i,active,offset,bad,q0,q1,q2,'s' if forward else 'r',False)
        b.free(q0);b.free(q1);b.free(q2)

    def _activity_offset(self,b,Y,i,c,offset,bad,active):
        T=2*self.t
        b.cswap(offset,Y[T-i],Y[T-i-1])
        if i==self.t:
            tmp=b.alloc_qubits(1)[0]
            self._xor_activity(b,Y,i,c,tmp)
            b.x(offset);b.ccx(offset,tmp,active);b.x(offset)
            self._xor_activity(b,Y,i,c,tmp);b.free(tmp)
        else:
            self._xor_activity(b,Y,i,c,active);b.cx(bad,active)
        b.cswap(offset,Y[T-i],Y[T-i-1])

    def _eea_offset(self,b,A,Bv,c,index,offset,forward):
        for i in (range(1,self.t+1) if forward else range(self.t,0,-1)):
            X,Y=(A,Bv) if i%2 else (Bv,A)
            work=b.alloc_qubits(2);active,bad=work
            if forward:
                if i<self.t:self._bad(b,Y,i,c,index,bad)
                self._activity_offset(b,Y,i,c,offset,bad,active)
                self._write_index(b,index,i,bad)
                self._step_offset(b,X,Y,i,active,offset,bad,True)
                b.cx(bad,offset)
                b.iadd(1,target=c,control=active);b.flip_if_equal(c,i,target=active)
                b.flip_if_equal(index,i,target=bad)
            else:
                b.flip_if_equal(index,i,target=bad);b.cx(bad,offset)
                b.flip_if_equal(c,i,target=active);b.isub(1,target=c,control=active)
                self._step_offset(b,X,Y,i,active,offset,bad,False)
                self._activity_offset(b,Y,i,c,offset,bad,active)
                self._write_index(b,index,i,bad)
                if i<self.t:self._bad(b,Y,i,c,index,bad)
            b.free(work)

    def build(self,b,syndromes,error_mbuc,chien=True):
        T=2*self.t
        A=[b.alloc_qubits(12) for _ in range(T+1)]
        Bv=list(syndromes)+[b.alloc_qubits(12)]
        c=b.alloc_qubits(6);index=b.alloc_qubits(6);offset=b.alloc_qubits(1)[0]
        b.x(A[T][0]);b.x(Bv[T][0])
        self._eea_offset(b,A,Bv,c,index,offset,True)
        for x,y in zip(A,Bv):b.cswap(c[0],x,y)
        if chien:self._chien(b,Bv,error_mbuc)
        for x,y in zip(A,Bv):b.cswap(c[0],x,y)
        self._eea_offset(b,A,Bv,c,index,offset,False)
        b.x(Bv[T][0]);b.x(A[T][0])
        b.free(c);b.free(index);b.free(offset);b.free(Bv[T])
        for r in A:b.free(r)


def build_circuit(m=4095,n=70,chien=True):
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder())
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(n)]
    mb=[b.create_classical_register(12,f'E{i}') for i in range(m)]
    OffsetDecoder(n,km.GF2Field(12)).build(b,syn,mb,chien)
    return b.finish_circuit()


if __name__=='__main__':
    import sys,time
    start=time.time();c=build_circuit(n=int(sys.argv[2]) if len(sys.argv)>2 else 70)
    print('built',c.num_qubits,c.operation_statistics().num_toffoli_like,time.time()-start,flush=True)
    c.to_file(sys.argv[1],format='kmb')
