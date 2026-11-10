"""Fixed EEA with one skipped degree and a quadratic quotient on the next step."""
import kickmix as km
from lean2_dqi import Lean2Decoder
from toom_dqi import plan, lin_xor


class RepairDecoder(Lean2Decoder):
    def _clear_expression(self,b,tmp,linear,products):
        meas=b.alloc_bits(12)
        b.hmr(tmp,meas)
        if linear is not None:
            for j in range(12):
                b.push_condition(meas[j]);b.z(linear[j]);b.pop_condition()
        for lhs,rhs in products:
            b.gf2_phase_by_product(lhs,rhs,mask=meas,field=self.gf)
        b.bit_store(meas,False);b.free(meas);b.free(tmp)

    def _relation(self,b,X,Y,i,active,pending,q0,q1,q2,which,compute):
        T=2*self.t
        # q1 holds the highest quotient term, q0 the next, q2 the optional
        # third term until the final conditional cyclic permutation.
        if not compute:
            b.cswap(pending,q1,q0);b.cswap(pending,q2,q1)
        if which=='r':
            top,nxt,third=X[T+1-i],X[T-i],X[T+2-i]
            b.cswap(pending,top,third);b.cswap(pending,nxt,third)
            lead,ynext,y2=Y[T-i],Y[T-i-1],Y[T-i-2]
            yregs=None
        else:
            top,nxt,third=X[T-i],X[T-i+1],X[T-i+2]
            yregs=[Y[T-i+1+j] for j in range(min(4,i))]
            for j in range(len(yregs)-1):b.cswap(pending,yregs[j],yregs[j+1])
            lead=yregs[0]
            ynext=yregs[1] if len(yregs)>1 else None
            y2=yregs[2] if len(yregs)>2 else None
        g=self._gate(b,lead,active)

        def highest():
            if compute:self._mul(b,top,g,q1)
            else:b.del_gf2_mul(top,g,target=q1,field=self.gf)

        def middle():
            tmp=b.alloc_qubits(12);b.cx(nxt,tmp)
            products=[]
            if ynext is not None:
                self._mul(b,q1,ynext,tmp);products.append((q1,ynext))
            if compute:self._mul(b,tmp,g,q0)
            else:b.del_gf2_mul(tmp,g,target=q0,field=self.gf)
            self._clear_expression(b,tmp,nxt,products)

        def lowest():
            tmp=b.alloc_qubits(12);b.cx(third,tmp)
            products=[]
            for lhs,rhs in [(q0,ynext),(q1,y2)]:
                if rhs is not None:
                    self._mul(b,lhs,rhs,tmp);products.append((lhs,rhs))
            prod=b.alloc_qubits(12);self._mul(b,tmp,g,prod)
            if compute:b.init_and(pending,prod,q2)
            else:b.del_and(pending,prod,q2)
            b.del_gf2_mul(tmp,g,target=prod,field=self.gf);b.free(prod)
            self._clear_expression(b,tmp,third,products)

        for fn in ([highest,middle,lowest] if compute else [lowest,middle,highest]):fn()
        self._ungate(b,lead,active,g)
        if which=='r':
            b.cswap(pending,nxt,third);b.cswap(pending,top,third)
        else:
            for j in reversed(range(len(yregs)-1)):b.cswap(pending,yregs[j],yregs[j+1])
        if compute:
            b.cswap(pending,q2,q1);b.cswap(pending,q1,q0)

    def _poly_mul2_xor(self,b,q0,q1,q2,ys,targets):
        n=len(ys)
        pts,E,out=plan(n+1)
        p=b.alloc_qubits(12)
        for j,a in enumerate(pts):
            ev=[(E[j][k],ys[k],ys[0]) for k in range(1,n)]
            qev=[(a,q1,q0),(self.gf.mul(a,a),q2,q0)]
            lin_xor(b,ev);lin_xor(b,qev)
            self._mul(b,q0,ys[0],p)
            lin_xor(b,[(out[j][k],p,targets[k]) for k in range(n+2)])
            b.del_gf2_mul(q0,ys[0],target=p,field=self.gf)
            lin_xor(b,qev);lin_xor(b,ev)
        self._mul(b,q2,ys[-1],p)
        lin_xor(b,[(out[n+1][k],p,targets[k]) for k in range(n+2)])
        b.del_gf2_mul(q2,ys[-1],target=p,field=self.gf);b.free(p)

    def _step2(self,b,X,Y,i,active,pending,forward):
        if i==1:
            return super()._step(b,X,Y,i,active,forward)
        T=2*self.t
        q0=b.alloc_qubits(12);q1=b.alloc_qubits(12);q2=b.alloc_qubits(12)
        self._relation(b,X,Y,i,active,pending,q0,q1,q2,'r' if forward else 's',True)
        self._poly_mul2_xor(b,q0,q1,q2,[Y[k] for k in range(T-i+1)],
                            [X[k] for k in range(T-i+3)])
        self._poly_mul2_xor(b,q0,q1,q2,[Y[T-d] for d in range(i)],
                            [X[T-d] for d in range(i+2)])
        self._relation(b,X,Y,i,active,pending,q0,q1,q2,'s' if forward else 'r',False)
        b.free(q0);b.free(q1);b.free(q2)

    def _bad(self,b,Y,i,c,index,bad):
        T=2*self.t
        work=b.alloc_qubits(6);e,z,nz,iz,a0,a1=work
        b.flip_if_equal(c,i-1,target=e)
        b.flip_if_equal(Y[T-i],0,target=z)
        b.flip_if_equal(Y[T-i-1],0,target=nz);b.x(nz)
        b.flip_if_equal(index,0,target=iz)
        b.init_and(e,z,a0);b.init_and(a0,nz,a1);b.ccx(a1,iz,bad)
        b.del_and(a0,nz,a1);b.del_and(e,z,a0)
        b.flip_if_equal(index,0,target=iz)
        b.x(nz);b.flip_if_equal(Y[T-i-1],0,target=nz)
        b.flip_if_equal(Y[T-i],0,target=z)
        b.flip_if_equal(c,i-1,target=e);b.free(work)

    def _write_index(self,b,index,i,bad):
        b.cx(bad,km.array([index[j] for j in range(len(index)) if i>>j&1]))

    def _eea2(self,b,A,Bv,c,index,forward):
        for i in (range(1,self.t+1) if forward else range(self.t,0,-1)):
            X,Y=(A,Bv) if i%2 else (Bv,A)
            work=b.alloc_qubits(3);active,bad,pending=work
            if forward:
                self._xor_activity(b,Y,i,c,active)
                self._bad(b,Y,i,c,index,bad)
                self._write_index(b,index,i,bad)
                if i>=2:b.flip_if_equal(index,i-1,target=pending)
                self._step2(b,X,Y,i,active,pending,True)
                for x,y in zip(A,Bv):b.cswap(bad,x,y)
                b.cx(bad,active);b.iadd(1,target=c,control=active)
                b.flip_if_equal(c,i,target=active)
                if i>=2:b.flip_if_equal(index,i-1,target=pending)
                b.flip_if_equal(index,i,target=bad)
            else:
                b.flip_if_equal(index,i,target=bad)
                if i>=2:b.flip_if_equal(index,i-1,target=pending)
                b.flip_if_equal(c,i,target=active)
                b.isub(1,target=c,control=active);b.cx(bad,active)
                for x,y in zip(A,Bv):b.cswap(bad,x,y)
                self._step2(b,X,Y,i,active,pending,False)
                self._xor_activity(b,Y,i,c,active)
                self._write_index(b,index,i,bad)
                self._bad(b,Y,i,c,index,bad)
                if i>=2:b.flip_if_equal(index,i-1,target=pending)
            b.free(work)

    def build(self,b,syndromes,error_mbuc,chien=True):
        T=2*self.t
        A=[b.alloc_qubits(12) for _ in range(T+1)]
        Bv=list(syndromes)+[b.alloc_qubits(12)]
        c=b.alloc_qubits(6);index=b.alloc_qubits(6)
        b.x(A[T][0]);b.x(Bv[T][0])
        self._eea2(b,A,Bv,c,index,True)
        for x,y in zip(A,Bv):b.cswap(c[0],x,y)
        if chien:self._chien(b,Bv,error_mbuc)
        for x,y in zip(A,Bv):b.cswap(c[0],x,y)
        self._eea2(b,A,Bv,c,index,False)
        b.x(Bv[T][0]);b.x(A[T][0])
        b.free(c);b.free(index);b.free(Bv[T])
        for r in A:b.free(r)


def build_circuit(m=4095,n=70,chien=True):
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder())
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(n)]
    mb=[b.create_classical_register(12,f'E{i}') for i in range(m)]
    RepairDecoder(n,km.GF2Field(12)).build(b,syn,mb,chien)
    return b.finish_circuit()


if __name__=='__main__':
    import sys,time
    start=time.time();c=build_circuit(n=int(sys.argv[2]) if len(sys.argv)>2 else 70)
    print('built',c.num_qubits,c.operation_statistics().num_toffoli_like,time.time()-start,flush=True)
    c.to_file(sys.argv[1],format='kmb')
