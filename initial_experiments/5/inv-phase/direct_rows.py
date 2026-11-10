"""Phase signed, independently rounded Booth rows without forming their sum."""
import kickmix as km
from short_adder import add
from components import shift,phase

def apply(b,source,coordinate,flags,exponent,denom=37):
    L=len(coordinate)
    assert L%2==0
    for g in range(L//2):
        cur=coordinate[2*g];prev=coordinate[2*g-1] if g else None
        sign=coordinate[2*g+1];diff,two=b.alloc_qubits(2)
        if prev is not None:b.cx(prev,cur)
        b.cx(sign,diff);b.cx(cur,diff)
        if prev is not None:b.cx(prev,diff)
        b.x(cur);b.init_and(cur,diff,two);b.x(cur)
        claims={};round_claims=[]
        for branch,ctrl in enumerate((cur,two)):
            e=exponent+2*g+branch
            for i,q in enumerate(source):
                if i+e>=0:claims.setdefault(i+e,[]).append((ctrl,q))
            rs=-e-1
            if e<0 and rs<len(source):round_claims.append((ctrl,source[rs]))
        core=max(claims,default=-1)+1+bool(round_claims)
        y=b.alloc_qubits(core+10)
        for j,terms in claims.items():
            for a,c in terms:b.ccx(a,c,y[j])
        carry=b.alloc_qubits(int(bool(round_claims)))
        if round_claims:
            for a,c in round_claims:b.ccx(a,c,carry[0])
            add(b,carry,y[:core],carry_workspace=y[core:])
        shift(b,y,core,flags);phase(b,y,denom,sign,coefficient=-1);shift(b,y,core,flags,True)
        if round_claims:add(b,carry,y[:core],subtract=True,carry_workspace=y[core:])
        def erase(q,terms):
            outcome=b.alloc_bits(1)[0];b.hmr(q,outcome);b.push_condition(outcome)
            for a,c in terms:b.cz(a,c)
            b.pop_condition();b.bit_store(outcome,False);b.free(km.array([outcome]))
        for j,terms in claims.items():erase(y[j],terms)
        if round_claims:erase(carry[0],round_claims)
        b.free(carry);b.free(y)
        b.x(cur);b.del_and(cur,diff,two);b.x(cur)
        if prev is not None:b.cx(prev,diff)
        b.cx(cur,diff);b.cx(sign,diff)
        if prev is not None:b.cx(prev,cur)
        b.free(km.array([diff,two]))
