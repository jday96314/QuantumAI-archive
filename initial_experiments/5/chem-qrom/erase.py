"""Streaming ANF phase repair with a care-prefix polynomial completion.

Only addresses < len(table) are promised. The ANF coefficient transform is
triangular in numeric order; coefficients beyond that prefix may be set to
zero. This gives an exact extension on every cared-for input and prunes the
monomial banks, including monomials introduced solely by zero padding.
"""
import kickmix as km


def coefficients(table, n, care=True):
    count=len(table);values=list(map(int,table))+[0]*((1<<n)-count)
    for i in range(n):
        for mask in range(1<<n):
            if mask>>i&1:values[mask]^=values[mask^(1<<i)]
    if care:values[count:]=[0]*((1<<n)-count)
    return values


def closure(masks):
    out=set(masks)|{0}
    pending=list(out)
    while pending:
        m=pending.pop();bits=m
        while bits:
            bit=bits&-bits;bits^=bit;rest=m^bit
            if rest not in out:out.add(rest);pending.append(rest)
    return out


def plan(table,n,split=None,care=True):
    coeff=coefficients(table,n,care)
    choices=[]
    for s in ([split] if split is not None else range(1,n)):
        low=closure(m&((1<<s)-1) for m,c in enumerate(coeff) if c)
        high=closure(m>>s for m,c in enumerate(coeff) if c)
        nl=sum(m.bit_count()>1 for m in low);nh=sum(m.bit_count()>1 for m in high)
        workspace=nl+max((m.bit_count()-1 for m in high),default=0)
        choices.append((nl+nh,workspace,s,low,high))
    return min(choices,key=lambda x:x[:3]),coeff


def erase(b,address,target,table,split=None,care=True):
    (gates,workspace,s,lmasks,hmasks),coeff=plan(table,len(address),split,care)
    assert workspace<=len(target),(workspace,len(target))
    outcomes=b.alloc_bits(len(target));b.hmr(target,outcomes)
    low={0:True};records=[];cursor=0
    for m in sorted(lmasks-{0}):
        bit=m&-m;rest=m^bit;q=address[bit.bit_length()-1]
        if not rest:low[m]=q
        else:
            out=target[cursor];cursor+=1;b.init_and(low[rest],q,out)
            low[m]=out;records.append((low[rest],q,out))
    flag=b.alloc_bits(1)[0]
    def apply(mask,high):
        for lm,lq in low.items():
            c=coeff[lm|(mask<<s)]
            if c:
                selected=km.array([outcomes[j] for j in range(len(outcomes)) if c>>j&1])
                b.cx(selected,flag);b.push_condition(flag);b.cz(lq,high);b.pop_condition();b.cx(selected,flag)
    apply(0,True)
    def visit(mask,start,parent,depth):
        for j in range(start,len(address)-s):
            mm=mask|(1<<j)
            if mm not in hmasks:continue
            q=address[s+j]
            if not mask:
                apply(mm,q);visit(mm,j+1,q,depth)
            else:
                out=target[cursor+depth];b.init_and(parent,q,out)
                apply(mm,out);visit(mm,j+1,out,depth+1);b.del_and(parent,q,out)
    visit(0,0,True,0)
    for a,c,q in reversed(records):b.del_and(a,c,q)
    b.bit_store(outcomes,False);b.free(outcomes);b.free(km.array([flag]))
    return dict(T=gates,workspace=workspace,split=s,coefficients=sum(bool(c) for c in coeff))
