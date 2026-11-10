"""Exact GF(2) row basis and nullspace utilities, using Python integers."""
def basis(vectors):
    pivots={}
    for x in vectors:
        while x:
            p=x.bit_length()-1
            if p in pivots:x^=pivots[p]
            else:pivots[p]=x;break
    return [pivots[p] for p in sorted(pivots,reverse=True)]

def coordinates(vectors,bs):
    out=[]
    for x in vectors:
        code=0
        for j,v in enumerate(bs):
            if x>>(v.bit_length()-1)&1:x^=v;code^=1<<j
        assert x==0;out.append(code)
    return out

def nullspace(vectors,width):
    bs=basis(vectors);pivots={x.bit_length()-1 for x in bs};out=[]
    for bit in range(width):
        if bit in pivots:continue
        v=1<<bit
        for row in reversed(bs):
            p=row.bit_length()-1
            if (row&v).bit_count()&1:v^=1<<p
        assert all((row&v).bit_count()%2==0 for row in bs)
        out.append(v)
    return out

def solve_many(vectors,values,width,outputs):
    pivots={}
    for x,y in zip(vectors,values):
        while x:
            p=x.bit_length()-1
            if p in pivots:
                xx,yy=pivots[p];x^=xx;y^=yy
            else:pivots[p]=(x,y);break
        if not x:assert not y,'inconsistent interpolation'
    solutions=[0]*outputs
    for p,(row,rhs) in sorted(pivots.items()):
        for j in range(outputs):
            if ((row&solutions[j]).bit_count()^(rhs>>j))&1:solutions[j]^=1<<p
    assert all(sum(((x&s).bit_count()%2)<<j for j,s in enumerate(solutions))==y for x,y in zip(vectors,values))
    return solutions

def row_reduction(matrix):
    rows=list(matrix);n=len(rows);ops=[]
    for bit in range(n):
        p=next(i for i in range(bit,n) if rows[i]>>bit&1)
        if p!=bit:rows[p],rows[bit]=rows[bit],rows[p];ops.append(('swap',p,bit))
        for i in range(n):
            if i!=bit and rows[i]>>bit&1:rows[i]^=rows[bit];ops.append(('cx',i,bit))
    assert rows==[1<<i for i in range(n)]
    return ops
