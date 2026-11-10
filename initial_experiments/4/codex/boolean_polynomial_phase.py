"""Combinational Boolean model of the selected inverse-square-root phase.

Research front end for phase-bit emission and measurement-based pebbling.
This does not build an accepted quantum circuit. It checks the arithmetic
against the independently implemented integer model before scheduling gates.
"""
import json
from pathlib import Path
import numpy as np
from boolean_lookup import Network,Signal,ZERO,ONE
import affine_invsqrt as affine
import binary_centered_invsqrt as poly


class Logic:
    def __init__(self,n):self.net=Network(n)
    def add(self,a,b,carry=ZERO):
        result=[]
        for x,y in zip(a,b):
            result.append(x^y^carry);carry=x^self.net.mul(x^y,x^carry)
        return result
    def resize(self,a,n):return (list(a)+[ZERO]*n)[:n]
    def rounded(self,coef,e,n):
        if e>=0:return self.resize([ZERO]*e+list(coef),n)
        r=coef[-e-1] if -e<=len(coef) else ZERO
        value=self.resize(coef[-e:],n)
        return self.add(value,[ZERO]*n,r) if r!=ZERO else value
    def horner(self,acc,coef,z,delta):
        n=len(acc);mul=self.net.mul
        for g in range((len(z)+1)//2):
            bit=lambda k:ZERO if k<0 else z[min(k,len(z)-1)]
            low,mid,high=bit(2*g-1),bit(2*g),bit(2*g+1)
            mag1=low^mid;mag2=mul(~mag1,mid^high)
            v1=self.rounded(coef,delta-21+2*g,n)
            v2=self.rounded(coef,delta-21+2*g+1,n)
            negative=mul(mag1^mag2,~high)
            row=[mul(mag1,a)^mul(mag2,b)^negative for a,b in zip(v1,v2)]
            acc=self.add(acc,row,negative)
        return acc


class CarrySaveLogic(Logic):
    def sum_rows(self,rows,n):
        columns=[[] for _ in range(n)]
        for row in rows:
            for i,value in enumerate(row):
                if value!=ZERO:columns[i].append(value)
        for i in range(n):
            while len(columns[i])>2:
                x,y,z=[columns[i].pop() for _ in range(3)]
                columns[i].append(x^y^z)
                if i+1<n:columns[i+1].append(x^self.net.mul(x^y,x^z))
        first=[c[0] if c else ZERO for c in columns]
        second=[c[1] if len(c)>1 else ZERO for c in columns]
        return self.add(first,second)
    def horner(self,acc,coef,z,delta):
        n=len(acc);mul=self.net.mul;rows=[acc]
        for g in range((len(z)+1)//2):
            bit=lambda k:ZERO if k<0 else z[min(k,len(z)-1)]
            low,mid,high=bit(2*g-1),bit(2*g),bit(2*g+1)
            mag1=low^mid;mag2=mul(~mag1,mid^high)
            e=delta-21+2*g
            def floor_round(e):
                if e>=0:return self.resize([ZERO]*e+list(coef),n),ZERO
                return self.resize(coef[-e:],n),coef[-e-1] if -e<=len(coef) else ZERO
            v1,r1=floor_round(e);v2,r2=floor_round(e+1)
            negative=mul(mag1^mag2,~high)
            rows.append([mul(mag1,a)^mul(mag2,b)^negative for a,b in zip(v1,v2)])
            rows.append([negative]+[ZERO]*(n-1))
            rounded=mul(mag1,r1)^mul(mag2,r2)
            round_negative=mul(rounded,~high)
            rows.append([rounded]+[round_negative]*(n-1))
        return self.sum_rows(rows,n)


def build(carry_save=False):
    table=affine.table();sizes=[int(c.max()).bit_length() for c in table.T]
    n=27+sum(sizes);L=(CarrySaveLogic if carry_save else Logic)(n);inputs=[Signal(frozenset([i])) for i in range(n)]
    x=inputs[:22];flags=inputs[22:27];regs=[];at=27
    for width in sizes:regs.append(inputs[at:at+width]);at+=width
    z=list(x[:14]);z[-1]=~z[-1]
    slope=L.horner(L.resize(regs[1],18),regs[2],z,18-8)
    residual=L.horner(list(regs[0]),slope,z,27-18)
    y=L.resize(residual,27)
    const=[ONE if (11<<23)>>i&1 else ZERO for i in range(27)]
    y[:]=L.add(y,const)
    y[25:]=L.add(y[25:],[flags[0],ZERO])
    q=L.resize(x[14:21],9);y[18:]=L.add(y[18:],[~v for v in q],ONE)
    y=L.resize(y,37)
    for j in range(4):
        shifted=[ZERO]*(1<<j)+y
        y=[v^L.net.mul(flags[j+1],v^s) for v,s in zip(y,shifted)]
    L.net.outputs=y
    return L.net,table,sizes


def evaluate(net,inputs):
    vals=list(inputs)
    def signal(s):
        out=np.full(inputs.shape[1],bool(s.const),dtype=bool)
        for i in s.terms:out^=vals[i]
        return out
    for a,b in net.nodes:vals.append(signal(a)&signal(b))
    out=np.zeros(inputs.shape[1],dtype=np.int64)
    for i,s in enumerate(net.outputs):out|=signal(s).astype(np.int64)<<i
    return out


def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--carry-save',action='store_true');args=p.parse_args()
    net,table,sizes=build(args.carry_save);rng=np.random.default_rng(82417)
    xs=np.unique(np.r_[rng.integers(1,3133453,32768),np.arange(1,256),2**np.arange(22),3133452])
    k=21-np.floor(np.log2(xs)).astype(np.int64);normalized=xs<<k
    address=((normalized>>13)&255)+((k&1)<<8);coefs=table[address]
    bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    inputs=np.concatenate([bits(normalized,22),bits(k,5),*[bits(coefs[:,j],s) for j,s in enumerate(sizes)]])
    got=evaluate(net,inputs);old=poly.table_for;raw=np.load('candidates/margin_curve8_t5.9e-08.npy');poly.table_for=lambda *args:raw
    try:want=poly.model(xs,7,(27,18,8),rounded=True,lookup_frac=8,booth=True)
    finally:poly.table_for=old
    assert np.array_equal(got,want),[(int(xs[i]),int(got[i]),int(want[i])) for i in np.where(got!=want)[0][:8]]
    parents=[sorted({j for s in pair for j in s.terms if j>=net.n}) for pair in net.nodes]
    depth=[];ancestors=[]
    for pp in parents:
        depth.append(1+max([depth[p-net.n] for p in pp],default=0));aa=set(pp)
        for p in pp:aa.update(ancestors[p-net.n])
        ancestors.append(aa)
    report=dict(status='MODEL_PASS',scope='combinational arithmetic equals selected integer phase model; no quantum circuit or accepted score',carry_save=args.carry_save,tested_x=len(xs),input_bits=net.n,coefficient_sizes=sizes,AND_nodes=len(net.nodes),maximum_AND_depth=max(depth),maximum_ancestor_AND_nodes=max(map(len,ancestors)),output_signal_sizes=[len(s.terms) for s in net.outputs])
    tag='boolean_polynomial_phase'+('_carry_save' if args.carry_save else '')
    Path(f'results/{tag}_model.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
