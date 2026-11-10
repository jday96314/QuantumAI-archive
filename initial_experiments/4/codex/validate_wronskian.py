"""Check the proposed template identity on the non-streaming q-form BM model.

This does not validate the borrowed-bank layout or its stopped branches.
"""
import argparse,json,random,time
from probe_wronskian import p,mul,model,evaluate
parser=argparse.ArgumentParser();parser.add_argument('--shots',type=int,default=1000);a=parser.parse_args()
rng=random.Random(481123);counts={};failures=[];began=time.time()
for shot in range(a.shots):
    weight=shot%36;positions=rng.sample(range(4095),weight);values=[rng.randrange(1,4096) for _ in positions]
    S=[]
    for k in range(1,71):
        value=0
        for pos,v in zip(positions,values):value^=mul(v,p.EXP[k*pos%4095])
        S.append(value)
    sig,bank,d=model(S);deff=d if d else -1;rotation=deff%37;bank=bank[rotation:]+bank[:rotation]
    H=[0]*36
    for i,x in enumerate(sig):
        for j,y in enumerate(bank):
            if (i+j)&1:H[(i+j-1)//2]^=mul(x,y)
    H=[p.EXP[p.LOG[x]*2048%4095] if x else 0 for x in H]
    good=bank[0]==0 and not any(H[weight+1:])
    for pos,value in zip(positions,values):
        x=p.EXP[-pos%4095];h=evaluate(H,x)
        hn=mul(h,p.EXP[weight*pos%4095])
        good &= evaluate(sig,x)==0 and mul(value,mul(hn,hn))==1
    counts[weight]=counts.get(weight,0)+1
    if not good:failures.append(dict(shot=shot,weight=weight,d=d))
result=dict(status='PASS' if not failures else 'FAIL',scope='non-streaming q-form BM classical model',shots=a.shots,weight_counts=counts,failures=failures,seconds=time.time()-began)
print(json.dumps(result,indent=2))
