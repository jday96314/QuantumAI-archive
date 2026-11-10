"""Scan the chemistry cost model (T*Q, lambda excluded) over factorization shapes."""
import math, itertools
from dfthcsos.verify import heuristic_costs
N=76
def angles(R,B,brot,J,temp=8):
    rows=N+R*B; sel=math.ceil(math.log2(rows))
    if sel+J*brot+temp>128: return None
    K=math.ceil(75/J)
    return dict(T=2*K*(rows-1), persistent=J*brot, temp=temp, q=sel+J*brot+temp)
def cost(R,B,C,brot,J):
    a=angles(R,B,brot,J)
    if a is None: return None
    try: h=heuristic_costs(N,R,B,C,brot)
    except Exception as e: return None
    T=a['T']+h['toffoli']; Q=h['persistent_qubits']+a['persistent']+max(h['temporary_qubits'],a['temp'])
    return T,Q,h
if __name__=='__main__':
    base=cost(1,719,350,8,13); print('base',base[0],base[1],base[0]*base[1])
    out=[]
    for R in [1,2,3,4,5,6,8]:
        for tot in [600,720,800,900,1000]:
            B=tot//R
            for Ctot in [300,350]:
                C=Ctot//R
                for brot in [7,8,9,10,12]:
                    for J in range(4,17):
                        r=cost(R,B,C,brot,J)
                        if r: out.append((r[0]*r[1],R,B,C,brot,J,r[0],r[1]))
    out.sort()
    seen=set()
    for o in out:
        key=o[1:5]
        if key in seen: continue
        seen.add(key)
        if len(seen)>60: break
        print(o)
