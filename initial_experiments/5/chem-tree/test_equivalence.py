"""Binary-partition equivalence and mutation controls on native words."""
from core import *
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

def main():
    z,u,w,inc,_=load(INC);rng=np.random.default_rng(90213)
    x=rng.normal(size=(43,N));x/=np.linalg.norm(x,axis=-1,keepdims=True)
    families0=families(u,w,inc);groups={};records=[]
    for name,tr in families0.items():
        words,v=quantize(x,tr);theta,s=encode(x,tr)
        assert np.max(abs(decode(theta,tr)-s[:,None]*x))<1e-13
        assert np.max(abs(v-reconstruct_u_vectors(words/128,tr)))<1e-13
        sig=binary_signature(tr)
        if sig in groups:
            prev=groups[sig];err=float(np.max(abs(abs((prev*v).sum(-1))-1)));assert err<1e-13
        else:groups[sig]=v
        records.append(dict(family=name,**summary(tr)))
    mutations=[]
    for mode in ['orbital','siblings','reparent']:
        tr=inc.copy();seen=set();equiv=0
        for k in range(100):
            new=mutate(tr,rng,mode);validate(new);words,v=quantize(x,new)
            assert np.max(abs(v-reconstruct_u_vectors(words/128,new)))<1e-13
            equiv+=binary_signature(new)==binary_signature(tr);tr=new;seen.add(binary_signature(tr))
        mutations.append(dict(mode=mode,steps=100,distinct_binary_classes=len(seen),equivalent_steps=equiv))
    report(HERE/'results/binary_equivalence_tests.json',dict(families=records,
        distinct_graphs=len({fingerprint(t) for t in families0.values()}),distinct_binary_classes=len(groups),mutations=mutations))
if __name__=='__main__':main()
