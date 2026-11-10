"""Route degree-selected template permutations with a small switch bank.

This is a research primitive, not yet integrated into a submitted decoder.
"""
import json
import kickmix as km
from probe_stream_template import normalize_bank
from monomial_unlookup import erase


def route(permutation,wires=None):
    """A recursive two-way rearrangeable network; permutation maps input to output."""
    n=len(permutation);wires=list(range(n)) if wires is None else list(wires)
    if n<=1:return []
    if n==2:return [[(wires[0],wires[1],int(permutation[0]==1))]]
    inverse=[0]*n
    for i,j in enumerate(permutation):inverse[j]=i
    colors=[None]*n
    def fill(start):
        colors[start]=0;stack=[start]
        while stack:
            i=stack.pop();neighbors=[]
            if (i^1)<n:neighbors.append(i^1)
            if (permutation[i]^1)<n:neighbors.append(inverse[permutation[i]^1])
            for j in neighbors:
                if colors[j] is None:colors[j]=1-colors[i];stack.append(j)
                else:assert colors[j]!=colors[i]
    if n%2:
        fill(n-1)
        assert colors[inverse[n-1]]==0
    for i in range(n):
        if colors[i] is None:fill(i)
    first=[(wires[i],wires[i+1],colors[i]) for i in range(0,n-1,2)]
    last=[(wires[j],wires[j+1],colors[inverse[j]]) for j in range(0,n-1,2)]
    sub=[[-1]*((n+1-c)//2) for c in [0,1]]
    for i,j in enumerate(permutation):sub[colors[i]][i//2]=j//2
    children=[route(sub[c],wires[c::2]) for c in [0,1]]
    inner=[]
    for depth in range(max(map(len,children))):inner.append([item for child in children if depth<len(child) for item in child[depth]])
    return [first]+inner+[last]


def plans():
    result=[]
    for degree in range(64):
        arrangement=normalize_bank(list(range(37)),degree) if degree<=35 else list(range(37))
        permutation=[0]*37
        for dest,src in enumerate(arrangement):permutation[src]=dest
        result.append(route(permutation))
    return result


def apply(b,words,selection,inverse=False):
    allplans=plans();layers=len(allplans[0]);maxwidth=max(len(layer) for layer in allplans[0]);controls=b.alloc_qubits(max(8,maxwidth))
    for depth in (reversed(range(layers)) if inverse else range(layers)):
        topology=[(i,j) for i,j,_ in allplans[0][depth]]
        assert all([(i,j) for i,j,_ in plan[depth]]==topology for plan in allplans)
        table=[sum(int(value)<<j for j,(_,_,value) in enumerate(plan[depth])) for plan in allplans]
        b.init_lookup(table=table,address=selection,target=controls)
        for index,(i,j) in enumerate(topology):b.cswap(controls[index],words[i],words[j])
        erase(b,selection,controls,table)
    b.free(controls)


def test():
    import numpy as np
    rng=np.random.default_rng(15732)
    for n in range(2,65):
        for _ in range(10):
            perm=rng.permutation(n).tolist();values=list(range(n))
            for layer in route(perm):
                for i,j,control in layer:
                    if control:values[i],values[j]=values[j],values[i]
            assert all(values[perm[i]]==i for i in range(n)),n
    b=km.CircuitBuilder();selection=b.create_quantum_register(6,'L');words=[b.create_quantum_register(12,f'w{i}') for i in range(37)]
    apply(b,words,selection);c=b.finish_circuit();N=512;degrees=np.arange(N)%64;data=rng.integers(0,4096,(37,N))
    sim=km.Simulator(batch_size=N,seed=5192);sim.use_same_registers_as(c)
    bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    sim.write_across_shots(selection,bits(degrees,6))
    for reg,value in zip(words,data):sim.write_across_shots(reg,bits(value,12))
    sim.do(c);expected=np.stack([normalize_bank(data[:,s].tolist(),int(degrees[s])) if degrees[s]<=35 else data[:,s] for s in range(N)],axis=1)
    for reg,value in zip(words,expected):assert np.array_equal(sim.read_across_shots(reg),bits(value,12))
    assert np.array_equal(sim.read_across_shots(selection),bits(degrees,6))
    public=set(selection)|{q for word in words for q in word};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any()
    assert not sim.read_phase_flipped_across_shots().any()
    print(json.dumps(dict(status='PASS',random_permutations=630,shots=N,layers=len(plans()[0]),Q=c.num_qubits,extra_qubits=c.num_qubits-450,T=c.operation_statistics().num_toffoli_like)))

if __name__=='__main__':test()
