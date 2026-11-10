"""Exact identical-payload subtree probe; no new oracle family is emitted.

The selected emitter is read-only round-4 Codex emit_qrom, reused by qfit.
Its prefix/branch counter is checked against all eleven accepted KMB counts.
Full payload and all-pairs XOR uniqueness give permutation-independent
certificates; nearest-neighbour/2-opt orders are secondary Hamming controls.
"""
from pathlib import Path
from collections import Counter
import hashlib,json,csv,time
import numpy as np
import kickmix as km
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

HERE=Path(__file__).resolve().parent;SOURCE=HERE/'candidates/banked';start=time.time()
z=dict(np.load(SOURCE/'submission.npz'));delta=batch_deltas(z['angle_words'],7)
A=len(z['angle_words']);N=76;B=A-N
assert A==795 and B==719 and delta.shape==(11,795,7)
assert np.all((delta>>7&1)==(delta>>6&1))
packed=np.stack([np.sum((d.astype(np.uint64)&127)<<np.arange(0,49,7,dtype=np.uint64),axis=1,dtype=np.uint64) for d in delta],axis=1)

def qrom_count(table,lowbits=8):
    """Exact Toffoli counter for selected emit_qrom on 9 address bits.

    Prefix has one bit: its compute/uncompute use only Clifford operations.
    Every nonconstant branching node contributes one Toffoli, including the
    two-leaf controlled-difference optimization. Invalid-only halves skip.
    """
    total=0;collapsed=0
    for offset in [0,256]:
        if offset>=len(table):continue
        def node(lo,hi):
            nonlocal collapsed
            if lo+offset>=len(table):return 0
            end=min(hi+offset,len(table));part=table[lo+offset:end]
            if hi-lo==1 or all(v==part[0] for v in part):
                if len(part)>1:collapsed+=1
                return 0
            mid=(lo+hi)//2
            if mid+offset>=len(table):return node(lo,mid)
            return 1+node(lo,mid)+node(mid,hi)
        total+=node(0,256)
    return total,collapsed

def aligned_duplicates(vals):
    """Nontrivial equal-payload aligned blocks, considering valid rows only."""
    out={}
    for size in [2,4,8,16,32,64,128,256,512,1024]:
        count=0
        for lo in range(0,len(vals),size):
            block=vals[lo:min(lo+size,len(vals))]
            if len(block)>1 and np.all(block==block[0]):count+=1
        out[str(size)]=count
    return out

certificates=[]
ix,jx=np.triu_indices(A,1)
for k in range(11):
    v=packed[:,k];pair_xors=v[ix]^v[jx];uniq,c=np.unique(pair_xors,return_counts=True)
    fixed=v[:N:2]^v[1:N:2]
    certificates.append(dict(k=k,valid_payloads=A,unique_payloads=int(len(np.unique(v))),
        distinct_input_pairs=len(pair_xors),unique_pair_xors=len(uniq),pair_xor_duplicate_excess=int((c-1).sum()),
        minimum_pair_hamming=int(np.bitwise_count(pair_xors).min()),
        fixed_onebody_pair_active_bits=int(np.bitwise_or.reduce(fixed)).bit_count(),
        current_aligned_equal_blocks=aligned_duplicates(v)))
assert all(c['unique_payloads']==A and c['pair_xor_duplicate_excess']==0 for c in certificates)

# These row orders are common across all eleven oracles. The first 76
# one-body rows stay fixed; every dictionary row travels with its U and w.
permutations={'current':np.arange(B,dtype=np.int64)}
for k in range(11):permutations[f'lex_oracle_{k:02}']=np.argsort(packed[N:,k],kind='stable')
permutations['lex_all']=np.lexsort(tuple(packed[N:,k] for k in reversed(range(11))))
v=packed[N:];distance=np.zeros((B,B),dtype=np.uint16)
for k in range(11):distance+=np.bitwise_count(v[:,None,k]^v[None,:,k]).astype(np.uint16)
anchor=np.bitwise_count(v^packed[N-1]).sum(1).astype(np.uint16)
def path_cost(order):return int(anchor[order[0]])+int(distance[order[:-1],order[1:]].sum())
rng=np.random.default_rng(20261024);best=None
for seed in [int(anchor.argmin()),*rng.choice(B,7,replace=False).tolist()]:
    path=[seed];used=np.zeros(B,dtype=bool);used[seed]=True
    for _ in range(B-1):
        ds=distance[path[-1]].astype(np.int32);ds[used]=100000
        j=int(ds.argmin());path.append(j);used[j]=True
    path=np.array(path,dtype=np.int64)
    for _ in range(30000):
        i,j=sorted(rng.integers(0,B,size=2).tolist())
        if i==j:continue
        before=int(anchor[path[i]]) if i==0 else int(distance[path[i-1],path[i]])
        after=int(anchor[path[j]]) if i==0 else int(distance[path[i-1],path[j]])
        if j+1<B:
            before+=int(distance[path[j],path[j+1]]);after+=int(distance[path[i],path[j+1]])
        if after<before:path[i:j+1]=path[i:j+1][::-1]
    if best is None or path_cost(path)<path_cost(best):best=path.copy()
permutations['hamming_nn_2opt']=best
np.savez(HERE/'parameters/subtree_permutations.npz',**permutations)

records=[];per_oracle=[]
for name,p in permutations.items():
    assert np.array_equal(np.sort(p),np.arange(B))
    inverse=np.argsort(p)
    assert np.array_equal(z['U'][:,p][:,inverse],z['U'])
    assert np.array_equal(z['W'][:,p][:,inverse],z['W'])
    assert np.array_equal(z['angle_words'][N:][p][inverse],z['angle_words'][N:])
    order=np.r_[np.arange(N),N+p];vv=packed[order]
    saving=[]
    for k in range(11):
        rows=list(map(int,vv[:,k]));rows.append(rows[-1])
        first=rows[::2];diff=[a^b for a,b in zip(first,rows[1::2])]
        # Bijective packing preserves equality; compression of active bit
        # coordinates is also injective, so no compressed table is needed.
        combined=[a|(d<<49) for a,d in zip(first,diff)]
        lookup,collapsed=qrom_count(combined)
        active=int(np.bitwise_or.reduce(np.array(diff,dtype=np.uint64))).bit_count()
        anf=37;projected=lookup+active+anf
        current=km.Circuit.from_file(str(SOURCE/f'angle_k_{k}.kmb'),format='kmb').operation_statistics().num_toffoli_like
        assert current==482 and lookup==396 and active==49 and collapsed==0
        saving.append(int(current-projected))
        per_oracle.append(dict(permutation=name,k=k,current_forward_T=int(current),projected_forward_T=projected,
            lookup_T=lookup,low_bit_selection_T=active,ANF_erasure_T=anf,collapsed_lookup_nodes=collapsed,
            forward_saving=int(current-projected),forward_reverse_saving=int(2*(current-projected))))
    records.append(dict(name=name,common_dictionary_permutation=True,onebody_rows_fixed=True,
        hamming_path_cost=path_cost(p),forward_family_saving=sum(saving),chemistry_T_saving=2*sum(saving)))

threshold=.03*30504
report=dict(status='PROBE_COMPLETE_NO_BUILD',source=str(SOURCE),source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCE.iterdir() if p.is_file()},
    certificates=certificates,permutations=records,per_oracle=per_oracle,current_angle_T=10604,current_chemistry_T=30504,
    optimized_projected_total_saving=0,optimized_projected_saving_fraction=0.,build_threshold_T=threshold,
    minimum_integer_build_threshold_T=int(np.ceil(threshold)),pruned_family_built=False,
    proof='Every full delta payload is unique; moreover all 315615 distinct-input pair XORs are unique in every oracle. Row permutations cannot create equal valid-payload subtrees or equal pair-difference leaves. Fixed one-body pairs force all 49 difference-bank bits active. The selected emitter already skips invalid-only ranges and collapses constant nodes.',
    scope='Exact equality of whole XOR payloads under address permutations; not a bound on other bitwise, gauge, linear or oracle architectures.',seconds=time.time()-start)
(HERE/'results/subtree_probe.json').write_text(json.dumps(report,indent=2)+'\n')
with (HERE/'results/subtree_probe.csv').open('w') as f:
    out=csv.DictWriter(f,fieldnames=list(per_oracle[0]));out.writeheader();out.writerows(per_oracle)
print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256','certificates','per_oracle']},indent=2))
