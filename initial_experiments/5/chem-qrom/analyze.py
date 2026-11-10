"""Reproducible content statistics, scoped to completed references."""
from common import *
from linear import basis
from erase import coefficients
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def pack(d):return [sum((int(v)&127)<<(7*j) for j,v in enumerate(row)) for row in d]

def distances(words):
    p=pack(words);return dict(mean_adjacent_hamming=float(np.mean([(a^b).bit_count() for a,b in zip(p,p[1:])])),
        mean_paired_hamming=float(np.mean([(a^b).bit_count() for a,b in zip(p[::2],p[1::2])])) )

def run():
    words=np.load(INC/'submission.npz')['angle_words'];rows=[]
    for k,d in enumerate(batch_deltas(words,7)):
        table=pack(d);table.append(table[-1]);diff=[a^b for a,b in zip(table[::2],table[1::2])]
        coeff=coefficients(table[:len(words)],10);hist={}
        for m,c in enumerate(coeff):
            if c:hist[m.bit_count()]=hist.get(m.bit_count(),0)+1
        ent=[]
        for j in range(7):
            _,cnt=np.unique(d[:,j],return_counts=True);p=cnt/cnt.sum();ent.append(float(-(p*np.log2(p)).sum()))
        # Complete multi-output subtrees, compared by their actual payload.
        subtree={}
        for size in [2,4,8,16]:
            nodes=[tuple(table[off:off+size]) for off in range(0,len(table),size)]
            subtree[size]=dict(nodes=len(nodes),unique=len(set(nodes)),constant=sum(len(set(x))==1 for x in nodes))
        rows.append(dict(k=k,entropy=ent,zero_fraction=float((d==0).mean()),difference_rank=len(basis(diff)),
                         active_difference_bits=int(np.bitwise_or.reduce(np.array(diff,dtype=object))).bit_count(),
                         unique_pair_differences=len(set(diff)),nonzero_ANF=len([c for c in coeff if c]),degree_histogram=hist,
                         repeated_subtrees=subtree))
    permutations={'baseline':distances(words)}
    for name in ['greedy7','lex7','parity0','hybrid4_o0']:
        f=HERE/'candidates'/name/'submission.npz'
        if f.exists():permutations[name]=distances(np.load(f)['angle_words'])
    save('content_analysis.json',dict(shape=list(words.shape),sign_extension_exceptions=int(np.count_nonzero((words>>7&1)!=(words>>6&1))),
                                   rows=rows,permutation_distances=permutations))

if __name__=='__main__':run()
