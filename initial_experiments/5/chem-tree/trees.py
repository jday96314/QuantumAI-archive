"""Legal shared ordered-tree representations; physical labels never change.

Host conventions: from_host/verifiers/dfthcsos/scoreinfra/dfthc_math.py
(Apache-2.0). All encoders and generators here are experiment-local.
"""
import hashlib
import heapq
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

N=76
def validate(tree,n=N):
    tree=np.asarray(tree)
    assert tree.shape==(n-1,2) and np.issubdtype(tree.dtype,np.integer)
    assert np.all((tree>=0)&(tree<n))
    seen={0}
    for p,c in tree:
        assert int(p) in seen and int(c) not in seen
        seen.add(int(c))
    assert len(seen)==n
    return tree

def decode(theta,tree):
    """Independent forward decoder: theta is radians, arbitrary batch shape."""
    tree=validate(tree,theta.shape[-1]+1)
    v=np.zeros((*theta.shape[:-1],len(tree)+1));v[...,0]=1
    for j,(p,c) in enumerate(tree):
        a=v[...,p].copy();b=v[...,c].copy()
        co=np.cos(theta[...,j]);si=np.sin(theta[...,j])
        v[...,p]=co*a-si*b;v[...,c]=si*a+co*b
    return v

def encode(u,tree,signed=True):
    """Exact continuous inverse, returning angles and global root sign.

    Both positive-radius and signed-radius gauges are legal. Branches differ by
    pi and induce only exact downstream sign flips on an even uniform grid.
    """
    v=np.array(u,dtype=float,copy=True);v/=np.linalg.norm(v,axis=-1,keepdims=True)
    theta=np.zeros((*v.shape[:-1],len(tree)))
    for j in reversed(range(len(tree))):
        p,c=tree[j];a=v[...,p].copy();b=v[...,c].copy()
        sign=np.where(a<0,-1.,1.) if signed else np.ones(a.shape)
        theta[...,j]=np.arctan2(b*sign,a*sign)
        v[...,p]=sign*np.hypot(a,b);v[...,c]=0
    return theta,np.where(v[...,0]<0,-1.,1.)

def words_from_angles(theta,bits=8):
    return np.rint(theta*(1<<(bits-1))/np.pi).astype(np.int64)%(1<<bits)

def quantize(u,tree,bits=8,signed=True):
    words=integer_encode(u,tree,bits,signed)
    return words,decode(words*np.pi/(1<<(bits-1)),tree)

def integer_encode(u,tree,bits=8,signed=True):
    # At a rounded +pi/2 boundary, canonicalize to -pi/2 AND propagate
    # the changed radius sign toward the root. Editing words afterward
    # without that propagation changes the physical projector.
    v=np.array(u,dtype=float,copy=True);words=np.zeros((*v.shape[:-1],len(tree)),dtype=np.int64)
    scale=1<<(bits-1)
    for j in reversed(range(len(tree))):
        p,c=tree[j];a=v[...,p].copy();b=v[...,c].copy()
        sign=np.where(a<0,-1.,1.) if signed else np.ones(a.shape)
        theta=np.arctan2(b*sign,a*sign);word=np.rint(theta*scale/np.pi).astype(np.int64)
        if signed:
            flip=word>=scale//2;word=np.where(flip,word-scale,word);sign=np.where(flip,-sign,sign)
            flip=word<-(scale//2);word=np.where(flip,word+scale,word);sign=np.where(flip,-sign,sign)
        words[...,j]=word%(2*scale);v[...,p]=sign*np.hypot(a,b);v[...,c]=0
    return words

def signed_words(words,tree,bits=8):
    """Canonical signed extension: rounded boundary +64 maps to -64."""
    return integer_encode(decode(words*np.pi/(1<<(bits-1)),tree),tree,bits,True)

def fingerprint(tree):return hashlib.sha256(np.asarray(tree,dtype=np.int64).tobytes()).hexdigest()

def binary_signature(tree):
    """Invariant of the induced binary amplitude-splitting tree.

    Representative relabels within a binary split use theta -> pi/2-theta
    and sign gauges. Those maps preserve an even uniform angle lattice.
    Distinct graph fingerprints can therefore describe equivalent lattices.
    """
    groups=[1<<i for i in range(len(tree)+1)];clades=[]
    for p,c in tree[::-1]:
        groups[p]|=groups[c];clades.append(groups[p]);groups[c]=0
    return hashlib.sha256(','.join(map(str,sorted(clades))).encode()).hexdigest()

def nested_tree(node):
    edges=[]
    def leaves(x):return [x] if isinstance(x,(int,np.integer)) else leaves(x[0])+leaves(x[1])
    def emit(x,rep):
        if isinstance(x,(int,np.integer)):return
        a,b=x
        if rep not in leaves(a):a,b=b,a
        new=leaves(b)[0];edges.append((rep,new));emit(a,rep);emit(b,new)
    emit(node,0)
    return validate(np.asarray(edges,dtype=np.int64))

def huffman(power):
    heap=[(float(p),i,i) for i,p in enumerate(power)];heapq.heapify(heap);uid=len(power)
    while len(heap)>1:
        pa,_,a=heapq.heappop(heap);pb,_,b=heapq.heappop(heap)
        heapq.heappush(heap,(pa+pb,uid,(a,b)));uid+=1
    return nested_tree(heap[0][2])

def binary_from_tree(tree):
    nodes=list(range(len(tree)+1))
    for p,c in tree[::-1]:nodes[p]=(nodes[p],nodes[c]);nodes[c]=None
    return nodes[0]

def nni_neighbors(tree,max_support=8):
    """Nearest-neighbor binary split rotations, affecting only local groups.

    A,(B,C) -> (A,B),C or (A,C),B. Unlike arbitrary orbital swaps, only
    two continuous angles per row change (up to exact split gauge maps).
    """
    root=binary_from_tree(tree);out={}
    def size(x):return 1 if isinstance(x,(int,np.integer)) else size(x[0])+size(x[1])
    def variants(x):
        if isinstance(x,(int,np.integer)):return []
        a,b=x;rows=[]
        if size(x)<=max_support:
            if not isinstance(a,(int,np.integer)):
                c,d=a;rows.extend([(c,(d,b)),(d,(c,b))])
            if not isinstance(b,(int,np.integer)):
                c,d=b;rows.extend([((a,c),d),((a,d),c)])
        rows.extend([(v,b) for v in variants(a)]);rows.extend([(a,v) for v in variants(b)])
        return rows
    old=binary_signature(tree)
    for node in variants(root):
        tr=nested_tree(node);sig=binary_signature(tr)
        if sig!=old:out[sig]=tr
    return list(out.values())

def summary(tree):
    depth={0:0};degree=np.zeros(len(tree)+1,dtype=int)
    for p,c in tree:depth[int(c)]=depth[int(p)]+1;degree[p]+=1
    return dict(sha256=fingerprint(tree),binary_signature=binary_signature(tree),max_depth=max(depth.values()),
                average_depth=float(np.mean(list(depth.values()))),max_children=int(max(degree)))

def ordered(order,kind='balanced'):
    order=list(map(int,order));assert order[0]==0 and set(order)==set(range(len(order)))
    if kind=='star':edges=[(0,c) for c in order[1:]]
    elif kind=='chain':edges=list(zip(order[:-1],order[1:]))
    elif kind=='heap':edges=[(order[(i-1)//2],order[i]) for i in range(1,len(order))]
    else:
        edges=[]
        def rec(group):
            if len(group)<2:return
            m=(len(group)+1)//2;left,right=group[:m],group[m:]
            edges.append((left[0],right[0]));rec(left);rec(right)
        rec(order)
    return validate(np.asarray(edges,dtype=np.int64),len(order))

def clustered(distance,method='average',reverse=False):
    """Binary orbital clusters become legal rotations between representatives."""
    root=to_tree(linkage(squareform(distance,checks=False),method=method));edges=[]
    def leaves(node):return [node.id] if node.is_leaf() else leaves(node.left)+leaves(node.right)
    def emit(node,rep):
        if node.is_leaf():return
        a,b=node.left,node.right
        if rep not in leaves(a):a,b=b,a
        choices=leaves(b);new=choices[-1] if reverse else choices[0]
        edges.append((rep,new));emit(a,rep);emit(b,new)
    emit(root,0)
    return validate(np.asarray(edges,dtype=np.int64))

def families(u,w,incumbent,seed=428):
    """All families share one tree over all vectors/blocks, no row-specific tree."""
    rng=np.random.default_rng(seed);u=u.reshape(-1,N)
    leverage=np.sum(w*w,axis=-1).reshape(-1);leverage/=leverage.mean()
    power=np.sum(leverage[:,None]*u*u,axis=0)
    order=np.r_[0,1+np.argsort(-power[1:])]
    out={'incumbent':incumbent,'balanced_identity':ordered(range(N)),
         'balanced_power':ordered(order),'heap_power':ordered(order,'heap'),
         'chain_power':ordered(order,'chain'),'star_power':ordered(order,'star')}
    out['huffman_amplitude']=huffman(power)
    P=u[:,np.tril_indices(N)[0]]*u[:,np.tril_indices(N)[1]]*np.where(np.tril_indices(N)[0]==np.tril_indices(N)[1],1,np.sqrt(2))
    # Z_i=sum_c w_ic F_c, with block membership retained. First isotropic
    # tensor-Jacobian term is 4||Z_i||^2 * subtree probability mass.
    uu=u.reshape(w.shape[0],w.shape[1],N);PP=P.reshape(w.shape[0],w.shape[1],-1)
    Z=w@(w.transpose(0,2,1)@PP);sensitivity=np.sum(Z*Z,axis=-1).reshape(-1)
    out['huffman_tensor']=huffman(np.sum(sensitivity[:,None]*u*u,axis=0))
    # Absolute amplitudes capture support; signed orbital correlations capture
    # common relative sign. Distances are physical-coordinate statistics only.
    for name,V in [('amplitude',np.abs(u)*np.sqrt(leverage[:,None])),
                   ('correlation',u*np.sqrt(leverage[:,None])),
                   ('squared',u*u*np.sqrt(leverage[:,None]))]:
        norm=np.linalg.norm(V,axis=0);corr=V.T@V/(norm[:,None]*norm[None,:]+1e-30)
        distance=np.sqrt(np.maximum(0,2-2*np.abs(corr)));np.fill_diagonal(distance,0)
        out['cluster_'+name]=clustered(distance)
        out['cluster_'+name+'_reverse']=clustered(distance,reverse=True)
    for k in range(4):
        perm=np.r_[0,rng.permutation(np.arange(1,N))]
        out[f'random_balanced_{k}']=ordered(perm)
    return out

def physical_families(H):
    """Reference-Hamiltonian coordinate statistics, never a tensor permutation."""
    ii,jj=np.tril_indices(N);diag=np.flatnonzero(ii==jj)
    V=H[:,diag];norm=np.linalg.norm(V,axis=0)
    corr=V.T@V/(norm[:,None]*norm[None,:]+1e-30)
    distance=np.sqrt(np.maximum(0,2-2*np.abs(corr)));np.fill_diagonal(distance,0)
    out={'host_density_clusters':clustered(distance)}
    exchange=np.zeros((N,N));exchange[ii,jj]=np.abs(np.diag(H))/np.where(ii==jj,1,2)
    exchange+=exchange.T;exchange[np.diag_indices(N)]*=.5
    scale=np.sqrt(np.maximum(1e-30,np.diag(exchange)))
    corr=np.minimum(1,exchange/(scale[:,None]*scale[None,:]))
    distance=np.sqrt(np.maximum(0,2-2*corr));np.fill_diagonal(distance,0)
    out['host_exchange_clusters']=clustered(distance)
    return out

def mutate(tree,rng,mode):
    """Coordinate permutation keeps root fixed; reparenting is re-topologized."""
    tree=tree.copy()
    if mode=='orbital':
        a,b=rng.choice(np.arange(1,N),2,replace=False)
        tree=np.where(tree==a,-1,tree);tree=np.where(tree==b,a,tree);tree[tree==-1]=b
    elif mode=='siblings':
        # Move any outgoing rotation relative to another outgoing rotation.
        parents,counts=np.unique(tree[:,0],return_counts=True);parents=parents[counts>1]
        if not len(parents):return mutate(tree,rng,'orbital')
        p=rng.choice(parents);a,b=rng.choice(np.where(tree[:,0]==p)[0],2,replace=False)
        # Reorder sibling subtrees by rebuilding a preorder walk with shuffled
        # children; a bare edge swap could introduce a child before its parent.
        children={i:list(tree[tree[:,0]==i,1]) for i in range(N)}
        ca,cb=children[p].index(tree[a,1]),children[p].index(tree[b,1])
        children[p][ca],children[p][cb]=children[p][cb],children[p][ca]
        rows=[]
        def walk(p):
            for c in children[p]:rows.append((p,c));walk(c)
        walk(0);tree=np.asarray(rows,dtype=np.int64)
    else:
        child=int(rng.integers(1,N));desc={child}
        for _ in range(N):
            new=desc|set(map(int,tree[np.isin(tree[:,0],list(desc)),1]))
            if new==desc:break
            desc=new
        parent=int(rng.choice([i for i in range(N) if i not in desc]))
        tree[tree[:,1]==child,0]=parent
        remaining=list(map(tuple,tree));rows=[];seen={0}
        while remaining:
            for j,(p,c) in enumerate(remaining):
                if p in seen:rows.append((p,c));seen.add(c);remaining.pop(j);break
        tree=np.asarray(rows,dtype=np.int64)
    return validate(tree)
