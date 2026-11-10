"""Independent NumPy/host-decoder audit of saved full-PSD dual certificates."""
from core import *
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

def main():
    H=target();rows=[]
    for tag in ['psd_repair_nni_b479','psd_repair_original_b479']:
        z=np.load(HERE/f'candidates/{tag}_dual.npz');d=z['D'];tr=z['tree_edges'];words=z['angle_words']
        u=reconstruct_u_vectors(words.reshape(-1,75)/128,tr).reshape(2,479,N)
        P=features(u);maximum=[];orth=[]
        for p in P:
            q,_=np.linalg.qr(p.T);orth.append(float(np.max(np.abs(q.T@q-np.eye(q.shape[1])))))
            maximum.append(float(np.linalg.eigvalsh(q.T@d@q)[-1]))
        lower2=float(2*np.sum(H*d)-np.sum(d*d));prior=json.loads((HERE/f'results/{tag}_final.json').read_text())
        assert max(maximum)<-9e-12 and lower2>.020**2
        assert abs(np.sqrt(lower2)-prior['dual_lower_bound'])<1e-9
        _,uu,ww,_,_=load(HERE/f'candidates/{tag}.npz');primal=physical_error(uu,ww,H)
        assert abs(primal-prior['full_PSD_error'])<1e-10 and np.sqrt(lower2)<=primal
        rows.append(dict(tag=tag,dual_npz_sha256=sha(HERE/f'candidates/{tag}_dual.npz'),
            block_max_eigenvalues=maximum,orthogonality_errors=orth,primal_error=primal,
            dual_lower_bound=np.sqrt(lower2),fixed_dictionary_excluded=True,status='PASS'))
    report(HERE/'results/independent_dual_audit.json',rows)
if __name__=='__main__':main()
