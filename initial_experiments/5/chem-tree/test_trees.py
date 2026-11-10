"""Independent decoder, encoder, host validation and gauge controls."""
from core import *
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.verify import evaluate_words
from scipy.linalg import expm

def main():
    rng=np.random.default_rng(2374);z,u,w,inc,words=load(INC)
    ts=families(u,w,inc);records=[]
    physical=rng.normal(size=(31,N));physical/=np.linalg.norm(physical,axis=1,keepdims=True)
    physical=np.r_[physical,np.eye(N)]
    for name,tree in ts.items():
        theta,sign=encode(physical,tree);back=decode(theta,tree)
        rt=float(np.max(np.abs(back-sign[:,None]*physical)))
        wordsq,q=quantize(physical,tree)
        host=reconstruct_u_vectors(wordsq/128,tree)
        hostdiff=float(np.max(np.abs(host-q)))
        # Independently multiply full elementary rotation matrices.
        mat=np.eye(N)
        for j,(p,c) in enumerate(tree):
            a=wordsq[3,j]*np.pi/128;R=np.eye(N)
            R[p,p]=R[c,c]=np.cos(a);R[p,c]=-np.sin(a);R[c,p]=np.sin(a)
            mat=R@mat
        dense=float(np.max(np.abs(mat[:,0]-q[3])))
        # Positive versus signed radius gauges are the SAME projector lattice.
        wp,qp=quantize(physical[:31],tree,signed=False)
        dots=np.abs((qp*q[:31]).sum(1));gauge=float(np.max(np.abs(dots-1)))
        canon=signed_words(wordsq,tree)
        qc=decode(canon*np.pi/128,tree)
        canonerr=float(np.max(np.abs(np.abs(np.sum(qc*q,axis=1))-1)))
        assert rt<3e-14 and hostdiff<3e-14 and dense<3e-14 and gauge<3e-14 and canonerr<3e-14
        assert np.all((canon>>7&1)==(canon>>6&1))
        # Invoke unchanged physical evaluator to exercise its tree validation.
        allwords=np.zeros((N+719,75),dtype=np.int64)
        evaluate_words(native(),allwords,tree,8)
        records.append(dict(name=name,**summary(tree),roundtrip_max=rt,host_decoder_max=hostdiff,
                            matrix_decoder_max=dense,gauge_projector_error=gauge,canonical_error=canonerr))
    # Accepted exact word calibration.
    raw=np.load(INC);q=decode(raw['angle_words']*np.pi/128,inc)
    recanon=signed_words(raw['angle_words'],inc)
    calibration=dict(native_metrics=evaluate_words(native(),raw['angle_words'],inc,8),
                     pair_error=physical_error(u,w),word_changes=int(np.sum(recanon!=raw['angle_words'])),
                     U_projector_max=float(np.max(np.abs(np.abs((q[N:].reshape(u.shape)*u).sum(-1))-1))))
    assert calibration['word_changes']==0
    assert abs(calibration['pair_error']-calibration['native_metrics']['error_2'])<1e-12
    report(HERE/'results/encoding_tests.json',dict(trees=records,calibration=calibration))

if __name__=='__main__':main()
