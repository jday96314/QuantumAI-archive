"""Metric and dual-bound controls, including native full-tensor residual."""
from core import *
from solve import solve
from dfthcsos.verify import evaluate_words

def main():
    inp=load_inputs(INC);z=dict(np.load(INC));u=reconstruct_u_vectors(z['angle_words'][76:]/128,z['tree_edges']).reshape(1,719,76)
    P=features(u);X=inp.w[0].T@P[0];approx=X.T@X;pair=float(np.linalg.norm(approx-target()))
    native=evaluate_words(inp,z['angle_words'],z['tree_edges'],8)
    full=unpack(X).reshape(350,76,76)
    h=np.einsum('apq,ast->pqst',full,full,optimize=True)
    ident=np.eye(76);h-=.5*(np.einsum('pq,st->pqst',inp.S,ident)+np.einsum('pq,st->pqst',ident,inp.S))
    direct=float(np.linalg.norm(h-inp.ref['tei']))
    assert abs(pair-direct)<1e-10 and abs(pair-native['error_2'])<1e-7
    gram=float(np.max(abs(P[0]@P[0].T-(u[0]@u[0].T)**2)))
    assert gram<1e-12
    # R1 has a closed-form optimum. The dual certificate must reproduce it.
    m=Covariance(u[:,:24]);S=project(m.A);diag=m.diagnostics(S,True)
    assert diag['primal_dual_gap_squared']<1e-7
    # Check derivative of the contracted objective in feasible coordinates.
    rng=np.random.default_rng(81);V=tensor(rng.normal(size=(1,24,24)));V=(V+V.transpose(1,2))/2
    eps=1e-5;fd=float((m.objective(S+eps*V)-m.objective(S-eps*V))/(2*eps));an=float(2*((m.apply(S)-m.A)*V).sum())
    assert abs(fd-an)<1e-6
    ev=np.linalg.eigvalsh(target())
    report(HERE/'results/calibration.json',dict(incumbent_native_metrics=native,pair_error_2=pair,full_tensor_error_2=direct,
       projector_gram_max_discrepancy=gram,reference_min_eigenvalue=float(ev.min()),negative_spectrum_norm=float(np.linalg.norm(ev[ev<0])),
       single_block_analytic_dual_control=diag,gradient_finite_difference=fd,gradient_analytic=an))
if __name__=='__main__':main()
