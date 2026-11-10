"""Bind the actual SOS/Lambda offset formulas to native baseline semantics."""
from physical import *

def main():
    inp=load_inputs(INC);z=np.load(INC);u=reconstruct_u_vectors(z['angle_words'][76:]/128,z['tree_edges']).reshape(1,719,76)
    X=reconstruct_x(u,inp.w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).ravel()
    h0=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],np.zeros_like(inp.W),X)
    const=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1;o=inp.W.ravel()
    gap,lam=offset_objective(o,xf,tr,h0,const,inp.w,False);native=evaluate_words(inp,z['angle_words'],z['tree_edges'],8)
    assert abs(gap-native['E_gap'])<1e-8 and abs(lam-native['Lambda'])<1e-8
    rng=np.random.default_rng(7041);_,grad=offset_objective(o,xf,tr,h0,const,inp.w);checks=[]
    for _ in range(3):
        v=rng.normal(size=len(o));v/=np.linalg.norm(v);eps=1e-5
        fd=(offset_objective(o+eps*v,xf,tr,h0,const,inp.w)[0]-offset_objective(o-eps*v,xf,tr,h0,const,inp.w)[0])/(2*eps)
        checks.append(dict(finite_difference=fd,analytic=float(grad@v),absolute_error=abs(fd-grad@v)))
    assert max(c['absolute_error'] for c in checks)<1e-6
    report(HERE/'results/offset_controls.json',dict(gap=gap,Lambda=lam,native_gap=native['E_gap'],native_Lambda=native['Lambda'],gradient_checks=checks))
if __name__=='__main__':main()
