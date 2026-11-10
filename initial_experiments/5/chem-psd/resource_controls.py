"""Compile and simulate constructive ANGLE bounds for the main row counts."""
from core import *
from signed_oracle import build,bound
from oracle_controls import check

def main():
    rng=np.random.default_rng(607103);records=[]
    for A in [796,972,778,844]:
        bits,J=8,7;ob=bound(A,bits,J);delta=rng.integers(0,128,(A,J));delta|=((delta>>6)&1)<<7
        c=build(delta,bits,ob['lowbits']);T=int(c.operation_statistics().num_toffoli_like)
        assert T==ob['T_per_oracle'] and c.num_qubits<=ob['Q_oracle']
        records.append(dict(A=A,bound=ob,simulations=[check(c,delta,bits,seed,2048) for seed in [1723+A,95171+A]]))
    report(HERE/'results/resource_controls.json',dict(scope='Synthetic XOR/phase/scratch controls, not chemistry submissions',controls=records))
if __name__=='__main__':main()
