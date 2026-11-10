"""Necessary matrix-rank bound for every native chemistry factor gauge.

For symmetric rank<=B matrices X, (Tr X)^2 <= B*||X||_F^2.
Summing over factor columns and allowing the full physical Frobenius error
gives a dictionary-independent necessary bound on B. This does not bound
circuit costs or rule out the competition target.
"""
import json
from pathlib import Path
import numpy as np


def main():
 sp=np.load('../../3/parameters/reference_spectrum.npz')
 eig,vec=sp['eigenvalues'],sp['eigenvectors']
 d=(sp['ii']==sp['jj']).astype(float);n=int(d.sum());D=len(d);epsilon=.02
 trace=float(eig.sum());overlap=float(np.sum(eig*(d@vec)**2))
 lower=(overlap-epsilon*n)/(trace+epsilon*np.sqrt(D))
 # Independently check both contractions against an explicitly formed G.
 G=(vec*eig)@vec.T
 assert abs(float(np.trace(G))-trace)<1e-9
 assert abs(float(d@G@d)-overlap)<1e-9
 report=dict(status='PASS',scope='necessary maximum physical matrix rank B for all factor column gauges, dictionaries, R and C; no circuit-cost or overall-score bound',epsilon=epsilon,n=n,packed_dimension=D,trace_G=trace,identity_overlap=overlap,zero_error_rank_bound=overlap/trace,error_aware_rank_bound=float(lower),necessary_integer_B=int(np.ceil(lower)),derivation='sum (Tr X_c)^2 <= B sum ||X_c||_F^2; numerator changes by at most epsilon*n, denominator by epsilon*sqrt(D)')
 Path('results/physical_factor_rank_bound.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report))


if __name__=='__main__':main()
