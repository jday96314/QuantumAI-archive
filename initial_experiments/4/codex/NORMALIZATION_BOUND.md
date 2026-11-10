# A global lower bound for the native normalization factor

The verifier fixes S and beta1 for every submission. Define

    c = E_upper - ecore + eta*beta1 = -559.1324731485038.

For each factor column write t=sum_b w_b and l=sum_b abs(w_b), so l>=abs(t).
Let h=h1_prime and u=2*Lambda-Egap. Directly from the native definitions,

    u = -c + 2 sum_j max(eigenvalue_j(h),0)
        + sum_columns [ (t²+l²)/2 + W*t + l*abs(W+t) ].

The bracket is convex in W and minimized at W=-t. Its minimum is
(l²-t²)/2>=0. Consequently u>=-c for every factorization, regardless of R,B,C
or the angle dictionary. The strict native gap gate gives

    lambda_eff² = Egap*u > 3.914898 * 559.1324731485038,
    lambda_eff > 46.78617959252637.

The selected normalization48.59417497035283 is within3.87% of this bound.
This bound limits the normalization factor. Gate and qubit costs still depend
on the factorization and circuit architecture and remain optimization targets.
`results/normalization_global_bound.json` records the constants derived from
the actual native fixed inputs.

`upper_gap_offsets.py` independently minimizes u for fixed selected factors.
The convex probe reaches u598.79455, but Egap rises to roughly169 and the
normalization product worsens. A601-point path from the original offsets to
and beyond the best upper-gap offsets has its best valid product at the original
point. See `results/upper_gap_offsets_probe.json`. This numerical probe concerns
those fixed factors; the analytic inequality above applies to all native factors.

There is also a necessary bound on the matrix rank of physical factor
columns, independent of their gauge or dictionary. Each column represents a
symmetric76-by-76 matrix X of rank at most B, hence
(Tr X)²<=B*||X||F². Let G be the reference covariance in the orthonormal packed
pair basis, d the packed identity and D=2926. Any fit within epsilon=.02 obeys

    B >= (d^T G d - epsilon*76)/(Tr G + epsilon*sqrt(D))
      = 9.56859572434093.

Thus every valid native factorization needs B>=10. This rules out the proposed
separate B2/B4/B8 factors even after jointly optimizing their column gauge.
It does not rule out other factorizations or bound the overall competition
score. `physical_factor_rank_bound.py` independently checks both reference
contractions; `results/physical_factor_rank_bound.json` stores the evidence.
