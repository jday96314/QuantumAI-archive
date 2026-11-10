# Phase-aware coordinates and reversible recurrences

All fields below have characteristic two. The independent reference uses
GF(16)/0x13 and host GF(4096)/0x1053. Native arithmetic uses the existing
read-only isomorphism to GF(4096)/0x10eb. Classical masks are runtime inputs;
no construction uses the verifier's error vector except as a test oracle.

## 1. A decoding tangent that is already the evaluator

Write the syndrome series as

\[
S(x)=\sum_{i=0}^{N-1}S_i x^i
 =\sum_j \frac{e_j a_j}{1+a_jx}\pmod{x^N}.
\]

For distinct nonzero positions and nonzero magnitudes, its exact proper Padé
pair is `Omega/B`, with **monic** B, deg B=d<=N/2, deg Omega<d. Define a
derivation on syndrome coordinates by

\[
DS=S^2\pmod{x^N},\quad DS_{2i}=S_i^2,\quad DS_{2i+1}=0.
\]

The initialization of this tangent is binary linear. This says nothing yet
about the cost of propagating it through a reversible decoder.

The candidate tangent `DB=Omega, D Omega=0` satisfies

\[
D(\Omega/B)=\Omega^2/B^2=S^2\pmod{x^N}.
\]

It respects the fixed monic leading coefficient because deg Omega<d. The
Padé equations have a unique tangent: after reversing its columns, the d by d
Hankel matrix is

\[
 [S_{r+s+1}]_{r,s=0}^{d-1}
 =V\operatorname{diag}(e_j a_j^2)V^T.
\]

It is nonsingular by the Vandermonde determinant. Thus `DB=Omega` is an
identity on this degree chart, including cases where an ordinary EEA has
degree gaps. This proof does **not** differentiate a Boolean pivot dispatch or
drop nilpotent high coefficients on a singular Euclidean chart.

An equivalent root interpretation is `D a_j=e_j a_j^2`, `D e_j=0`.
Consequently `D beta_j=e_j`, where beta=1/a. Differentiating a monic locator
gives `DB(beta)=e B'(beta)`, the standard evaluator identity.

The 350 exact tangent checks solve the differentiated Hankel equations, rather
than invoking an existing decoder's derivative implementation.

## 2. Wronskian coordinates and quotient transport

In the compact convention, let A include the inverse normalizer so that

\[
 e(\beta)=\frac{\beta^{70}}{A(\beta)B'(\beta)}.
\]

Define `H^2=(AB)'=A B'+A' B`. For A=C^2+xD^2 and B=E^2+xO^2,

\[
H=CO+DE.
\]

At a locator root, H^2=A B', so `e=beta^70/H^2`. The Wronskian target was
already identified in round 4. The implementations here contribute a different
in-place Euclidean factorization and the completed 46T cubic point body.

If adjacent cofactors update by

\[
(A,B)\longmapsto(B,qB+\lambda A),
\]

then, for **any** polynomial quotient q,

\[
H_{new}=\sqrt{q'}\,B+\sqrt\lambda\,H.
\]

Here `sqrt(q')=sum sqrt(q_(2i+1)) x^i`; it is a polynomial. For degree-one or
degree-two quotients it is constant. For a degree-three quotient it can have
degree one: treating all degree gaps as constant q' would be incorrect.
With normalization `(Anew,Bnew)=(uB,v(qB+lambda A))`, multiply both terms by
sqrt(uv). Common padding P multiplies H by P.

`recurrence_components.py --kind transport` emits the linear/quadratic case
in place, including scale inversion and restoration. At k=3 it measures
407T forward and 814T for the phase sandwich, 158Q including 120 public bits.
Its phase checks the constant coefficient; the independent algebra test checks
the entire polynomial for general quotient degrees.

This transport is **not yet a decoder in (B,H) coordinates**: the next locator
update still needs A. Retaining A, B, and H costs at least 840+420 bits at
full size in this storage policy. Reconstructing A from (B,H) with the measured
Sylvester inverse before every update costs approximately 660,008T in one
direction across degrees 2..35, before the rest of decoding.

There is a precisely scoped algebraic obstruction on nonsquarefree B. For
example, B=x^2+1 has B'=0, and A=0 and A=1 produce the same H=0. The map
A->H therefore cannot be an unconditional replacement on those states without
retaining information. On squarefree B it is injective: H=0 implies
`A B'=0 mod B`, hence B divides A, hence A=0 when deg A<deg B. Equal finite
dimensions then give bijectivity. Intermediate decoder cofactors need not be
squarefree even when the final locator is. This is not a global decoder bound.

## 3. In-place Sylvester factorization

For deg f=m>=deg g=n, let c have n coefficients and d have m coefficients.
The target map is `h=c*f+d*g`, with exactly m+n output coefficients in the
same c/d slots. Divide f=qg+r, where deg r=k<n. Then

\[
h=c r+(d+qc)g.
\]

Set d'=d+qc and split `d'=d_low+x^k d_high`. The first n+k output
coefficients are obtained by the smaller Sylvester map

\[
(d_{low},c;g,r)\longmapsto d_{low}g+cr,
\]

followed by the low part of `x^k d_high*g`. The remaining m-k coefficients
are an upper triangular transform of d_high with diagonal lead(g). The
diagonal is nonzero. Each operation is a reversible shear or a charged
in-place nonzero field scaling.

Polynomial division stores its quotient in the original high f slots and its
remainder in the low slots. The source division is reversed after the payload
transform. Sources f/g are restored on every recursion return. Inverse field
words are released before recursion; only one zero-predicate bit per level
remains live. Recomputing the predicate after restoring the source avoids
holding all conjunction ladders simultaneously.

The classical reference permits arbitrary k. The emitted circuit currently
uses k=n-1. It arms zero leading coefficients to one, defining a reversible
permutation even on unsupported charts. On those charts it **does not** claim
to compute the desired Sylvester map. Native tests exercise both composition
orders with zero remainder chains; those are restoration checks, not evidence
that the missing gap chart is solved.

For odd B degree 2t+1, O is monic degree t. Separate the top C coefficient:
apply Sylvester(E,O) to D and the lower C coefficients, then add the top-C
cross term. For even B degree 2t, E is monic degree t and O has declared
degree t-1. Separate the top D coefficient instead. Native degrees 6 and 7
check both layouts. Degree 35 is native checked too.

For t>=1, measured generic Sylvester forward cost is
`126*t^2+790*t-518`. Including the odd H adapter gives
`126*t^2+832*t-518`; the even-degree H cost is
`126*t^2+664*t-582`. These formulas count the implemented elementary field
products, inverses, zero tests, and source restoration, and reproduce the
native examples. Degree-35 H costs 50,040T forward, 100,080T build/erase,
897Q with all 840 public bits. The native nine-point phase sandwich has 918Q.

## 4. Runtime-mask phase without a quantum product

For a fixed public scan point beta and a runtime coordinate mask M, define

\[
 M'_i=\operatorname{parity}\big(M\mathbin{\&}
                         \beta^{70}(2^i)^2\big).
\]

This is an invertible binary linear map of the **runtime classical** mask.
Then `parity(M & beta^70/H^2)=parity(M' & inverse(H))`. A classical trace-dual
conversion/scaling allows the existing cubic absolute-trace inverse circuit
to supply this phase. The implemented body costs 46T, 66Q with 24 public
quantum bits. It restores H, the root accumulator, its mask, and its scratch.
The nine-point native tests include beta!=1 and the beta^70 factor; no
syndrome-dependent preprocessing is hidden in a classical constant.

## 5. Fixed-degree squarefree padding, with charged adapters

The final common x^m padding is not squarefree when m>1. Replace it with

\[
P=x^{m\bmod2}\prod_{k\text{ odd},\ d\leq k<35}
       (x^2+x+c_k).
\]

The c_k are the existing distinct irreducible-quadratic constants from the
completed round-4 reference, converted through the field isomorphism for the
independent model. Their absolute trace is one. P has no physical nonzero
field roots, and is coprime to the split true locator. Thus P*Btrue is
squarefree of fixed degree 35, including variable true degree d.

`padding_adapter.py` implements replacement of paired x factors using
controlled binary transvections, with an explicit quantum A/H leading
coefficient. A pair adapter costs 7,344T per direction, and an H-only adapter
costs 3,672T. All flags, including non-unary flags, are quantum inputs in the
native tests. No branch frequency discounts their static gate cost.

A straightforward charged assembly is:

1. Decode to common-monomial A/B and correct the previous-cofactor frame.
2. Pair repad A/B to squarefree P padding.
3. Transform A to H, restoring the padded B.
4. Inverse pair repad H/B back to common-monomial padding.
5. Remove common monomials by the existing controlled rotation method.
6. Apply the 46T point scan; reverse all five preparations and the decoder.

Four pair adapters cost 29,376T over the lifetime. The model also charges the
existing entire 10,790T unpadding/preparation budget. Combined with the
824,740T incumbent core, H, and scan, this estimates 1,153,356T at a target
width of 1004. Integration of the routing and gap-core flags remains unbuilt;
1004 is a feasibility assumption, not an emitted full-circuit width.

The independent completed-pipeline model includes the actual armed Sylvester
map and the literal inverse padding/rotation network. On its two fixed
20,000-shot host-law streams it has 81/78 wrong magnitude vectors and 41/40
random-mask phase errors. Native primitives restoring individually do not
turn that model into native full-oracle acceptance.

## 6. A different decoder: C-fraction coefficients in syndrome storage

For a formal series f with nonzero f0, one C-fraction peel is

\[
a=f(0),\quad g=(a/f+1)/x,\quad f=a/(1+xg).
\]

Two peels can be emitted particularly simply. Let a=f0, u=f1, b=u/a,
T=(f-a)/x. Then

\[
g=((u/T)+1)/x+b,\qquad
f=\frac{a}{1+b x/(1+xg)}.
\]

Normalize the existing tail by u, apply the read-only optimized unit-series
reciprocal **in that tail**, normalize u by a, and XOR b into the first tail
coefficient. The first two words now hold a,b; the remaining words hold g.
The unit reciprocal ignores its implicit constant slot, which can retain b.
The inverse performs these operations in reverse order. No separate original
syndrome or locator/template bank is retained.

Given a rational tail g=P/Q with Q0=1, folding the heads gives

\[
\Omega=a(Q+xP),\qquad B=(1+b x)Q+xP.
\]

This too is in place. First replace each P_(i-1) by R_i=P_(i-1)+Q_i.
Multiply Q_i by b and add R_(i+1), producing B_(i+1). Replace b by b+R_1,
producing B_1; multiply the R words by a, producing Omega_i. The retained a
is Omega_0. All scalings have explicitly computed inverses, and all inverse
scratch is erased. Recursing over pairs yields the free locator coefficients
in the odd slots and evaluator coefficients in the even slots.

The recurrence also clarifies the Riccati tangent. If `Df=u*f+v*f^2`, a
single peel gives `Da=u*a+v*a^2`, `Dg=v*a*g`. Starting with u=0,v=1, only
the first two C-fraction heads have nonzero tangents:

\[
Da_0=a_0^2,\quad Da_1=a_0a_1,\quad Da_i=0\ (i\geq2).
\]

After monic normalization, this is the same `DB=Omega, D Omega=0` identity.
Forty independent dual-series checks test this localization directly.

Zero tail after completion is benign: all later heads and data remain zero
under the totalized scale operations. A zero head with a nonzero tail needs
another continued-fraction chart. Replacing it by one is a reversible
implementation convention, not exact decoding. Independent host-law models
find approximately 1.4% such wrong pairs.

The current serial extraction costs 503,655T each way; synthesis costs 57,664T
each way. Core peak is 878Q. The phase uses the standard 80T evaluator division
body. B1 can hold B'(beta), B2 can hold B(beta), and Omega0 can hold
Omega(beta), through an invertible linear coordinate change for beta!=0.
This removes three accumulator words: predicted full width 884Q, cross-checked
by the native N6 and N14 complete small oracles. Total full-size gate estimate
is 1,450,238T. No full-size decoder binary was built or scored from that model.

## 7. Direct moment functional

Convert each runtime mask to its field trace-dual mu_j. Interpolate a classical
polynomial `h(a)=sum_(k=1)^4095 h_k a^k` with h(a_j)=mu_j at every physical
position. This basis is invertible on the nonzero field elements. The required
phase is

\[
\operatorname{Tr}\sum_{k=1}^{4095}h_k S_{k-1}.
\]

For constant-normalized locator B=1+b1*x+...+bd*x^d, moments satisfy
`S_r=sum_(i=1)^d b_i S_(r-i)`. A rolling bank can overwrite its oldest word
by multiplying it by the nonzero leading bd and adding the other d-1 products.
`recurrence_components.py --kind moment` measures 239T forward for d=3 and
1,583T forward for d=35, including one shared inverse and a zero test.
Runtime-mask phase, inverse, source and scratch restoration pass natively.

Even granting 70 initial moments and sharing every pivot inverse, the serial
implementation needs at least `2*4025*35*42=11,833,500T` for the extra moments
and their reverse. This is a concrete rejection of this rolling implementation,
not a lower bound on all transposed or direct phase algorithms. Interpolated
mask coefficients depend only on runtime masks, never on quantum syndromes.

## What remains unresolved

The missing high-value construction is a low-width blocked extraction or
different recurrence that retains the C-fraction/tangent information without
repeating full tail reciprocals, and handles degree gaps cheaply. At Q884, a
1.5x result with the current 80T scan/synthesis would require forward
extraction <=292,699T, versus 503,655T now. A 2x result requires <=164,158T.
Those budgets assume the zero-error confidence floor; gap repair has to fit
inside them. Neither abstract bijectivity nor the tangent identity provides
that circuit automatically.
