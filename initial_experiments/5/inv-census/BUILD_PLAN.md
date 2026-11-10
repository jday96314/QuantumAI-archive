# Phase 3 — plans and qualification gates

No complete compiled design in this census projects to the 1.3× cutoff
140735.7901. None warrants a full verifier run or selection change.

Two Phase 1 hypotheses crossed that cutoff under explicit assumptions. The
exponent-Walsh hypothesis is now disqualified: it requires 307 rotations rather
than <=35. The fused square hypothesis remains a conditional development plan;
its numerical uncertainty has narrowed, but its T/Q budget is unproved.

## Fused square / middle product: conditional plan

The exact all-input model now establishes a useful smaller arithmetic task:

1. Normalize as before, with k and s=floor(k/2).
2. Set signed z=(M mod2^14)-8192. Use the original full-precision z in the
   linear term. Only the quadratic term gets rounded.
3. Let u=floor((z+16)/32), and
   h=round((32*u)^2/2^16). Thus h=round(u^2/64).
4. Accumulate the phase numerator
   `c0*2^23 - c1*z*2^11 + c2*h*2^16`, then scale by2^(s+4) in2^-64 turns.

All 3,133,452 inputs pass: worst error/tolerance0.8046964363863656 at x2424816;
margin 137902011 units=7.475682995815996e-12 turns. This uses the same coefficient
table without any new coefficient fitting. A six-bit coordinate rounding fails
(ratio 1.1173053262473323); the tested five-bit rounding is not a generic bound.

**Rounding guard:** u includes both -256 and+256. A signed 9-bit register is
insufficient. Use unsigned abs(u), which needs 9 bits including endpoint256,
or signed 10 bits. The significant h word needs 11 bits. The full signed-coordinate
and square sizes are recorded in results/curvature_frontier.json.

One's-complement folding helps: for z<0, v=-z-1; otherwise v=z. Then
`abs(floor((z+16)/32)) = floor((v+16)/32)` in both branches. This permits
Clifford sign folding of the low13 coordinate bits and one rounded increment;
keep the rounding guard until the coordinate is reversed.

### Concrete implementation sequence

- Build a standalone rounded-absolute-coordinate kernel, exhaust every 14-bit z
  and repeated native measurement seeds. Check the endpoint guard explicitly.
- Build `h=round(abs(u)^2/64)` directly into11 bits. A full square is17 bits,
  not the prior28-bit square. Start with an exact square kernel as a reference;
  replace it with truncated-column accumulation only after its carries are
  checked on all 257 unsigned magnitudes. A sum of individually truncated
  monomials is not automatically the rounded high part of the full square.
- Implement the7-by11 curvature middle product directly into the low constant
  accumulator with signed overflow/carry repair. Exhaust all 128 curvature words
  and 257 magnitudes. Do not materialize a wider c2*square result.
- Combine with the full14-bit linear product and a single shared phase stage.
  Refit table constants only if the actual omitted-column pipeline requires it;
  rerun the full integer domain for every numerical change.
- Measure each table erasure in X, retain outcomes, and repair its Boolean phase
  using the original address after restoring any coordinate permutation. Reuse
  the erased bank for carry workspace. Erasure is not free without phase repair.
- Compile complete liveness and count both forward and reverse work, coefficient
  replay, normalization, exponent application and rotation synthesis.

### Resource qualification, not an achieved projection

For the proved margin, a hypothetical Q 80,T 1400,R 50 circuit costs about 140.2k
and barely qualifies at 1.3×. Q 80,T 1400,R 34 costs about 131k. Achieving1.5× atQ 70
andR 34 requires mean Toffolis <= about 1505. These are engineering budgets.

The difficult live stage is the curvature update: input+flags27, c0 remainder23,
curvature7, square11, rounded-coordinate guard/workspace, and possibly slope17.
Keeping the slope resident can exceed80 before the adder is allocated. Loading
it later shortens the live set but can add a second large lookup. The previous
fused-square132k–140k estimate is therefore **conditional**, not a defensible
complete circuit projection at this point. Likewise, phasing curvature separately
can duplicate exponent routing and rotations, exactly the failure exposed by
inv-phase; it must not be omitted from the budget.

If the first complete draft loses, diagnose whether the cost is the exact square,
rounded-coordinate guard, table replay, carry repair, or repeated exponent phase.
Develop the corresponding primitive before discarding the architecture. Current
results leave this architecture open; they do not establish its resource target.

Only a complete candidate with an exhaustive integer margin and sampled native
model matches should proceed to full verification. Use the unmodified default
200540928-shot host workload under `flock /tmp/r6-fullverify.lock`, at most3 CPU
workers, before replacing the selection. Do not tune hashes or salts.

## Other architectures

The packed-exponent map is injective on every original input, but is not an
operand-ready normalization: its marker must be removed or absorbed into the
fit, and decode work must be counted. It is a possible later complement to a
new arithmetic architecture, not a five-qubit saving already obtained here.

Coefficient-scaled output, cubic Horner, compressed line lookup, and truncated
Newton remain the Phase 1 research alternatives. None currently has a supported
complete1.3× projection, so no acceptance/build commitment is claimed for them.
