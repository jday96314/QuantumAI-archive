# Fixed-boundary dispatch width plan (before emission)

Status: allocation contract to validate, not a measured full decoder width.
At initial planning dqi-core's full builder was not accepted, so only its two
checked components were adopted. Its report later changed to ACCEPTED; the
verified trimmed builder/KMB are now the read-only base reference. The protected
round-4 asset remains intact. No exception-correct replacement was accepted here.

The physical polynomial state stays in 70 field words (840 qubits). At nominal
boundary k the fixed remainder tail has exactly 70-2k words and pinv occupies
its fixed word. The cofactor prefix has 2k-1 words. A branch changes only how
that prefix encodes the two cofactors and its zero exceptional coordinates;
it does not shift the remainder tail. These abstract word positions map to the
incumbent frames' fixed physical handles, rather than contiguous syndrome order.

Persistent metadata is six completion bits, two six-bit gap-stage tags (zero
means unused), and a three-bit mode. A cubic quotient consumes two gap tags;
more than two missing degree slots is an explicit model failure, counted in
all censuses. A temporary selector contributes one further bit.

| Mode | Pending slots | Previous quotient degree | Current pending degree |
| --- | ---: | ---: | ---: |
| 0 | 0 | 1 | none |
| 1 | 0 | 2 | none |
| 2 | 0 | 3 | none |
| 3 | 1 | 1 | 2 |
| 4 | 1 | 2 | 2 |
| 5 | 1 | 1 | 3 |
| 6 | 2 | 1 | 3 |
| 7 | unsupported | failure extension | not accepted |

A one-slot quadratic deferral stores beta*Bprev, B, inverse(p*beta), and a
remainder band shortened by two words. The omitted high words are (0,r) with
r=inverse(pinv). A two-slot cubic deferral stores
A=beta*Bprev+q2*x^2*B and pinv/beta, and omits (0,0,r,s), with
s=(q2+B_sublead)*r. This encodes its partial arithmetic in the public state;
there is no matrix, extra moment history bank, or Bennett bridge. One-slot
cubic deferral omits two known zeros. The model packs and unpacks solely from
70 public words plus the metadata, and asserts the reconstructed EEA state.

At a block with pending carry the logical moment prefix can have ten words,
but at most six are fresh: the carry's omitted words are reconstructed from
the public partial cofactor and pinv. Only six fresh words may be copied.

| Stage | State | Metadata/control | Fresh head | Arithmetic scratch | Planned peak |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial k=0 slot | 840 | 22 | 0 | <=39 | <=901 |
| Block k=1 | 840 | 22 | 72 | 39 | **973** |
| Block k=4 | 840 | 22 | 72 | 39 | **973** |
| Block k=7 | 840 | 22 | 72 | 39 | **973** |
| Block k=10 | 840 | 22 | 72 | 39 | **973** |
| Block k=13 | 840 | 22 | 72 | 39 | **973** |
| Block k=16 | 840 | 22 | 72 | 39 | **973** |
| Block k=19 | 840 | 22 | 72 | 39 | **973** |
| Serial k=22..34 | 840 | 22 | <=24 | <=51 | <=937 |
| Scan, explicit leading B | 840 | 22 | 0 | 12+36+42 | <=952 |

The 39 arithmetic bits budget the ordinary 27 (including three emitter bits)
plus one reusable 12-bit coefficient/cache word. Conditional head arithmetic,
masked quotient coefficients, and reconstructed pending pivots must reuse this
word in disjoint lifetimes. They may not each allocate another twelve bits.
The inverse policy must account for these live bits; new contexts exceeding
the component's checked 12-bit body allowance must use the low-workspace path.

The unit contract is now asserted inside blocks as well as at their boundaries:
a quadratic slot first scales only, then performs its two coefficient MACs;
a cubic slot first does nothing, next scales and adds q2*x^2*B, then finishes
with q1*x*B+q0*B. Sparse history records both known missing slots as soon as a
cubic is identified, including a future slot. Mode is asserted against a function
of that history and the completion count; it is not additional hidden metadata.

Native arithmetic embeddings reach 973 Q for the Schur head and 960 Q for
the cofactor, with operands retained. These do not validate complete pending-prefix
placement. In particular, a two-slot cubic has two nonzero omitted moments;
one must use a public carry coefficient and one the reusable cache. A live cache
overlap would add twelve bits (985 rather than 973). No full integration or
native width acceptance is claimed. The realization was stopped on its separate
moment-chart T projection, above the revised +300k limit.

No acceptance claim follows from this table. Any extra live word invalidates
the 973-Q plan; a measured peak above 990 triggers the new kill condition.
The target <=965 needs at least eight further bits of saving; <=950 is stretch.
Routing, field arithmetic, head selection and cleanup are charged on all inputs.
