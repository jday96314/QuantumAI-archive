# Candidate DQI scan: BM-template Wronskian

Follow-up cost correction: the current selected inverse-product kernel is88T,
while an inverse-only adaptation of its cubic trace formula would use46T.
Furthermore, `profile_direct_repad_transvections.py` shows that applying the
existing controlled-transvection helper reduces common-to-squarefree cofactor
padding from 28069T to 7140T per direction at N70. The completed padding adapter
in `transvection_repad.py` passes 4096 cases each forward, inverse and both fresh
composition orders with arbitrary coefficients and flags, clean scratch and
correct phase. Its compiled width is 875 public qubits; the resource-only
profile included three unused scratch qubits. This invalidates the older
padding-cost premise for rejecting that route. The inverse-only cubic kernel
is not yet built, H division remains to be integrated, and the generic Sylvester
map still excludes degree gaps. A full decoder resource or correctness
improvement remains unproven. Evidence:
`results/transvection_repad_n70_check.json` and
`results/direct_repad_transvection_profile.json`.

This research is **not in the submission**. The current packaged score remains
9.367722153134604. The quantum inverse-mask primitive is implemented; the full
reversible template-to-Wronskian transform is not.

The earlier scan discussed below used 80 Toffolis per field position. The primitive in
`inverse_mask_kernel.py` applies the phase given by the parity of a runtime
classical mask against a field inverse, gated by a twelve-bit root test. It uses
**50 Toffolis**, with 56 total qubits including 24 quantum inputs. Its test covers
all 4,096 divisors, root/non-root inputs, random masks, input restoration,
zero scratch, and phase. See `results/inverse_mask_kernel.json`.

Let sigma be the final q-form BM locator and tau its template. In the ordinary
(non-streaming) classical BM model, the final template is cyclically rotated.
Let d be the final signed counter, and set d_eff = d except that d=0 maps to -1.
Undo d_eff cyclic positions in the 37-word template. For tested valid errors of
weight L, the resulting template has tau[0]=0 and degree L+1, and at each actual
error root x,

    error_value(x) * sigma'(x) * tau(x) = x^(2L).

Define

    H(x)^2 = (sigma(x) * tau(x))'
           = sigma'(x) * tau(x) + sigma(x) * tau'(x).

In characteristic two only even powers survive this derivative, so H is a
polynomial of degree at most L. At an error root, this gives

    error_value(x) = 1 / (H(x) * x^(-L))^2.

A 36-word barrel shift by 35-L can store H*x^(35-L); evaluation with fixed
powers x^(i-35) yields H*x^(-L). Squaring is a free linear map on field bits, so
the classical mask can absorb the final square when applying inverse parity.
This would remove the independent quantum omega operand from each scan point.

`probe_wronskian.py` checks the exponent/rotation convention on seven weights.
`validate_wronskian.py` checks the full polynomial identity for arbitrary sampled
positions and nonzero values. The first 1,000 shots, spanning all weights 0..35,
pass. A larger run probes rare degree/window cases; this remains a classical
model check, not the host DQI verifier and not a submission score.

The remaining implementation problems are substantial:

* Compute the Wronskian **in place** in the BM-template bank. A fresh 432-qubit
  output would erase the gate savings because the DQI score squares width.
* For a fixed locator, this is a structured linear map on the template. Writing
  sigma=E^2+x O^2 and tau=U^2+x V^2 gives H=EV+OU. The operator resembles a
  Sylvester map. Generic quantum Gaussian elimination is too expensive; an
  in-place Euclidean factorization or comparable structured method is needed.
* In the syndrome-consuming decoder, inactive paths may leave original syndrome
  values in unused template slots. The current identity tests do **not** model
  that borrowed storage. After computing omega from the retained tail, its
  coefficients above L are zero on a correct decode. These slots might park
  the unused template garbage: template index k>=L+2 pairs with omega index k-2.
  This possible layout needs an exact classical streaming model first.
* Degree/window exceptions and zero pivots must restore all inputs. Existing
  accepted phase errors may be priced, but syndrome corruption is unacceptable.

A direct out-of-place Toom product for the Wronskian itself needs only about two
36-point field-product evaluations. Storage and in-place invertibility, rather
than that arithmetic count, are the main unresolved issues. No whole-decoder
resource or correctness claim follows from the 50-Toffoli primitive.

The inherited public decoder already uses a different reciprocal-polynomial
method, described in its [primary documentation](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html).
It stores a Bezout pair in a zero syndrome tail; that tail is occupied by the
current streaming decoder's backup bank. The [Sylvester-map description](https://www.math.ias.edu/~avi/PUBLICATIONS/AndrewsWi_Constant-Depth%20Arithmetic_2024.pdf)
provides background for the structured linear-map formulation. The Wronskian
layout above is a local proposed adaptation supported by the listed probes.

## Larger checks and actual streaming layout

The 100,000-shot unconsumed model test found 28 failures. It must not be
reported as a universal identity proof. The corrected streaming model performs
better: its 100,000 uniform-weight samples found six failures, all of which
already had an incorrect locator. No additional Wronskian failures appeared
among its 99,994 correct locators. These are diagnostic model results, not the
host's weighted default verification.

The actual streaming bank uses growing-prefix rotations, so the simple cyclic
counter correction described above is insufficient for it. The function
`probe_stream_template.normalize_bank` reverses those rotations after completion.
It depends only on L. The new quantum `template_permutation.py` realizes that
37-word permutation in 11 layers with **2,738 Toffolis and 24 extra qubits**.
Its quantum test passes; it fits within the decoder's existing 49 free wires.

`capture_template.py` records the template and counter wire locations without
changing the existing BM tape. All old hook fields match experiment 3 exactly.
`check_template_prefix.py` compares the actual compiled BM prefix with the
classical streaming model on 1,024 inputs, including all 70 syndrome words,
36 locator words, 37 template words, and the counter. These value checks pass.
They do not assert phase correctness of an unfinished full decoder.

The models' reported counter is **one larger** than the native final counter:
the model retains the otherwise omitted round-68 increment. For native dW,
L = floor((69-dW)/2) on correct paths. This distinction is verified by the
compiled-prefix comparison.

The degree-trimmed template check does not yet implement reversible garbage
parking. Nor is the in-place Sylvester/Wronskian transform implemented. Those
remain the next substantive tasks; adding a fresh full polynomial bank would
not establish a useful decoder improvement.

## Subsequent generic transform experiment

An in-place Sylvester transform has now been implemented and quantum-tested
on512 generic degree-17 pairs with E[0]=1 and no Euclidean degree gaps. It uses
46 extra qubits and62,098 Toffolis per direction. Two directions already cost
slightly more than the theoretical30*4095 scan saving, before layout overhead.
It is not a complete decoder or a score improvement. See SERIES_RESEARCH.md
and results/sylvester_quantum_check.json for scope and evidence. Handling degree
gaps, garbage parking, and a substantially cheaper transform remain unresolved.
