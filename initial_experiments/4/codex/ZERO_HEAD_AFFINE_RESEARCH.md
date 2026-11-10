# Fixed top bit and measured phase headroom

The fully native-verified and independently audited archive now scores
**9.367722153134604**. Experiments are stopped at the user's request; the 14.0
target remains unmet. See [EXPERIMENT_REPORT.md](EXPERIMENT_REPORT.md). ZIP SHA256:
`9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68`.
The previous9.3326564887archive and both audits are preserved under
`revisions/cubic_trace_product/`. DQI and chemistry are unchanged.

The new candidate is `candidates/zero_head_affine.kmb`, independently rebuilt
byte-for-byte as `candidates/zero_head_affine_rebuilt.kmb`. SHA256:
`27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1`.
Its source manifest is `parameters/zero_head_affine_source_manifest.json`.
The coefficient table, numerical approximation,90-qubit workspace schedule
and34 rotations are unchanged. Static Toffolis fall1816->1793.

Two changes account for the reduction. First, the selected integer model's
normalized mantissa is always in [67108897,134217726], strictly within
[2^26,2^27). Its bit26 is one before and after affine restoration. Thus the
offset can be restored modulo2^26 while leaving the borrowed normalized input
bit untouched. All3133452 values were checked in
`results/affine_restore_top_bit_bound.json`.

The constant3*2^23 modulo2^26 uses a one-AND three-bit increment: complement
the low bit a, XOR it into the next bit c, then flip the high bit unless a*c.
The AND is erased by measurement against the unchanged low bits. The parity
restoration becomes a CNOT and the coarse-address subtraction is one bit
shorter. These changes save12T in both directions combined.

Second, output barrel-shift moves into initially clean private headroom use

    t=c*s; s^=t,

and reverse via

    s^=t; erase_AND(c,s,t).

The reverse AND uses measurement and a Clifford phase correction, saving one
Toffoli per clean destination. There are ten such destinations across the
four stages. The first top-bit move starts from the known one, so two CNOTs
replace its forward AND as well. Normalization still uses its original swaps;
private phase headroom is the only storage measured by this change.

The first three earlier variations remain separate evidence:
`constant3_affine` passes all200540928 native shots at184666.4402837133 cost,
90Q/1812staticT. `fixed_top_affine` passes the same full test at
183946.526175108 cost,90Q/1804staticT. Neither is the final combined candidate.
The combined `zero_head_affine` passes a one-shot native superset test over
x=1..3133696, including the full requested domain plus244 values, with
unchanged worst error/tolerance0.9796714395489131 at958102. This one-shot
screen is not a substitute for its separate default64-shot verification.

That full verification is now complete:200540928shots PASS,90Q,
1782.000090874218meanT,1793staticT,34rotations,
cost182956.5271606425. The worst error/tolerance stays0.9796714395489131.
Evidence: `results/zero_head_affine_verify.json`. Both archive audits pass,
including independent ZIP reopening, all14members, canonical chemistry digest,
and binding to the full native default shot counts for every task.

Rebuild and run the full required verification:

```bash
source env.sh
"$EXP_PYTHON" zero_head_affine_invsqrt.py
"$EXP_PYTHON" ../../../from_host/verifiers/invsqrt/verify_submission_invsqrt.py \
  candidates/zero_head_affine.kmb --max-workers 4
```
