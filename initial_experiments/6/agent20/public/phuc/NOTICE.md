# Circuit attribution

The inverse-square-root circuit is Phuc Thai's EXP225a, SHA256
0f84f5bab3a02e6b357d4656553607174345a33d252cef643f9136db13212141.
The decoder is Phuc Thai's EXP119a, SHA256
c5da9cd4e2d73dd41417e539f2659009d646b3ee30ddb5ca494261378a502856.
Both are reproduced without changes from the public notebook output:
https://www.kaggle.com/code/phucthaiv02/qcc-verified-portfolio-df-thc-retraining-study
The accompanying Apache-2.0 LICENSE is retained for the credited code lineage.

Phuc Thai credits Dmitrii Gluzdov's 84-qubit inverse-square-root construction,
built on Sumit Chongder, anazemcev, Kozykorpesh and Zongao Bian's arithmetic,
and marsyas's fused phase-stage idea. Phuc Thai adds carry-save rows, a
two-bit counter, parity refitting, Cuccaro prefixes, normalization changes,
exact rewrites and phase fusion.

The decoder builds on SomeJoint / Agniv Banerjee's in-place J-fraction
Reed–Solomon decoder and finite-field arithmetic, with manwithacat's
classical-mask projection. Phuc Thai adds scratch lifetime changes, a lean
Forney scan, two gap-position pointers and paired extraction through a shared
series reciprocal. The underlying public code is credited as Apache-2.0 in
the notebook. See the downloaded notebook snapshot and release manifest in
the research workspace for the full provenance.

Chemistry uses the existing round-five structured banked factorization and
ANGLE compiler. This local submission combines independently verified
components; it does not claim the public circuits as original work.
The unchanged organizer verifiers were rerun locally at full default shots
and an independently fixed salt for both new components.
