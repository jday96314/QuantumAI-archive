# Circuit attribution

The inverse-square-root circuit is Phuc Thai's EXP222a, SHA256
58c7ccb60f61088c9c3d14a562f6c1d795bcb8ac12b364cb16941b54b74f5cf1,
reproduced without changes from Ahmed Berat Özer's public notebook output:
https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression
Phuc Thai's work builds on Dmitrii Gluzdov's 84-qubit arithmetic and the
Sumit Chongder, anazemcev, Kozykorpesh and Zongao Bian lineage.

The decoder is Ahmed Berat Özer's merged J-fraction circuit, SHA256
a91c66ac9d0a720e112976ae13a0c4b855223c6605c42b0809de2ea7e6ad5a25,
from the same public notebook, reproduced without changes. Its original
NOTICE is reproduced below. The accompanying Apache-2.0 LICENSE is retained.

Chemistry uses the existing round-five structured banked factorization and
ANGLE compiler. This local submission combines independently verified
components; it does not claim the public circuits as original work.
The unchanged organizer verifiers were rerun locally at full default shots
and an independently fixed salt for both new components.

Original decoder NOTICE:

# vqcc_dqi_src

SPDX-License-Identifier: Apache-2.0 (see LICENSE).

* SomeJoint / Agniv Banerjee: the in-place J-fraction decoder, its finite-field
  arithmetic, the two-slot 3+1 gap repair, the conversion and the Forney scan
  (public notebook "Quantum Circuits | 8.95 LB: Cost Breakdown", Apache-2.0).
  All files except the three below are SomeJoint's, unchanged.
* manwithacat: the classical mask projection (`mask_projection.py`, as
  republished by anazemcev, Apache-2.0).
* anazemcev (Artem1410): the ideas of reusing the HMR-cleared word as the AND
  mask, sharing one word between sigma(x) and the Forney inverse, discarding
  the inverse scaffold, and the 4-to-3 gap-flag codec ("Compact DQI and
  Carry-Save AR", Apache-2.0). Re-implemented in `jfraction_merged.py` and
  `jfraction_merged_decoder.py`.
* New in `jfraction_merged.py`: the merged two-slot extraction (one unit-series
  reciprocal per two slots, five gap modes) and the scale/premultiply lifetime
  changes.
