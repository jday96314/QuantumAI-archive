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
