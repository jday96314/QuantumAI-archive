# Circuit attribution

This archive combines unchanged, independently verified components from the
round-six agent14 and agent20 workspaces. It does not claim third-party
components as original work.

The inverse-square-root circuit is Phuc Thai's EXP225a, SHA256
`0f84f5bab3a02e6b357d4656553607174345a33d252cef643f9136db13212141`,
reproduced from the saved public notebook output:
https://www.kaggle.com/code/phucthaiv02/qcc-verified-portfolio-df-thc-retraining-study

Phuc Thai credits Dmitrii Gluzdov's 84-qubit construction, built on Sumit
Chongder, anazemcev, Kozykorpesh and Zongao Bian's arithmetic, and marsyas's
fused phase-stage idea. Phuc Thai adds carry-save rows, a two-bit counter,
parity refitting, Cuccaro prefixes, normalization changes, exact rewrites
and phase fusion. The accompanying Apache-2.0 LICENSE is retained from the
agent20 public-component package for the credited code lineage.

The DQI circuit is agent14's raw continued-fraction quadratic-extension
Schur decoder with amplitude-polynomial scanning, a three-field terminal
state, and remainder-consuming quotient recovery. Its SHA256 is
`9e36454e027e14963dee25ed7b2595b2ab752df3d6982dd231ca32d93eaf8fd6`.
It reuses the repository's finite-field arithmetic and compiler lineage,
including andigu/uncompute, Apache-2.0, pinned revision
`64b5a13aa249f1360846657392e6c66817b26d0f`, retained under
`initial_experiments/3/vendor/uncompute` in the research repository.
This archive does not contain Phuc Thai's EXP119a decoder.

Chemistry retains the round-five structured banked physical factorization.
Agent20's `paired_chemistry_angles.py` adapts the paired bilinear combine
credited to Phuc Thai's fouraddr2 compiler to the existing compressed
eight-bit bank. Two low-address terms share one Toffoli; the known bank
product is compiled into the base table. Physical NPZ bytes and loaded
angle words are unchanged from agent14's selected chemistry. The new
ANGLE circuits passed every valid address with arbitrary incoming output
words and both full native verification runs.

All three components have existing full default and independent native
verification reports. Integration reuses that evidence after checking
circuit hashes and the canonical chemistry digest before packaging and
after extraction. The DQI decoder retains rare unsupported extension-pole
cases: each saved full run has five value errors and zero phase errors in
102,400 shots, charged by the unchanged verifier's repetition penalty.
