# Follow-up: emit and bank the trimmed decoder

Bank your 25,596-T saving as a complete circuit:
1. Emit the full selected decoder with only the trimmed middle products and budgeted inverse retention (no schedule change). Compact wires and confirm Q=961 and static T=994,358.
2. Run the incumbent and the trimmed circuit through the full verifier with the **same explicit `--seed`**. The decoding logic is unchanged, so failure counts must match exactly; any difference is a bug.
3. Run the trimmed circuit at default settings and with one independent salt (102,400 shots each, under the flock). Report R, mean T and cost for both. Accept on lower T at equal Q with matched decoding; do not select on failure-count luck.
4. Write the accepted builder, KMB and verifier JSONs here and put "ACCEPTED: trimmed decoder" at the top of RESULTS.md so dqi-dispatch can build on it. Then stop.
CPU cap 2.