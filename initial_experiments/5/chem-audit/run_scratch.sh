#!/usr/bin/env bash
set -euo pipefail
source ./env.sh
# One GPU job at a time, at most two CPU threads. All seeds are explicit.
for shape in '1 639 350 12' '2 447 175 13' '2 463 175 14' '2 479 175 15' '4 180 88 16'; do
    read -r rank bases columns seed <<< "$shape"
    tag="scratch_r${rank}b${bases}_l300"
    if "$AUDIT_PY" -c 'import json,sys; from pathlib import Path; p=Path(sys.argv[1]); sys.exit(0 if p.exists() and json.loads(p.read_text()).get("finished") else 1)' "results/${tag}.json"; then continue; fi
    if [ "$rank" = 1 ]; then norm_cap=5; else norm_cap=100; fi
    "$AUDIT_PY" frontier_fit.py --R "$rank" --B "$bases" --C "$columns" --seed "$seed" \
        --lambda-reg 300 --steps 20000 --checkpoint 2000 --precision float32 --weight-cap "$norm_cap" --tag "$tag" > "results/${tag}.log" 2>&1
done

# Sweep normalization pressure from each scratch-trained basin. These are
# continuations, explicitly recorded as such, rather than independent starts.
for tag in scratch_r1b719 scratch_r1b639 scratch_r2b447 scratch_r2b463 scratch_r2b479 scratch_r4b180; do
    AUDIT_DEVICE=cpu "$AUDIT_PY" gauge.py --fit "parameters/${tag}_l300_last.npz" \
        --out "parameters/${tag}_gauged.npz" --sweeps 20 > "results/${tag}_gauged.log" 2>&1
    for regularization in 30 3000; do
        name="${tag}_l${regularization}"
        "$AUDIT_PY" frontier_fit.py --init "parameters/${tag}_gauged.npz" --seed 23 \
            --lambda-reg "$regularization" --steps 6000 --lr .0005 --final-lr .000001 \
            --checkpoint 2000 --precision float32 --weight-cap 100 --tag "$name" > "results/${name}.log" 2>&1
    done
done
