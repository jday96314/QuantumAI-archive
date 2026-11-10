#!/usr/bin/env bash
set -euo pipefail
source ./env.sh
for config in '31 2 .0002' '32 4 .00005'; do
    read -r seed samples rate <<< "$config"
    tag="robust_banked_s${seed}"
    "$AUDIT_PY" frontier_fit.py --init candidates/banked/submission.npz --seed "$seed" \
        --lambda-reg 300 --steps 8000 --lr "$rate" --final-lr .000001 \
        --noise "$samples" --weight-cap 5 --lambda-cap 315 --enforce-lambda \
        --precision float32 --checkpoint 1000 --tag "$tag" > "results/${tag}.log" 2>&1
    "$AUDIT_PY" integer_refine.py --init "parameters/${tag}_rounded.npz" \
        --out "parameters/${tag}_bounded_integer.npz" --steps 40 --max-rows 64 --cycles 2 \
        > "results/${tag}_bounded_integer.log" 2>&1
    "$AUDIT_PY" integer_refine.py --init "parameters/${tag}_integer.npz" \
        --out "parameters/${tag}_unbounded_integer.npz" --steps 40 --max-rows 64 --cycles 2 --unbounded \
        > "results/${tag}_unbounded_integer.log" 2>&1
done
AUDIT_DEVICE=cpu "$AUDIT_PY" integer_grow.py --start 511 --stop 639 --batch 32 --lambda-cap 320 \
    --tag integer_growth > results/integer_growth.log 2>&1
