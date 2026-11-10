#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Wait for the already-live expected-noise comparison, without restarting it.
while kill -0 2246 2>/dev/null; do sleep 10; done
source env.sh
export PYTHONPATH="$PROJECT_ROOT/initial_experiments/5/chem-qfit:$PYTHONPATH"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 AUDIT_THREADS=2
/mnt/data01/miniconda3/envs/torch2.13/bin/python quantized_fixed_fit.py --ste --init parameters/robust_B511_last.npz --seed 142039 --steps 16000 --lr .0003 --final-lr .000001 --lambda-reg 300 --epsilon-reg .0025 --weight-cap 5 --lambda-cap 330 --hard-lambda-penalty 10 --precision float32 --checkpoint 1000 --tag ste_B511 > results/ste_B511.log 2>&1
