#!/bin/bash
PYT=/mnt/data01/miniconda3/envs/torch2.13/bin/python
run(){ $PYT init_eig.py --R $1 --B $2 --C $3 --interleave --out e_r$1b$2c$3.npz; $PYT shapefit.py --R $1 --B $2 --C $3 --steps 100000 --init e_r$1b$2c$3.npz --out c_r$1b$2c$3.npz > c_r$1b$2c$3.log 2>&1; }
run 2 360 175 &
run 4 180 128 &
wait
run 4 256 88 &
run 1 511 350 &
wait
