#!/bin/bash
PYT=/mnt/data01/miniconda3/envs/torch2.13/bin/python
for s in "1 719 350" "2 360 175" "4 180 88" "4 180 128" "8 90 44"; do
  set -- $s
  $PYT shapefit.py --R $1 --B $2 --C $3 --steps 200000 --out s_r$1b$2c$3.npz > s_r$1b$2c$3.log 2>&1
done
