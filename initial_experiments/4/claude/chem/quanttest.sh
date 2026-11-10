#!/bin/bash
# 8-bit quantization test: lattice anneal then curvature/pair integer descent.
set -e
in=$1; tag=$2
export PYTHONPATH=/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers:/media/shared/ExpansionDrive2/QuantumAI/from_host/kickmix
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib
PYT=/mnt/data01/miniconda3/envs/torch2.13/bin/python
$PYT blockanneal.py --init $in --out an_$tag.npz --columns 175 --steps 4000 --lr 0.02 > an_$tag.log 2>&1
$PYT block_coordinate.py --init an_$tag.npz --out qd_$tag.npz --columns 175 --steps 2000 --moves 256 > qd_$tag.log 2>&1
