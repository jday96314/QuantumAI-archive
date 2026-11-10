#!/usr/bin/env bash
# Fresh searches may change float64 low-order digits; delivered checkpoints
# and native digests, rather than a rerun's timing, bind the reported results.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source env.sh
"$TREE_PY" test_trees.py > results/encoding_tests.log 2>&1
"$TREE_PY" test_equivalence.py > results/binary_equivalence_tests.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" screen.py --cycles 12 > results/screen_reproduced.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" topology_search.py --control b351_float --tag topology_b351 --steps 240 --cycles 12 > results/topology_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" topology_search.py --control b351_float --tag random_b351 --mode random --steps 240 --cycles 12 > results/random_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" phase_grid.py --control b351_float --tag grid9_b351 --passes 2 --cycles 12 > results/grid9_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" topology_search.py --mode joint --control b479_exact --tag joint_b479 --steps 320 --cycles 8 --temperature 0.00015 > results/joint_b479.log 2>&1
TREE_DEVICE=cpu "$TREE_PY" local_nni.py --control b479_exact --tag nni_b479 --rounds 2 --support 8 --cycles 12 > results/nni_b479.log 2>&1
TREE_DEVICE=cpu "$TREE_PY" local_nni.py --control b639_exact --tag nni_b639 --rounds 2 --support 8 --cycles 12 > results/nni_b639.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" descend.py --init candidates/topology_b351.npz --tag repair_topology_b351 --steps 120 --cycles 6 --refit-budget 4000 > results/repair_topology_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" descend.py --init candidates/b351_float__incumbent.npz --tag repair_incumbent_b351 --steps 120 --cycles 6 --refit-budget 4000 > results/repair_incumbent_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" refit_control.py --init candidates/topology_b351.npz --tag fixed_topology_b351 --budget 4600 > results/fixed_topology_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" refit_control.py --init candidates/b351_float__incumbent.npz --tag fixed_incumbent_b351 --budget 4600 > results/fixed_incumbent_b351.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" descend.py --init candidates/nni_b479_top0.npz --tag repair_nni_b479 --steps 60 --cycles 6 --refit-budget 1600 > results/repair_nni_b479.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" descend.py --init candidates/b479_exact__incumbent.npz --tag repair_original_b479 --steps 60 --cycles 6 --refit-budget 1600 > results/repair_original_b479.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" psd_capacity.py --init candidates/repair_nni_b479.npz --tag psd_repair_nni_b479 --steps 2200 > results/psd_repair_nni_b479.log 2>&1
TREE_DEVICE=cuda "$TREE_PY" psd_capacity.py --init candidates/repair_original_b479.npz --tag psd_repair_original_b479 --steps 2200 > results/psd_repair_original_b479.log 2>&1
"$TREE_PY" physical.py candidates/topology_b351_polished.npz candidates/grid9_b351.npz candidates/b639_float__huffman_tensor.npz candidates/joint_b479_polished.npz candidates/repair_topology_b351.npz candidates/repair_incumbent_b351.npz candidates/repair_nni_b479.npz candidates/repair_original_b479.npz --jacobi 8 > results/physical_reproduced.log 2>&1
"$TREE_PY" compare_scores.py > results/surrogate_comparison.log 2>&1
"$TREE_PY" span_robustness.py > results/span_robustness.log 2>&1
"$TREE_PY" plot_results.py > results/plot_results.log 2>&1
