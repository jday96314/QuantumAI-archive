# QuantumAI

## Circuit verification and scoring

Run from this directory. The bundled Linux x86-64 `kickmix.so` requires
**Python 3.13** and a C++ runtime providing **GLIBCXX_3.4.32** (GCC 13 or newer).

```bash
python3.13 -m venv /tmp/quantumai-venv
source /tmp/quantumai-venv/bin/activate
python -m pip install numpy==2.5.4 scipy==1.18.1
export PYTHONPATH="$PWD/from_host/kickmix:$PWD/from_host/verifiers${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
python -c 'import kickmix, numpy, scipy'

python from_host/verifiers/invsqrt/verify_submission_invsqrt.py from_host/sample_submission/invsqrt.kmb --max-workers 8
python from_host/verifiers/dqi/verify_submission_dqi.py from_host/sample_submission/dqi.kmb --max-workers 4
python from_host/verifiers/dfthcsos/verify.py from_host/sample_submission/dfthcsos
```

If importing `kickmix` reports a missing `GLIBCXX` version, set
`LD_LIBRARY_PATH` to the directory containing a compatible `libstdc++.so.6`.
The venv lives outside the shared drive because that filesystem rejects the
symlinks needed by `venv`.

Replace the sample paths with your circuit files (`.kmb` or `.kmx`) or chemistry
submission directory (`angle_k_0`, … plus `submission.npz`). Each command checks
correctness and prints JSON containing `status`, `score`, and resource counts;
a failed check exits nonzero. **Lower raw task scores are better.** The competition
combines improvements relative to the baselines using a geometric mean.

Keep the default problem sizes and shot counts for comparable scores. Inverse
square root checks all 3,133,452 inputs 64 times; DQI uses 102,400 shots; chemistry
uses 64 shots per address. Reduced checks change confidence or synthesis margins
and can change the score. `--salt` controls reproducible hashed seeds; `--seed`
overrides them. Run each command with `--help` for options.

See [baseline scores and profiling](docs/baseline-profiling.md) for measured
costs, runtime bottlenecks, and profiling commands, and the host's
[chemistry scoring guide](from_host/verifiers/dfthcsos/SCORING_README.md) for
its submission contract and cost model.
