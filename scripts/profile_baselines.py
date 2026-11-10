"""Run the unchanged host verifiers, optionally profiling parent and fork workers.

Run from the repository root with from_host/kickmix and from_host/verifiers
on PYTHONPATH. Outputs include the complete verifier result and wall time.
"""

import argparse
import cProfile
import json
import os
from pathlib import Path
import platform
import time
from functools import wraps

import kickmix as km
import numpy as np
import scipy


TIMINGS = {}


def summarize_workers(directory):
    summary = {"wall_seconds": 0.0, "tasks": 0, "timings": {}}
    for path in directory.glob("worker-*.json"):
        worker = json.loads(path.read_text())
        summary["wall_seconds"] += worker["wall_seconds"]
        summary["tasks"] += 1
        for name, value in worker["timings"].items():
            record = summary["timings"].setdefault(name, {"calls": 0, "seconds": 0.0})
            record["calls"] += value["calls"]
            record["seconds"] += value["seconds"]
    return summary


def timed(name, fn):
    @wraps(fn)
    def measured(*args, **kwargs):
        start = time.perf_counter()
        try:
            return fn(*args, **kwargs)
        finally:
            record = TIMINGS.setdefault(name, {"calls": 0, "seconds": 0.0})
            record["calls"] += 1
            record["seconds"] += time.perf_counter() - start
    return measured


def instrument(verifier, task):
    """Time native boundaries explicitly: cProfile omits some nanobind calls."""
    original_simulator = km.Simulator

    class TimedSimulator:
        def __init__(self, *args, **kwargs):
            self._sim = timed("Simulator.__init__", original_simulator)(*args, **kwargs)

        def __getattr__(self, name):
            member = getattr(self._sim, name)
            if callable(member):
                member = timed(f"Simulator.{name}", member)
                setattr(self, name, member)
            return member

    km.Simulator = TimedSimulator
    names = (
        "_circuit_sha256", "validate_registers", "reference_values", "sample_shots",
        "to_bits", "from_bits", "load_inputs", "load_tree_edges", "load_circuits",
        "submission_sha256", "recover_words_from_circuits", "evaluate_words",
        "resource_costs", "load_reference", "check_fcidump", "load_factorization",
        "load_fixed_shifts",
    )
    for name in names:
        if hasattr(verifier, name):
            setattr(verifier, name, timed(name, getattr(verifier, name)))
    if task == "dfthcsos":
        from dfthcsos.scoreinfra import circuits
        circuits._reject_rotations = timed("_reject_rotations", circuits._reject_rotations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("invsqrt", "dqi", "dfthcsos"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--timings", action="store_true")
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--invsqrt-max-x", type=int)
    parser.add_argument("--dqi-num-shots", type=int)
    parser.add_argument("--summarize-only", action="store_true",
                        help="Aggregate existing worker JSON into result.json without rerunning")
    args = parser.parse_args()
    if args.summarize_only:
        path = args.output / "result.json"
        document = json.loads(path.read_text())
        document["worker_summary"] = summarize_workers(args.output)
        path.write_text(json.dumps(document, indent=2) + "\n")
        return
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("--output must be a new or empty directory to avoid mixing runs")
    args.output.mkdir(parents=True, exist_ok=True)
    baseline = Path(__file__).resolve().parents[1] / "from_host/sample_submission"
    if args.task == "invsqrt":
        from invsqrt import verify_submission_invsqrt as verifier
        spec = verifier.ProblemSpecInvSqrt.for_n_p(10)
        if args.invsqrt_max_x is not None:
            from dataclasses import replace
            spec = replace(spec, max_x=args.invsqrt_max_x)
    elif args.task == "dqi":
        from dqi import verify_submission_dqi as verifier
        spec = verifier.ProblemSpecDQI.benchmark_instance()
        if args.dqi_num_shots is not None:
            from dataclasses import replace
            spec = replace(spec, num_shots=args.dqi_num_shots)
    else:
        from dfthcsos import verify as verifier

    # The host forks closures. Wrap each task before passing it to the host's
    # existing dispatcher so child-process work appears in separate profiles.
    if args.timings:
        instrument(verifier, args.task)
    if (args.profile or args.timings) and args.task != "dfthcsos":
        original_parallel = verifier._run_parallel

        def profiled_parallel(fn, tasks, max_workers):
            def profiled_task(task):
                # Python 3.13 permits only one active cProfile instance per
                # thread; forked workers also inherit the parent's profiler.
                if profile:
                    profile.disable()
                task_profile = cProfile.Profile() if args.profile else None
                saved_timings = dict(TIMINGS)
                TIMINGS.clear()
                start = time.perf_counter()
                try:
                    return task_profile.runcall(fn, task) if task_profile else fn(task)
                finally:
                    wall = time.perf_counter() - start
                    stem = args.output / f"worker-{os.getpid()}-{task[0]}"
                    if task_profile:
                        task_profile.dump_stats(str(stem) + ".prof")
                    if args.timings:
                        Path(str(stem) + ".json").write_text(json.dumps({"wall_seconds": wall, "timings": TIMINGS}, indent=2) + "\n")
                    TIMINGS.clear()
                    TIMINGS.update(saved_timings)
                    if profile:
                        profile.enable()
            return original_parallel(profiled_task, tasks, max_workers)

        verifier._run_parallel = profiled_parallel

    profile = cProfile.Profile() if args.profile else None
    start = time.perf_counter()
    if profile:
        profile.enable()
    if args.task == "dfthcsos":
        result = verifier.verify(baseline / args.task)
    else:
        circuit = timed("circuit_load", km.Circuit.from_file)(baseline / f"{args.task}.kmb", format="kmb")
        result = verifier.verify(circuit, spec, max_workers=args.max_workers)
    if profile:
        profile.disable()
    wall_seconds = time.perf_counter() - start
    if profile:
        profile.dump_stats(str(args.output / "parent.prof"))
    document = {
        "task": args.task,
        "wall_seconds": wall_seconds,
        "profiled": args.profile,
        "timings": TIMINGS,
        "max_workers": args.max_workers if args.task != "dfthcsos" else 1,
        "invsqrt_max_x": args.invsqrt_max_x,
        "dqi_num_shots": args.dqi_num_shots,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "result": result,
    }
    if args.timings and args.task != "dfthcsos":
        document["worker_summary"] = summarize_workers(args.output)
    (args.output / "result.json").write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps(document, indent=2), flush=True)


if __name__ == "__main__":
    main()
