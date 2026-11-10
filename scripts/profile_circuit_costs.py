"""Attribute baseline gate cost to the host's INFO_TAG subroutine annotations.

Static counts include every branch. Executed counts use a reproducible batch
of 256 inputs and preserve simulator state between complete top-level stages.
Checks the segmented simulation against the original circuit before reporting.
"""

import argparse
from collections import Counter
import gc
import io
import json
from pathlib import Path
import tempfile
import time

import kickmix as km
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("invsqrt", "dqi"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    circuit = km.Circuit.from_file(root / f"from_host/sample_submission/{args.task}.kmb", format="kmb")
    num_qubits = circuit.num_qubits
    rng = np.random.default_rng(20261022)
    if args.task == "invsqrt":
        from invsqrt import verify_submission_invsqrt as v
        spec = v.ProblemSpecInvSqrt.for_n_p(10)
        x, _ = v.validate_registers(circuit, spec)
        inputs = [(x, v.to_bits(rng.integers(1, spec.max_x + 1, size=256, dtype=np.uint64), spec.width))]
    else:
        from dqi import verify_submission_dqi as v
        spec = v.ProblemSpecDQI.benchmark_instance()
        s, e, _ = v.validate_registers(circuit, spec)
        syndromes, mbuc, _ = v.sample_shots(spec, v.FieldTables(km.GF2Field(spec.degree)), rng, 256)
        inputs = [(km.array(s), v.to_bits(syndromes, spec.degree)), (km.array(e), v.to_bits(mbuc, spec.degree))]

    def simulator():
        sim = km.Simulator(batch_size=256, counting_operations=True, seed=20261022)
        sim.use_same_registers_as(circuit)
        for register, bits in inputs:
            sim.write_across_shots(register, bits)
        return sim

    reference, segmented = simulator(), simulator()
    reference.do(circuit)
    expected_stats = dict(reference.operation_statistics())
    expected_bits = reference.read_across_shots(km.array([km.q(i) for i in range(num_qubits)]))
    expected_phases = [reference.read_shot_phase(i) for i in range(256)]
    static = Counter()
    stages = {}
    stack = []
    lines = io.StringIO()
    active = "outside_marks"

    def flush():
        nonlocal lines
        text = lines.getvalue()
        if not text:
            return
        fragment = km.Circuit(text)
        start = time.perf_counter()
        segmented.do(fragment)
        seconds = time.perf_counter() - start
        counts = dict(segmented.operation_statistics())
        # operation_statistics is cumulative, so record the delta per stage.
        previous = flush.previous
        delta = {op: (n - previous.get(op, 0)) / 256 for op, n in counts.items()}
        record = stages.setdefault(active, {"executed_operations_per_shot": {}, "simulation_seconds": 0.0})
        for op, n in delta.items():
            record["executed_operations_per_shot"][op] = record["executed_operations_per_shot"].get(op, 0) + n
        record["simulation_seconds"] += seconds
        flush.previous = counts
        lines = io.StringIO()

    flush.previous = {}
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8") as source:
        circuit.to_file(source, format="kmx")
        del circuit
        gc.collect()
        source.seek(0)
        for line in source:
            if line.startswith('INFO_TAG "@enter:'):
                stack.append(line.split(":", 1)[1].rstrip().removesuffix('"'))
            elif line.startswith('INFO_TAG "@exit:'):
                name = line.split(":", 1)[1].rstrip().removesuffix('"')
                assert stack.pop() == name
            label = "/".join(stack[:2] if args.task == "invsqrt" else stack[:1]) or "outside_marks"
            if label != active:
                flush()
                active = label
            lines.write(line)
            if line.startswith(("CCX ", "CCZ ", "CSWAP ")):
                static[label] += 1
                # Inclusive nested paths identify the arithmetic hot spots.
                for depth in range(1, len(stack) + 1):
                    static["inclusive:" + "/".join(stack[:depth])] += 1
        flush()
    assert not stack
    assert dict(segmented.operation_statistics()) == expected_stats, "Segmented operation counts differ"
    assert np.array_equal(segmented.read_across_shots(km.array([km.q(i) for i in range(num_qubits)])), expected_bits)
    assert [segmented.read_shot_phase(i) for i in range(256)] == expected_phases
    for record in stages.values():
        ops = record["executed_operations_per_shot"]
        record["toffoli_per_shot"] = sum(ops.get(op, 0) for op in ("CCX", "CCZ", "CSWAP"))
    document = {"task": args.task, "shots": 256, "input_and_simulator_seed": 20261022,
                "matches_original_simulation": True, "num_qubits": num_qubits,
                "static_toffoli_by_path": dict(static), "stages": stages}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps({k: v for k, v in document.items() if k != "static_toffoli_by_path"}, indent=2))


if __name__ == "__main__":
    main()
