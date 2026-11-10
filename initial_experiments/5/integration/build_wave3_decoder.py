"""Integrate the verified trimmed middle-product emitter with the gap decoder.

Keep the frozen gap schedule, low-workspace inversions, scan, and dispatch.
Only the bilinear middle-product algorithm changes. Native acceptance follows
in validate_wave3_decoder.py; source experiments remain read-only.
"""

import hashlib
import importlib.util
import json
import os
from pathlib import Path

from probe_gap_trimming import GAP, HERE, ROUND5, decoder, frozen_frames


def sha(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def main():
    destination = HERE / "wave3"
    for name in ("candidates", "results", "parameters"):
        (destination / name).mkdir(parents=True, exist_ok=True)
    original = json.loads((GAP / "parameters/best_candidate_manifest.json").read_text())
    for name, digest in original["source_sha256"].items():
        if sha(GAP / name) != digest:
            raise ValueError(f"Gap source changed since acceptance: {name}")
    trimmed_path = ROUND5 / "dqi-core/components/trimmed_middle.py"
    spec = importlib.util.spec_from_file_location("wave3_trimmed_middle", trimmed_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    decoder.BoundedArithmetic = module.TrimmedArithmetic
    decoder.frames = frozen_frames
    decoder.build.slot_counts = json.loads((GAP / "parameters/best_schedule.json").read_text())["counts"]
    profile = json.loads((destination / "results/gap_trimmed_profile.json").read_text())
    if profile["Q"] > original["Q"] or profile["saved_core_T"] <= 0:
        raise ValueError("Resource probe does not justify building a replacement")
    sources = {str((GAP / path).relative_to(ROUND5)): digest
               for path, digest in original["source_sha256"].items()}
    for path in (Path(__file__), HERE / "probe_gap_trimming.py", trimmed_path,
                 GAP / "parameters/best_schedule.json"):
        sources[str(path.relative_to(ROUND5))] = sha(path)
    manifest = {
        "status": "FROZEN_BEFORE_NATIVE_VERIFICATION",
        "scope": "fixed accepted gap schedule, trimmed middle products only; original inverse strategy and phase kernel",
        "sources_sha256_relative_to_round5": sources,
        "original_circuit_sha256": original["sha256"],
        "expected_Q": profile["Q"],
        "expected_core_T": profile["T"],
        "expected_static_T": original["static_Toffolis"] - profile["saved_core_T"],
        "host_verifier_sha256": sha(ROUND5.parents[1] / "from_host/verifiers/dqi/verify_submission_dqi.py"),
    }
    os.chdir(destination)
    Path("parameters/gap_trimmed_pre_native.json").write_text(json.dumps(manifest, indent=2) + "\n")
    path = Path("candidates/gap_trimmed_n70_m4095.kmb")
    circuit = decoder.build(70, 3, 100, 4095, path, linear_padding=True,
                            virtual_zero=True, virtual_lead=True, low_inverse=True,
                            lazy_flags=True)
    stats = circuit.operation_statistics()
    if circuit.num_qubits != manifest["expected_Q"] or stats.num_toffoli_like != manifest["expected_static_T"]:
        raise ValueError("Compiled circuit differs from resource-only profile")
    manifest.update(status="BUILT_REQUIRES_NATIVE_VERIFICATION", circuit=str(path),
                    circuit_sha256=sha(path), Q=circuit.num_qubits,
                    static_Toffolis=stats.num_toffoli_like)
    Path("parameters/gap_trimmed_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
