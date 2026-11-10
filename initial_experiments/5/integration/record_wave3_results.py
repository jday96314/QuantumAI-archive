"""Write the report-review decision and final native-bound integration result."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROUND5 = HERE.parent
OUT = HERE / "wave3"


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def main():
    summary = read(OUT / "submission.summary.json")
    old = read(HERE / "submission.summary.json")
    acceptance = read(OUT / "results/gap_trimmed_acceptance.json")
    metric = read(OUT / "results/metric_preflight_and_score.json")
    audit = read(OUT / "results/archive_audit.json")
    comparison = read(OUT / "results/new_default_stream_comparison.json")
    assert summary["status"] == metric["status"] == audit["status"] == "PASS"
    assert acceptance["status"] == "VERIFIED_IMPROVEMENT"
    assert comparison["status"] == "PASS"
    assert summary["sha256"] == metric["archive_sha256"] == audit["archive_sha256"]
    assert summary["overall_score"] > old["overall_score"]
    assert sha(HERE / "submission.zip") == old["sha256"], "Previously submitted archive changed"
    review = {
        "claude/chem_audit_pruning_and_submission_report.md": "Zero subtree-pruning savings; its banked physical improvement is already present. Its packaged trimmed DQI is slower than the accepted gap decoder.",
        "claude/chem_census_factorization_and_compilation_report.md": "Selected: verified structured compiler applied to the banked physical fit, 696942846.9115301 chemistry cost.",
        "claude/dqi_core_trimmed_decoder_report.md": "Ported only the accepted trimmed middle-product emitter to the gap decoder; directly replacing the gap decoder with the old trimmed architecture would regress cost.",
        "claude/dqi_dispatch_exception_correction_report.md": "No complete accepted replacement; resource projections omit required work.",
        "claude/dqi_census_decoder_architecture_report.md": "Working research components and conditional budgets; no full verified replacement.",
        "claude/inv_census_inverse_architecture_report.md": "No lower-cost verified inverse; retain the current circuit.",
        "gpt/chem_qrom_content_aware_compilation_report.md": "Verified 697327965.8381622 compiler gain; already composed with better banked normalization in chem-census, so do not count it twice.",
        "gpt/chem_tree_topology_and_quantization_report.md": "No physically valid cheaper replacement.",
        "gpt/dqi_core_phase_aware_decoder_report.md": "Incomplete gap handling and conversion in candidate architectures; no accepted full-task gain.",
        "gpt/inv_phase_inverse_square_root_report.md": "Both full native candidates pass but cost more than the current inverse.",
    }
    report_root = ROUND5 / "brainstorming/results/wave3"
    decisions = [{"report": str((report_root / path).relative_to(ROUND5)),
                  "sha256": sha(report_root / path), "decision": decision}
                 for path, decision in review.items()]
    (OUT / "results/report_selection_review.json").write_text(json.dumps(decisions, indent=2) + "\n")
    native = {task: read(OUT / f"results/{task}_verify.json") for task in summary["task_scores"]}
    saved = acceptance["saved_static_Toffolis"]
    gain = (summary["overall_improvement"] - 1) * 100
    reductions = {task: (1 - 1 / ratio) * 100 for task, ratio in summary["task_improvements"].items()}
    lines = [
        "# Round 5 wave-three integrated submission", "",
        f"**Verified local score: {summary['overall_score']:.12f}**, up from **{old['overall_score']:.12f}** ({gain:.6f}% higher).",
        "The deliverable is [submission.zip](submission.zip). The previously submitted archive and its report remain unchanged one directory above.", "",
        "## Integrated improvements", "",
        "| Task | Previous task cost | New native task cost | Cost reduction |",
        "| --- | ---: | ---: | ---: |",
    ]
    for task in summary["task_scores"]:
        lines.append(f"| {task} | {old['task_scores'][task]:,.10f} | {summary['task_scores'][task]:,.10f} | {reductions[task]:.6f}% |")
    lines += [
        "", f"The accepted wave-three middle-product emitter removes **{saved:,} static Toffolis** from the existing gap-aware DQI decoder, at the same **1,004-qubit peak**. The gap-aware schedule, dispatch, original low-workspace inverse, and phase scan are retained. The complete compiled count is **1,172,538 static Toffolis** (801,388 core + 371,150 scan).",
        "", "Chemistry uses `chem-census/candidates/structured`: the content-aware QROM compiler composed with the banked physical fit. It uses 19 ANGLE circuits, 31,608 host Toffolis and 454 host qubits. This already incorporates the overlapping QROM and normalization gains.",
        "", "The inverse remains the round-four `zero_head_affine.kmb`; neither inverse report supplies a cheaper accepted replacement.",
        "", "## Review of all ten reports", "",
        "| Report | Integration decision |", "| --- | --- |",
    ]
    for path, decision in review.items():
        lines.append(f"| [{Path(path).stem}](../../brainstorming/results/wave3/{path}) | {decision} |")
    lines += [
        "", "## Verification and score", "",
        "The DQI candidate was fixed from resource counts before any native sampling. Its source manifest records the accepted gap sources and the exact emitter from the accepted trimmed-decoder build. No budgeted inverse-retention change was ported.",
        "", "| DQI run | Shots | Value / phase errors | Mean executed Toffolis | Native cost |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for name in ("default", "independent", "paired"):
        report = read(OUT / f"results/gap_trimmed_verify_{name}.json")
        lines.append(f"| {name} | {report['verified_shots']:,} | {report['value_errors']} / {report['phase_errors']} | {report['toffoli_count']:,.9f} | {report['score']:,.6f} |")
    chemistry = native["dfthcsos"]
    lines += [
        "", f"The paired run uses the original gap decoder's default hash-derived seed. It reproduces the same decoding counters, repetition factor, and width, with exactly **{saved:,} fewer mean Toffolis**. Default and independent-salt runs pass separately; the default supplies the score. Confidence remains delta = 2^-128 and each full run uses 102,400 shots.",
        "", f"The new default stream has **{comparison['phase_errors']} phase error and {comparison['value_errors']} value errors**. A further full 102,400-shot replay of the original gap decoder on that same stream reproduces those counters and the repetition factor, with exactly {saved:,} more mean Toffolis. The reported score uses the default stream with its actual confidence penalty. See the [controlled stream comparison](results/new_default_stream_comparison.json).",
        "", "The full-mask adversarial test passes 1,024 cases covering all gap stages, immediate completion, zero/one-error, weight-35, and random controls, including twice-executed restoration. All five known residual examples pass 640 independent-mask input/scratch and twice-execution restoration checks, retaining the original phase-error counts. The decoder still approximates the rare degree-three and adjacent-double-gap cases.",
        "", f"Fresh chemistry verification passes all **{chemistry['verified_shots']:,} shots** at default settings, matching the two independent native passes in its source experiment. The inverse's unchanged bytes reuse their **200,540,928-shot PASS** evidence. No reduced workloads or favorable-salt selection are used.",
        "", "The package audit reopens the ZIP, reads every member to EOF, verifies hashes against the selected source files and extracted artifacts, and matches the native DQI/inverse circuit hashes and canonical chemistry digest. The unmodified bundled metric's circuit and NPZ readers pass, and its `score` function reproduces the combined score from the hash-bound native costs. This is a local default-stream score using the same baseline costs as the prior integration.",
        "", f"Archive size: **{summary['archive_bytes']:,} bytes**. Archive SHA256:", "", "```text", summary["sha256"], "```",
        "", "Machine-readable records:", "",
        "- [Submission summary](submission.summary.json)",
        "- [Native DQI acceptance](results/gap_trimmed_acceptance.json)",
        "- [Adversarial restoration](results/gap_trimmed_adversarial.json) and [known residuals](results/gap_trimmed_residual.json)",
        "- [Archive audit](results/archive_audit.json) and [unmodified metric readers/score](results/metric_preflight_and_score.json)",
        "- [All report decisions and source hashes](results/report_selection_review.json)",
        "", "## Reproduction", "", "From the repository root:", "", "```bash",
        "source initial_experiments/5/integration/env.sh",
        '"$EXP_PYTHON" initial_experiments/5/integration/probe_gap_trimming.py',
        '"$EXP_PYTHON" initial_experiments/5/integration/build_wave3_decoder.py',
        '"$EXP_PYTHON" initial_experiments/5/integration/validate_wave3_decoder.py',
        '"$EXP_PYTHON" initial_experiments/5/integration/check_wave3_default_stream.py',
        'flock /tmp/r6-fullverify.lock "$EXP_PYTHON" from_host/verifiers/dfthcsos/verify.py initial_experiments/5/chem-census/candidates/structured > initial_experiments/5/integration/wave3/results/chemistry_verify_default.json',
        "bash initial_experiments/5/integration/package_wave3.sh", "```", "",
        "`package_wave3.sh` packages the accepted native-verified inputs and audits the exact archive. It reuses the already completed full native checks after binding their digests to the packaged bytes. Re-running that script alone does not repeat simulation. Reused field arithmetic retains its upstream Apache-2.0 attribution in the original experiment reports.", "",
    ]
    (OUT / "RESULTS.md").write_text("\n".join(lines))
    print(f"Final wave-three report written; score {summary['overall_score']:.12f}, gain {gain:.6f}%.")


if __name__ == "__main__":
    main()
