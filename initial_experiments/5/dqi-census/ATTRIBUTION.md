# Sources and attribution

All original host files and earlier experiment files are read-only. The field
emitter, polynomial evaluation/interpolation products, allocator, inverse-chain
terms and MBU conventions are reused from andigu/uncompute, Apache-2.0,
revision64b5a13aa249f1360846657392e6c66817b26d0f. LICENSE-uncompute is copied
here; the complete original source remains at ../../3/vendor/uncompute.

native_components.py adapts the modular XOR-MAC construction in
../dqi-phase/fast_modular.py and its explicit inversion chain in
../dqi-phase/modular_schedule.py. The new square specialization performs
Frobenius as a linear operand map, then pays the full quantum-modulus reduction.
The Frobenius inversion route was already known to be wide; this census's
new optimized counts do not claim discovery of that architecture.

The rational-coordinate reference and quadratic direct-F' norm circuits are
implemented in this directory. Unit triangular arithmetic and low-memory
scalar inversion import ../../4/codex/{mobius_series,series_inverse,
low_workspace_inverse}.py. The native compiler is the unchanged round-4
translate_prefix. Host-basis references, polynomial arithmetic and trace duals
import ../dqi-phase/oracle.py; no earlier writer/main function is invoked.
Field-basis conversion imports ../../3/published_dqi.py.

Primary reading:

- Khattar et al., Verifiable Quantum Advantage via Optimized DQI Circuits,
  https://arxiv.org/html/2510.10967v1. SectionsII.2,II.4, AppendixC/D inform
  the state-sharing/direct-division/MBU discussion. The stated asymptotic
  bounds do not certify this host's finite-size RQ²T objective.
- https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md.
  The fixed-mask sweep and published TQ count require adaptation to this host.
  A downloaded snapshot and hash are retained. It is not a new host-native
  verification result.

Prior result attribution: read-only round-5 repair, phase, hybrid, core and
concurrently changing dispatch RESULTS; round-4 REPORT and Codex continuation.
Liveness/splice counts are inherited exact emitter evidence, not new BM runs.
No dqi-dispatch implementation file is imported or built. Source/result hashes
are recorded in results/source_manifest.json.
