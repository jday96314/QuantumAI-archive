"""Retain inversion scratch only across source-preserving decoder contexts.

Gap opening and selected cofactor paths exchange source and reciprocal. Their
low-workspace inverse supports that exchange; chain_live retains intermediates
for the original source and requires source preservation. Keep the original
primitive at those sites. The caller inspection is compile-time only, over
hash-bound Python emitters, and never reads circuit inputs or verifier seeds.
"""
from contextlib import contextmanager
import inspect
from pathlib import Path

from low_workspace_inverse import bind as low_bind


def bind(module,limit=1004):
    low=low_bind(module)
    @contextmanager
    def retained(E,source,reserve):
        if E.b.active_qubits+30+reserve<=limit:
            with module.opened(E,'safe_retained_inverse_output',12) as out:
                with module.chain_live(E,source,out):yield out
        else:
            with low(E,source) as out:yield out

    def inverse(E,source):
        caller=inspect.currentframe().f_back
        filename=Path(caller.f_code.co_filename).name
        reserve=12
        if filename=='quantum_control.py':
            reserve=36
            caller=caller.f_back
            filename=Path(caller.f_code.co_filename).name
        if filename in ('serial_opening.py','gap_serial.py'):
            return low(E,source)
        if filename=='onepass_gap.py' and caller.f_code.co_name=='boundary':
            reserve=36
        return retained(E,source,reserve)
    return inverse
