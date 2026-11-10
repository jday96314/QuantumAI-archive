"""Keep inversion intermediates only when the selected peak has room."""
from contextlib import contextmanager
from low_workspace_inverse import bind as low_bind

def bind(m,limit=961):
    low=low_bind(m)
    @contextmanager
    def inverse(E,source):
        # Callers in the selected decoder need at most twelve additional
        # body scratch bits while holding this inverse (checked by replay).
        if E.b.active_qubits+42<=limit:
            with m.opened(E,'retained_inverse_output',12) as out:
                with m.chain_live(E,source,out):yield out
        else:
            with low(E,source) as out:yield out
    return inverse
