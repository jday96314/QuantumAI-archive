"""Coherent choice of r/(1+x^2*r) or r/(1+x^3*r).

Insert one clean coefficient by a controlled rotation. The ordinary branch
omits writes to its extra last coefficient; triangularity ensures that this
does not affect lower coefficients. Only that coefficient's writes need a
quantum control, rather than controlling the whole reciprocal algorithm.
"""
from contextlib import contextmanager


@contextmanager
def controlled_targets(m,E,targets,control):
    selected=set(targets);original=E.emit_mac
    def emit(terms):
        written={q for term in terms for q in term[-1]}
        special=sorted(written&selected)
        if not special:original(terms);return
        read={q for term in terms for q in (term[1]+(term[3] if term[0]=='P' else ()))}
        if not (read&written):
            # The polynomial middle product reads a lower band and writes
            # a disjoint upper band. Aggregate its selected output before
            # applying the control, sharing all original product fanouts.
            with m.opened(E,'conditional_tail_output',len(special)) as bank:
                mapping=dict(zip(special,bank))
                routed=[(*term[:-1],tuple(mapping.get(q,q) for q in term[-1])) for term in terms]
                repayment=[(*term[:-1],tuple(mapping[q] for q in term[-1] if q in mapping)) for term in terms]
                repayment=[term for term in repayment if term[-1]]
                original(routed)
                for q,t in mapping.items():E.b.ccx(control,t,q)
                E.emit_free(repayment)
            return
        for term in terms:
            special=tuple(q for q in term[-1] if q in selected)
            ordinary=tuple(q for q in term[-1] if q not in selected)
            if not special:original([term]);continue
            if term[0]=='L':
                if ordinary:original([(*term[:-1],ordinary)])
                conditional=m.norm_term((control,),0,term[1],term[2],special)
                if conditional:original([conditional])
            else:
                def fanout(fan):
                    for q in ordinary:E.b.cx(fan,q)
                    with m.opened(E,'conditional_tail_fan',1) as aux:
                        E.b.ccx(control,fan,aux[0])
                        for q in special:E.b.cx(aux[0],q)
                        measurement=E.b.hmr(aux[0]);E.b.cz(control,fan,cond=measurement)
                E.fan_pair(term[1],term[2],term[3],term[4],fanout)
    E.emit_mac=emit
    try:yield
    finally:E.emit_mac=original


def apply(m,E,math,R,shift_extra):
    if not R:return
    with m.opened(E,'variable_shift_extra_coefficient',12) as zero:
        array=[*R,zero]
        swaps=[pair for i in reversed(range(len(R))) for pair in zip(array[i],array[i+1])]
        E.swaps(shift_extra,swaps)
        present=() in math.known_zero;math.known_zero.add(())
        with controlled_targets(m,E,zero,shift_extra):math.inverse([(),(),*array])
        if not present:math.known_zero.remove(())
        E.swaps(shift_extra,list(reversed(swaps)))
