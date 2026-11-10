"""Generic unit-series Schur transform: F -> normalized Euclid quotients.

Each pair holds beta, alpha. Nonzero pivots are required. Research only.
"""
def apply(m,E,math,F,inverse=False,flags=None,record=True):
    if not F:return
    assert len(F)%2==0
    if inverse:apply(m,E,math,F[2:],True,None if flags is None else flags[1:],record)
    if not inverse and flags is not None:
        if record:
            with m.eq_test(E,F[0],0) as top:E.cx(top,flags[0])
        E.cx(flags[0],F[0][0])
    with m.inv_word(E,F[0]) as iv,m.opened(E,'cf_product',12) as product:
        if not inverse:
            for word in F[1:]:m.mon_slot(E,word,iv,F[0],product)
        # The reciprocal emitter treats F[0] as implicit one and never reads it.
        math.inverse(F)
        if inverse:
            for word in F[1:]:m.mon_slot(E,word,F[0],iv,product)
    if inverse and flags is not None:
        E.cx(flags[0],F[0][0])
        if record:
            with m.eq_test(E,F[0],0) as top:E.cx(top,flags[0])
    if not inverse:apply(m,E,math,F[2:],False,None if flags is None else flags[1:],record)
