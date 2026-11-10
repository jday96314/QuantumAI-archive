"""Retain one zero-predicate bit instead of its entire conjunction ladder."""
from contextlib import contextmanager


@contextmanager
def armed(m,E,word):
    with m.opened(E,'compact_zero_predicate',1) as flag:
        with m.eq_test(E,word,0) as test:E.cx(test,flag[0])
        E.cx(flag[0],word[0])
        yield
        E.cx(flag[0],word[0])
        with m.eq_test(E,word,0) as test:E.cx(test,flag[0])


def install(cfraction):cfraction.armed=armed
