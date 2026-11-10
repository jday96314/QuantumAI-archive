"""Recover J-fraction gap flags locally instead of retaining a flag bank.

Research approximation: g_j is recovered as [c_j=0, a_j!=0, c_(j-1)!=0].
On ordinary and completed physical charts this equals the extraction flag.
A gap followed by a zero pending header and another gap is unsupported, as
are higher-order zero pivots. Guard HMR clears residual flag bits on those
exceptional charts; their error must be measured by the unchanged verifier.
This is not an all-input literal-inverse claim.

Base decoder: SomeJoint / Agniv Banerjee; merged extraction and lifetime
changes: Ahmed Berat Özer; mask projection: manwithacat; compact storage
ideas: anazemcev. Apache-2.0; see public/ahmed_v3/dqi_root/NOTICE.md/LICENSE.
This file adds local flag recovery and is a modified implementation.
"""
from contextlib import contextmanager
from pathlib import Path
import argparse
import hashlib
import inspect
import json
import sys
import time

HERE = Path(__file__).resolve().parent
BASE = HERE/'public/ahmed_v3/dqi_root/dqi'
sys.path[:0] = [str(BASE/'continued_fraction'), str(BASE/'continued_fraction/keyeq'),
               str(BASE), str((HERE/'../../../from_host/verifiers/dqi').resolve())]
from common import build, km
from static_divisions import StaticArithmeticBuilder, StaticDivisionDecoder
import jfraction_merged as merged
from jfraction_decoder import conversion_layout, scale_words
from jfraction_gap_decoder import rotate, conversion_core, toggle_header
from jfraction_merged_decoder import forney_shared


def recycled_pair():
    """Release the source-derived pre flag during the two scalar inversions."""
    source = inspect.getsource(merged._pair)
    marker = '    # ---- post-processing (A and B share it; B still shifted) ----'
    helper = '''    def recycle_pre(callback):
        nonlocal pre
        flip_pre()
        b.free(km.array([pre]))
        callback()
        pre = b.alloc_qubits(1)[0]
        flip_pre()

'''
    assert source.count(marker) == 1
    source = source.replace(marker,helper+marker)
    for old,new in [
        ('        window(forward=True)','        recycle_pre(lambda: window(forward=True))'),
        ('        window(forward=False)','        recycle_pre(lambda: window(forward=False))'),
        ('        _scale_by(b, act[0], act[2], divide=divide)',
         '        recycle_pre(lambda: _scale_by(b, act[0], act[2], divide=divide))')]:
        assert source.count(old) == 1
        source = source.replace(old,new)
    namespace = dict(merged.__dict__)
    namespace['km'] = km
    exec(source,namespace)
    (HERE/'parameters/recovered_gap_pair_build_source.py').write_text(source)
    return namespace['_pair']


pair = recycled_pair()


@contextmanager
def zero_chain(b, words, start=None):
    """Retain an AND of complemented source bits; erase it by measurement."""
    bits = [q for word in words for q in word]
    assert len(bits) >= 2
    chain = b.alloc_qubits(len(bits) - (1 if start is None else 0))
    b.x(km.array(bits))
    if start is None:
        b.init_and(bits[0], bits[1], chain[0])
        first = 2
    else:
        b.init_and(start, bits[0], chain[0])
        first = 1
    for index, source in enumerate(bits[first:], 1):
        b.init_and(chain[index-1], source, chain[index])
    yield chain[-1]
    for index in reversed(range(1, len(chain))):
        b.del_and(chain[index-1], bits[index+first-1], chain[index])
    if start is None:
        b.del_and(bits[0], bits[1], chain[0])
    else:
        b.del_and(start, bits[0], chain[0])
    b.x(km.array(bits))
    b.free(chain)


def xor_predicate(b, c, a, target, previous_c=None, pending=None):
    """XOR zero(c) AND nonzero(a), with the specified extra exclusion."""
    with zero_chain(b, [c]) as zc, zero_chain(b, [a], zc) as zca:
        if previous_c is None and pending is None:
            b.cx(zc, target)
            b.cx(zca, target)
            return
        predicate = b.alloc_qubits(1)[0]
        b.cx(zc, predicate)
        b.cx(zca, predicate)
        if previous_c is not None:
            with zero_chain(b, [previous_c]) as zp:
                b.x(zp)
                b.ccx(predicate, zp, target)
                b.x(zp)
        else:
            b.x(pending)
            b.ccx(predicate, pending, target)
            b.x(pending)
        b.cx(zca, predicate)
        b.cx(zc, predicate)
        b.free(km.array([predicate]))


def xor_head(b, words, j, target):
    xor_predicate(b, words[2*j], words[2*j+1], target,
                  previous_c=words[2*j-2] if j else None)


@contextmanager
def head_flag(b, words, j):
    if j < 0:
        yield False
        return
    flag = b.alloc_qubits(1)[0]
    xor_head(b, words, j, flag)
    yield flag
    xor_head(b, words, j, flag)
    b.free(km.array([flag]))


def guard_clear(b, flags):
    """Reset exceptional residual flags; any resulting error is benchmarked."""
    out = b.alloc_bits(len(flags))
    b.hmr(km.array(flags), out)
    b.bit_store(out, False)
    b.free(out)
    b.free(km.array(flags))


def extract(b, words, reverse=False):
    slots = len(words)//2
    steps = list(range(0, slots, 2))
    for j in reversed(steps) if reverse else steps:
        with head_flag(b, words, j-1) as previous:
            flags = list(b.alloc_qubits(min(2, slots-j)))
            if reverse:
                for k, flag in enumerate(flags):
                    xor_head(b, words, j+k, flag)
            if len(flags) == 2 and len(words)-2*j >= 6:
                pair(b, words[2*j:], previous, *flags, reverse=reverse)
            else:
                # A final two-slot window is handled by serial gap steps.
                for k in reversed(range(len(flags))) if reverse else range(len(flags)):
                    merged._single_slot(b, words, previous if k == 0 else flags[k-1], flags[k], j+k, reverse=reverse)
            if not reverse:
                for k, flag in enumerate(flags):
                    xor_head(b, words, j+k, flag)
            guard_clear(b, flags)


def convert(b, words, reverse=False):
    p, q, stages = conversion_layout(words)
    stages = [(len(words)//2-1-i, *stage) for i, stage in enumerate(stages)]
    for j, x, y, old_p, old_q in reversed(stages) if reverse else stages:
        with head_flag(b, words, j-1) as pending:
            gap = b.alloc_qubits(1)[0]
            new_p = [x] + old_q[1:]
            def numerator_flag():
                if len(new_p) > 1:
                    xor_predicate(b, new_p[0], new_p[1], gap,
                                  pending=None if pending is False else pending)
            if reverse:
                numerator_flag()
                rotate(b, new_p, gap, left=True)
                scale_words(b, x, old_q[1:], divide=True)
                toggle_header(b, y, old_q, pending)
            else:
                xor_head(b, words, j, gap)
                if j:
                    b.cswap(pending, x, y)
                rotate(b, [x,y]+old_q[1:], gap, left=True)
            conversion_core(b, y, old_p, old_q, pending)
            if reverse:
                rotate(b, [x,y]+old_q[1:], gap, left=False)
                if j:
                    b.cswap(pending, x, y)
                xor_head(b, words, j, gap)
            else:
                toggle_header(b, y, old_q, pending)
                scale_words(b, x, old_q[1:], divide=False)
                rotate(b, new_p, gap, left=False)
                numerator_flag()
            guard_clear(b, [gap])
    return p, q


class RecoveredGapDecoder(StaticDivisionDecoder):
    def build_syndrome_decoding_mbuc(self, builder, *, syndromes, error_mbuc, n_iters):
        b, words = builder, list(syndromes)
        from mask_projection import projected_record
        error_mbuc = projected_record(b, words, error_mbuc)
        b.stage_records = []
        self._record(b, 'initialization')
        extract(b, words)
        self._record(b, 'moment_to_fraction')
        numerator, denominator = convert(b, words)
        self._record(b, 'fraction_to_polynomials')
        forney_shared(b, numerator, denominator[1:], error_mbuc, field=self.field)
        self._record(b, 'phase_correction')
        convert(b, words, True)
        self._record(b, 'polynomials_to_fraction')
        extract(b, words, True)
        self._record(b, 'fraction_to_moments')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=70)
    parser.add_argument('--full', action='store_true')
    args = parser.parse_args()
    b = StaticArithmeticBuilder()
    start = time.time()
    if args.full:
        circuit = build(b, iterations=35, decoder_type=RecoveredGapDecoder)
        path = HERE/'candidates/recovered_gap_dqi.kmb'
    else:
        words = [b.create_quantum_register(12, name=f'S{i}') for i in range(args.n)]
        extract(b, words)
        p,q = convert(b, words)
        convert(b, words, True)
        extract(b, words, True)
        circuit = b.finish_circuit()
        path = HERE/f'candidates/recovered_gap_n{args.n}_both.kmb'
    circuit.to_file(str(path), format='kmb')
    with path.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    report = dict(status='BUILT_APPROXIMATION_REQUIRES_NATIVE_CHECKS', path=str(path),
                  Q=circuit.num_qubits, T=circuit.operation_statistics().num_toffoli_like,
                  sha256=sha, seconds=time.time()-start, stages=getattr(b,'stage_records',[]))
    (HERE/'results'/f'{path.stem}_build.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
