"""J-fraction gap decoder with the merged two-slot extraction.

Identical to jfraction_gap_decoder except that:
  * moment_to_fraction and its inverse use one unit-series reciprocal per two
    slots (jfraction_merged);
  * the Forney scan shares one 12-qubit word between the locator value and
    the inverse of its derivative, so no constant-one register is needed
    (lifetime idea credited to anazemcev's public compact DQI notebook);
  * optionally, manwithacat's classical mask projection zeroes 70 masks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE / "keyeq"), str(HERE.parent)]
from jfraction_decoder import StaticArithmeticBuilder, build, km
from jfraction_gap_decoder import GapJFractionDecoder, convert_gap
from jfraction_merged import FlagBank, convert_banked, extract_merged
from compact_euclid import _EnabledZeroTest

PROJECT = os.environ.get("VQCC_DQI_PROJECTION", "1") == "1"
SHARED_FORNEY = os.environ.get("VQCC_DQI_FORNEY", "shared") == "shared"
BANKED = os.environ.get("VQCC_DQI_BANK", "1") == "1"


def forney_shared(b, numerator, sigma_tail, outcomes, *, field):
    """Forney phase scan; sigma(x) and 1/(x sigma'(x)) share the word `work`.

    numerator holds Omega[0..34], sigma_tail holds sigma[1..35] (sigma[0] = 1).
    Omega[0]'s word accumulates Omega(x) and sigma[1]'s word accumulates
    x sigma'(x); after the root test the value word is cleared and reused
    for the controlled inverse, then rebuilt to erase the retained AND chain.
    """
    degree = field.degree
    assert len(numerator) == len(sigma_tail) == 35 and degree == 12
    work = b.alloc_qubits(degree)
    flag = b.alloc_qubits(1)[0]
    enabled = b.alloc_bits(1)
    zero_test = _EnabledZeroTest(b, work, True)
    alpha = field.primitive_element
    order = (1 << degree) - 1
    omega, derivative = numerator[0], sigma_tail[0]

    def accumulate_omega():
        for word in numerator[1:]:
            b.cx(word, omega)

    def accumulate_derivative():
        for i in range(2, 35, 2):
            b.cx(sigma_tail[i], derivative)

    def accumulate_sigma():
        b.x(work[0])
        for i in range(1, 35, 2):
            b.cx(sigma_tail[i], work)
        b.cx(derivative, work)

    positions = min(order, len(outcomes))
    for j in range(positions):
        b.bit_store(enabled[0], True, control=outcomes[j])
        b.push_condition(enabled[0])
        accumulate_omega()
        accumulate_derivative()
        accumulate_sigma()
        zero_test.flip_if_equal(work, 0, target=flag)
        accumulate_sigma()
        _, scaffold = b.init_controlled_inverse(derivative, flag, target=work)
        x = field.pow(alpha, -j % order)
        b.gf2_imul(x, target=omega, field=field)
        b.gf2_phase_by_product(omega, work, mask=outcomes[j], field=field)
        b.gf2_idiv(x, target=omega, field=field)
        b.del_controlled_inverse(derivative, flag, target=work, scaffold=scaffold)
        b.free(scaffold)
        accumulate_sigma()
        zero_test.flip_if_equal(work, 0, target=flag)
        accumulate_sigma()
        accumulate_derivative()
        accumulate_omega()
        b.pop_condition()
        b.bit_store(enabled[0], False)
        if j + 1 < positions:
            for i, word in enumerate(sigma_tail, 1):
                b.gf2_imul(field.pow(alpha, -i % order), target=word, field=field)
            for i, word in enumerate(numerator[1:], 1):
                b.gf2_imul(field.pow(alpha, -i % order), target=word, field=field)
    if positions > 1:
        for i, word in enumerate(sigma_tail, 1):
            b.gf2_imul(field.pow(alpha, i * (positions - 1) % order), target=word, field=field)
        for i, word in enumerate(numerator[1:], 1):
            b.gf2_imul(field.pow(alpha, i * (positions - 1) % order), target=word, field=field)
    assert zero_test.chain is None and zero_test.calls == 2 * positions
    b.free(work + km.array([flag]) + enabled)


class MergedJFractionDecoder(GapJFractionDecoder):
    def build_syndrome_decoding_mbuc(self, builder, *, syndromes, error_mbuc, n_iters):
        b, words = builder, list(syndromes)
        if PROJECT:
            from mask_projection import projected_record
            error_mbuc = projected_record(b, list(syndromes), error_mbuc)
        b.stage_records = []
        flags = FlagBank(b, 35) if BANKED else b.alloc_qubits(35)
        convert = convert_banked if BANKED else convert_gap
        self._record(b, "initialization")
        extract_merged(b, words, flags)
        self._record(b, "moment_to_fraction")
        numerator, denominator = convert(b, words, flags)
        assert len(numerator) == 35 and len(denominator) == 36 and denominator[0] == 1
        self._record(b, "fraction_to_polynomials")
        if SHARED_FORNEY:
            forney_shared(b, numerator, denominator[1:], error_mbuc, field=self.field)
        else:
            from compact_euclid import CompactEuclidDecoder
            one = b.alloc_qubits(12)
            b.x(one[0])
            denominator[0] = one
            CompactEuclidDecoder(70, self.field)._forney(b, numerator, denominator, True, error_mbuc)
            b.x(one[0])
            b.free(one)
        self._record(b, "phase_correction")
        convert(b, words, flags, reverse=True)
        self._record(b, "polynomials_to_fraction")
        extract_merged(b, words, flags, reverse=True)
        if BANKED:
            assert flags.empty()
        else:
            b.free(flags)
        self._record(b, "fraction_to_moments")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    builder = StaticArithmeticBuilder()
    circuit = build(builder, iterations=35, decoder_type=MergedJFractionDecoder)
    circuit.to_file(str(args.output), format="kmb")
    with args.output.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    stats = circuit.operation_statistics()
    print(json.dumps({"candidate": "continued fraction, merged two-slot extraction",
                      "projection": PROJECT,
                      "num_qubits": circuit.num_qubits, "num_bits": circuit.num_bits,
                      "toffoli_count": stats.num_toffoli_like,
                      "t_count": stats.num_t_like, "num_ops": len(circuit),
                      "sha256": digest, "stages": builder.stage_records,
                      "build_seconds": time.monotonic() - start}, indent=2))


if __name__ == "__main__":
    main()
