"""Adversarial tests for the rigorous finite-matrix validator."""

import json
import sys
from decimal import Decimal
from fractions import Fraction
from random import Random

import mpmath as mp
import pytest

pytest.importorskip("flint")

from connes_cvs.validation import arb_eigenpair_residual_bound


def _mpf_from_dyadic(metadata):
    mantissa = int(metadata["mantissa_hex"], 16)
    if metadata["sign"] < 0:
        mantissa = -mantissa
    if mantissa == 0:
        return mp.mpf(0)
    sign = 1 if mantissa < 0 else 0
    magnitude = abs(mantissa)
    return mp.mpf(
        (sign, magnitude, metadata["exponent_base_2"], magnitude.bit_length())
    )


def _fraction_from_dyadic(metadata):
    mantissa = int(metadata["mantissa_hex"], 16) * metadata["sign"]
    exponent = metadata["exponent_base_2"]
    if exponent >= 0:
        return Fraction(mantissa << exponent, 1)
    return Fraction(mantissa, 1 << (-exponent))


def _relative_residual(Q, vector, eigenvalue):
    residual = Q * vector - eigenvalue * vector
    return mp.norm(residual) / mp.norm(vector)


def test_one_third_endpoint_is_outward_and_center_is_exact():
    with mp.workdps(100):
        one_third = mp.mpf(1) / 3
        result = arb_eigenpair_residual_bound(
            mp.matrix([[0]]), mp.matrix([[1]]), one_third, dps=80
        )

        exact_endpoint = _mpf_from_dyadic(result["residual_bound_dyadic"])
        exact_center = _mpf_from_dyadic(result["lambda_center_dyadic"])
        display_endpoint = mp.mpf(result["residual_bound"])

        assert exact_center == one_third
        assert exact_endpoint >= abs(one_third)
        assert display_endpoint >= exact_endpoint
        assert display_endpoint < abs(one_third) * (1 + mp.mpf("1e-20"))
        assert result["lambda_center"].startswith("0x")


def test_exact_eigenpair_has_exact_zero_residual():
    result = arb_eigenpair_residual_bound(
        mp.matrix([[1, 0], [0, 2]]),
        mp.matrix([[1], [0]]),
        mp.mpf(1),
        dps=30,
    )
    assert result["residual_bound"] == "0"
    assert result["residual_bound_exact"] == "0x0*2**0"
    assert result["residual_bound_dyadic"] == {
        "sign": 0,
        "mantissa_hex": "0",
        "exponent_base_2": 0,
    }


def test_random_exact_dyadic_residuals_are_below_certificate_endpoint():
    """Check soundness against exact rational arithmetic, not floating point."""

    rng = Random(20260812)
    for _ in range(75):
        n = rng.randrange(1, 5)
        denominator_shift = rng.randrange(0, 10)
        denominator = 1 << denominator_shift
        Q_fraction = [[Fraction(0) for _ in range(n)] for _ in range(n)]
        Q = mp.matrix(n)
        for i in range(n):
            for j in range(i, n):
                value = rng.randrange(-20, 21)
                Q_fraction[i][j] = Q_fraction[j][i] = Fraction(
                    value, denominator
                )
                Q[i, j] = Q[j, i] = mp.ldexp(mp.mpf(value), -denominator_shift)

        vector_integers = [rng.randrange(-10, 11) for _ in range(n)]
        if not any(vector_integers):
            vector_integers[0] = 1
        vector_fraction = [
            Fraction(value, denominator) for value in vector_integers
        ]
        vector = mp.matrix(
            [mp.ldexp(mp.mpf(value), -denominator_shift) for value in vector_integers]
        )
        lambda_integer = rng.randrange(-20, 21)
        lambda_fraction = Fraction(lambda_integer, denominator)
        eigenvalue = mp.ldexp(mp.mpf(lambda_integer), -denominator_shift)

        result = arb_eigenpair_residual_bound(
            Q, vector, eigenvalue, dps=15
        )
        endpoint = _fraction_from_dyadic(result["residual_bound_dyadic"])

        norm_sq = sum(value * value for value in vector_fraction)
        residual_sq = Fraction(0)
        for i in range(n):
            residual_i = sum(
                Q_fraction[i][j] * vector_fraction[j] for j in range(n)
            ) - lambda_fraction * vector_fraction[i]
            residual_sq += residual_i * residual_i
        assert endpoint * endpoint * norm_sq >= residual_sq


@pytest.mark.parametrize("certificate_dps", [15, 20, 30])
def test_low_certificate_precision_handles_high_precision_eigenpair(
    certificate_dps,
):
    with mp.workdps(120):
        Q = mp.matrix([[1, mp.mpf("0.3")], [mp.mpf("0.3"), 2]])
        eigenvalues, eigenvectors = mp.eigsy(Q)
        eigenvalue = eigenvalues[0]
        eigenvector = eigenvectors[:, 0]
        actual = _relative_residual(Q, eigenvector, eigenvalue)

        result = arb_eigenpair_residual_bound(
            Q, eigenvector, eigenvalue, dps=certificate_dps
        )
        exact_endpoint = _mpf_from_dyadic(result["residual_bound_dyadic"])
        display_endpoint = mp.mpf(result["residual_bound"])

        assert mp.isfinite(exact_endpoint)
        assert exact_endpoint >= actual
        assert display_endpoint >= exact_endpoint


def test_jordan_and_other_nonsymmetric_matrices_are_rejected():
    with pytest.raises(ValueError, match="real-symmetric"):
        arb_eigenpair_residual_bound(
            mp.matrix([[0, mp.mpf("1e100")], [0, 0]]),
            mp.matrix([[1], [0]]),
            0,
            dps=30,
        )
    with pytest.raises(ValueError, match="real-symmetric"):
        arb_eigenpair_residual_bound(
            mp.matrix([[1, 2], [mp.mpf("2.125"), 3]]),
            mp.matrix([[1], [1]]),
            1,
            dps=30,
        )


@pytest.mark.parametrize(
    "Q, vector, eigenvalue, error, match",
    [
        (mp.matrix(0, 0), mp.matrix(0, 1), 0, ValueError, "non-empty square"),
        (mp.matrix([[1, 2]]), mp.matrix([[1]]), 0, ValueError, "square"),
        (
            mp.matrix([[1, 0], [0, 2]]),
            mp.matrix([[1, 0]]),
            1,
            ValueError,
            "n-by-1",
        ),
        (
            mp.matrix([[1, 0], [0, 2]]),
            mp.matrix([[0], [0]]),
            1,
            ValueError,
            "nonzero",
        ),
        (mp.matrix([[mp.nan]]), mp.matrix([[1]]), 0, ValueError, "finite"),
        (mp.matrix([[1]]), mp.matrix([[mp.inf]]), 1, ValueError, "finite"),
        (mp.matrix([[1]]), mp.matrix([[1]]), mp.nan, ValueError, "finite"),
        (
            mp.matrix([[mp.mpc(1, 1)]]),
            mp.matrix([[1]]),
            1,
            TypeError,
            "real number",
        ),
    ],
)
def test_shapes_zero_nonfinite_and_nonreal_are_rejected(
    Q, vector, eigenvalue, error, match
):
    with pytest.raises(error, match=match):
        arb_eigenpair_residual_bound(Q, vector, eigenvalue, dps=30)


def test_input_hash_is_canonical_deterministic_and_precision_independent():
    with mp.workdps(100):
        Q = mp.matrix([[1, mp.mpf("0.125")], [mp.mpf("0.125"), 2]])
        vector = mp.matrix([[mp.mpf("0.75")], [mp.mpf("-0.5")]])
        eigenvalue = mp.mpf("1.25")

        first = arb_eigenpair_residual_bound(Q, vector, eigenvalue, dps=30)
        second = arb_eigenpair_residual_bound(Q, vector, eigenvalue, dps=80)
        changed = arb_eigenpair_residual_bound(
            mp.matrix([[1, mp.mpf("0.25")], [mp.mpf("0.25"), 2]]),
            vector,
            eigenvalue,
            dps=30,
        )

    assert first["input_sha256"] == second["input_sha256"]
    assert first["input_sha256"] == (
        "ab146b2625b0916b6f15af3a79d3ddb281d7d836e3488df40808a9614f01eb10"
    )
    assert first["input_sha256"] != changed["input_sha256"]
    assert first["input_hash_schema"].endswith("sha256.v1")
    assert len(first["input_sha256"]) == 64


def test_dps_5000_serializes_without_changing_process_integer_limit():
    get_limit = getattr(sys, "get_int_max_str_digits", None)
    before = get_limit() if get_limit is not None else None

    with mp.workdps(5000):
        one_third = mp.mpf(1) / 3
        result = arb_eigenpair_residual_bound(
            mp.matrix([[one_third]]),
            mp.matrix([[1]]),
            one_third,
            dps=5000,
        )

    encoded = json.dumps(result, sort_keys=True)
    after = get_limit() if get_limit is not None else None

    assert before == after
    assert result["residual_bound"] == "0"
    assert len(result["lambda_center_dyadic"]["mantissa_hex"]) > 4000
    assert int(
        result["lambda_center_dyadic"]["mantissa_hex"], 16
    ).bit_length() > 15000
    assert len(result["lambda_center_decimal_display"]) < 100
    assert "mantissa_hex" in encoded


@pytest.mark.parametrize("binary_exponent", [1_000_000, -1_000_000])
def test_extreme_exponent_uses_compact_sound_outward_display(binary_exponent):
    huge = mp.ldexp(mp.mpf(1), binary_exponent)
    result = arb_eigenpair_residual_bound(
        mp.matrix([[0]]), mp.matrix([[1]]), huge, dps=15
    )
    display_numerator, display_denominator = Decimal(
        result["residual_bound"]
    ).as_integer_ratio()
    display = Fraction(display_numerator, display_denominator)
    endpoint = _fraction_from_dyadic(result["residual_bound_dyadic"])

    assert display >= endpoint
    assert len(result["residual_bound"]) < 50
