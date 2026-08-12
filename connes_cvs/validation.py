"""Rigorous finite-matrix validation helpers.

The public function in this module certifies a supplied approximate eigenpair
of a *finite*, real-symmetric matrix.  It deliberately does not make a claim
about an infinite-dimensional operator or a truncation limit.
"""

import hashlib
from numbers import Integral
from typing import Optional

import mpmath as mp


_HASH_SCHEMA = b"connes-cvs.mpf-eigenpair.sha256.v1\0"
_LOG10_2_NUMERATOR = 30102999566398119521373889472449
_LOG10_2_UPPER_NUMERATOR = _LOG10_2_NUMERATOR + 1
_LOG10_2_DENOMINATOR = 10**32
_DISPLAY_DIGITS = 25
_MAX_EXACT_DISPLAY_BINARY_EXPONENT = 250_000


def _validated_dps(dps: Optional[int]) -> int:
    if dps is None:
        return max(50, int(mp.mp.dps))
    if isinstance(dps, bool) or not isinstance(dps, Integral):
        raise TypeError("dps must be an integer or None")
    dps_int = int(dps)
    if dps_int < 15:
        raise ValueError("dps must be at least 15")
    return dps_int


def _finite_mpf(value, name: str) -> mp.mpf:
    try:
        converted = mp.mpf(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not mp.isfinite(converted):
        raise ValueError(f"{name} must be finite")
    return converted


def _mpf_to_arb(value: mp.mpf, arb_type):
    """Convert an mpmath value to Arb without a decimal round trip."""

    sign, mantissa, exponent, _ = value._mpf_
    signed_mantissa = -int(mantissa) if sign else int(mantissa)
    if signed_mantissa == 0:
        return arb_type(0)
    return arb_type(signed_mantissa) * (arb_type(2) ** int(exponent))


def _signed_mantissa(value: mp.mpf):
    sign, mantissa, exponent, _ = value._mpf_
    signed = -int(mantissa) if sign else int(mantissa)
    return signed, int(exponent)


def _dyadic_metadata(mantissa: int, exponent: int):
    """Return a JSON-safe exact description of ``mantissa * 2**exponent``."""

    if mantissa == 0:
        sign = 0
        exponent = 0
    else:
        sign = -1 if mantissa < 0 else 1
    return {
        "sign": sign,
        "mantissa_hex": format(abs(mantissa), "x"),
        "exponent_base_2": int(exponent),
    }


def _dyadic_literal(mantissa: int, exponent: int) -> str:
    """Return an unambiguous exact (integer-mantissa) dyadic literal."""

    if mantissa == 0:
        return "0x0*2**0"
    sign = "-" if mantissa < 0 else ""
    return f"{sign}0x{abs(mantissa):x}*2**{int(exponent)}"


def _dyadic_num_den(mantissa: int, exponent: int):
    if exponent >= 0:
        return mantissa << exponent, 1
    return mantissa, 1 << (-exponent)


def _compare_positive_dyadic_to_power10(
    mantissa: int, exponent: int, decimal_exponent: int
) -> int:
    """Compare positive ``mantissa * 2**exponent`` with ``10**k``."""

    numerator, denominator = _dyadic_num_den(mantissa, exponent)
    if decimal_exponent >= 0:
        left = numerator
        right = denominator * (10**decimal_exponent)
    else:
        left = numerator * (10 ** (-decimal_exponent))
        right = denominator
    return (left > right) - (left < right)


def _positive_dyadic_decimal_exponent(mantissa: int, exponent: int) -> int:
    """Compute floor(log10(mantissa * 2**exponent)) exactly.

    A fixed-point binary-magnitude estimate locates the answer to within a
    tiny neighbourhood; exact integer comparisons then correct it.  No large
    integer is ever converted to a base-10 string.
    """

    binary_exponent = mantissa.bit_length() - 1 + exponent
    estimate = (
        binary_exponent * _LOG10_2_NUMERATOR // _LOG10_2_DENOMINATOR
    )
    while (
        _compare_positive_dyadic_to_power10(mantissa, exponent, estimate) < 0
    ):
        estimate -= 1
    while (
        _compare_positive_dyadic_to_power10(
            mantissa, exponent, estimate + 1
        )
        >= 0
    ):
        estimate += 1
    return estimate


def _positive_dyadic_outward_decimal(
    mantissa: int, exponent: int, digits: int = _DISPLAY_DIGITS
) -> str:
    """Return a short decimal guaranteed to be at least the exact dyadic.

    The conversion uses integer ceiling division, so this display remains an
    outward upper bound rather than a round-to-nearest approximation.
    """

    if mantissa < 0:
        raise ValueError("outward upper-bound formatting requires nonnegative data")
    if mantissa == 0:
        return "0"

    # Exact decimal ceiling conversion below is linear in the magnitude of
    # the base-10 exponent.  Avoid allowing an otherwise compact mpf such as
    # 2**1_000_000_000 to force a giant temporary integer.  The fallback is
    # still rigorous: the two fixed-point constants bracket log10(2), with
    # the lower one used for negative exponents where the inequality reverses.
    binary_exponent = mantissa.bit_length() - 1 + exponent
    if abs(binary_exponent) > _MAX_EXACT_DISPLAY_BINARY_EXPONENT:
        strict_upper_binary_exponent = binary_exponent + 1
        if strict_upper_binary_exponent >= 0:
            scaled = (
                strict_upper_binary_exponent * _LOG10_2_UPPER_NUMERATOR
            )
        else:
            scaled = strict_upper_binary_exponent * _LOG10_2_NUMERATOR
        decimal_exponent = -(
            (-scaled) // _LOG10_2_DENOMINATOR
        )
        coefficient = "1" + ("0" * (digits - 1))
        return f"{coefficient[0]}.{coefficient[1:]}e{decimal_exponent:+d}"

    decimal_exponent = _positive_dyadic_decimal_exponent(mantissa, exponent)
    scale = digits - 1 - decimal_exponent
    numerator, denominator = _dyadic_num_den(mantissa, exponent)
    if scale >= 0:
        numerator *= 10**scale
    else:
        denominator *= 10 ** (-scale)
    rounded_up = (numerator + denominator - 1) // denominator

    if rounded_up == 10**digits:
        rounded_up //= 10
        decimal_exponent += 1

    coefficient = str(rounded_up)
    if len(coefficient) != digits:
        raise ArithmeticError("internal outward-decimal formatting failure")
    return f"{coefficient[0]}.{coefficient[1:]}e{decimal_exponent:+d}"


def _update_integer(hasher, value: int) -> None:
    """Hash an integer with an unambiguous sign and binary length prefix."""

    integer = int(value)
    magnitude = abs(integer)
    payload = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    hasher.update(b"\x01" if integer < 0 else b"\x00")
    hasher.update(len(payload).to_bytes(8, "big"))
    hasher.update(payload)


def _update_mpf(hasher, value: mp.mpf) -> None:
    sign, mantissa, exponent, bitcount = value._mpf_
    hasher.update(b"\x01" if sign else b"\x00")
    _update_integer(hasher, int(mantissa))
    _update_integer(hasher, int(exponent))
    _update_integer(hasher, int(bitcount))


def _input_sha256(Q_values, vector_values, lambda_value: mp.mpf) -> str:
    """Hash the exact validated mpf payload using a canonical binary schema."""

    hasher = hashlib.sha256()
    hasher.update(_HASH_SCHEMA)
    hasher.update(b"Q\0")
    _update_integer(hasher, len(Q_values))
    _update_integer(hasher, len(Q_values[0]))
    for row in Q_values:
        for value in row:
            _update_mpf(hasher, value)
    hasher.update(b"v\0")
    _update_integer(hasher, len(vector_values))
    for value in vector_values:
        _update_mpf(hasher, value)
    hasher.update(b"lambda\0")
    _update_mpf(hasher, lambda_value)
    return hasher.hexdigest()


def arb_eigenpair_residual_bound(
    Q: mp.matrix,
    eigenvector: mp.matrix,
    eigenvalue,
    dps: Optional[int] = None,
) -> dict:
    """Certify a finite real-symmetric matrix eigenpair with Arb intervals.

    Returns a rigorous upper bound on

    ``||Q v - lambda v||_2 / ||v||_2``.

    For a finite real-symmetric matrix, the standard residual theorem then
    places at least one eigenvalue of ``Q`` within this bound of the exact
    dyadic center supplied in ``lambda_center_dyadic``.

    ``residual_bound`` is a short decimal rounded *upward*.  The exact Arb
    endpoint is in ``residual_bound_dyadic`` using a hexadecimal mantissa, so
    serialization remains reliable at precisions far above Python's optional
    decimal-integer conversion limit.  The same convention is used for the
    exact input center.

    This is a finite-matrix certificate only.  It makes no claim about an
    infinite-dimensional operator or a truncation limit.
    """

    dps_int = _validated_dps(dps)
    try:
        from flint import arb, ctx as flint_ctx
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "arb_eigenpair_residual_bound requires optional python-flint>=0.5.0"
        ) from exc

    if not isinstance(Q, mp.matrix):
        raise TypeError("Q must be an mpmath matrix")
    if Q.rows == 0 or Q.cols == 0 or Q.rows != Q.cols:
        raise ValueError("Q must be a non-empty square matrix")
    if not isinstance(eigenvector, mp.matrix):
        raise TypeError("eigenvector must be an mpmath matrix")
    if eigenvector.rows != Q.rows or eigenvector.cols != 1:
        raise ValueError("eigenvector must be an n-by-1 column matrix")

    n = Q.rows
    Q_values = [
        [_finite_mpf(Q[i, j], f"Q[{i},{j}]") for j in range(n)]
        for i in range(n)
    ]
    for i in range(n):
        for j in range(i + 1, n):
            if Q_values[i][j] != Q_values[j][i]:
                raise ValueError("Q must be exactly real-symmetric")

    vector_values = [
        _finite_mpf(eigenvector[i, 0], f"eigenvector[{i},0]")
        for i in range(n)
    ]
    if all(value == 0 for value in vector_values):
        raise ValueError("eigenvector must be nonzero")
    lambda_value = _finite_mpf(eigenvalue, "eigenvalue")

    input_sha256 = _input_sha256(Q_values, vector_values, lambda_value)
    lambda_mantissa, lambda_exponent = _signed_mantissa(lambda_value)

    old_prec = flint_ctx.prec
    old_threads = flint_ctx.threads
    certificate_bits = max(64, 4 * dps_int)
    try:
        flint_ctx.prec = certificate_bits
        flint_ctx.threads = 1

        Q_arb = [
            [_mpf_to_arb(Q_values[i][j], arb) for j in range(n)]
            for i in range(n)
        ]
        vector_arb = [_mpf_to_arb(value, arb) for value in vector_values]
        lambda_arb = _mpf_to_arb(lambda_value, arb)

        residual_sq_upper_sum = arb(0)
        norm_sq_enclosure = arb(0)
        for i in range(n):
            residual_i = sum(
                (Q_arb[i][j] * vector_arb[j] for j in range(n)), arb(0)
            ) - lambda_arb * vector_arb[i]
            residual_abs_upper = residual_i.abs_upper()
            residual_sq_upper_sum += residual_abs_upper * residual_abs_upper
            norm_sq_enclosure += vector_arb[i] * vector_arb[i]

        residual_sq_upper = residual_sq_upper_sum.upper()
        norm_sq_lower = norm_sq_enclosure.lower()
        if norm_sq_lower <= 0:
            raise ArithmeticError(
                "certificate precision did not establish a positive vector norm; "
                "increase dps"
            )

        if residual_sq_upper == 0:
            upper = arb(0)
        else:
            ratio_sq_upper = (residual_sq_upper / norm_sq_lower).upper()
            upper = ratio_sq_upper.sqrt().upper()
        bound_mantissa_raw, bound_exponent_raw = upper.man_exp()
        bound_mantissa = int(bound_mantissa_raw)
        bound_exponent = int(bound_exponent_raw)
        if bound_mantissa < 0:
            raise ArithmeticError("internal certificate endpoint was negative")
    finally:
        flint_ctx.prec = old_prec
        flint_ctx.threads = old_threads

    return {
        "schema": "connes-cvs.finite-matrix-eigenpair-certificate.v1",
        "residual_bound": _positive_dyadic_outward_decimal(
            bound_mantissa, bound_exponent
        ),
        "residual_bound_exact": _dyadic_literal(
            bound_mantissa, bound_exponent
        ),
        "residual_bound_dyadic": _dyadic_metadata(
            bound_mantissa, bound_exponent
        ),
        "residual_bound_display_digits": _DISPLAY_DIGITS,
        "lambda_center": _dyadic_literal(lambda_mantissa, lambda_exponent),
        "lambda_center_decimal_display": mp.nstr(
            lambda_value, min(50, dps_int)
        ),
        "lambda_center_dyadic": _dyadic_metadata(
            lambda_mantissa, lambda_exponent
        ),
        "matrix_dimension": n,
        "certificate_dps": dps_int,
        "certificate_bits": certificate_bits,
        "input_sha256": input_sha256,
        "input_hash_schema": _HASH_SCHEMA.rstrip(b"\0").decode("ascii"),
        "python_flint_version": getattr(__import__("flint"), "__version__", "unknown"),
        "finite_matrix_only": True,
        "rigorous_scope": "exact supplied finite symmetric mpf matrix only",
        "implication": (
            "For this finite real-symmetric matrix, at least one eigenvalue lies "
            "within residual_bound of the exact lambda_center_dyadic value. "
            "No infinite-dimensional or truncation-limit claim is implied."
        ),
    }
