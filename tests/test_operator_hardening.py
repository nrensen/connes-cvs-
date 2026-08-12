"""Identity and precision tests for v0.3 operator hardening."""

from __future__ import annotations

import mpmath as mp
import pytest


def _legacy_prime_powers(c):
    is_pp = [False] * (c + 1)
    is_prime = [True] * (c + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, c + 1):
        if is_prime[i]:
            for j in range(i * i, c + 1, i):
                is_prime[j] = False
    primes = [i for i in range(2, c + 1) if is_prime[i]]
    for p in primes:
        pk = p
        while pk <= c:
            is_pp[pk] = True
            if pk > c // p:
                break
            pk *= p
    result = []
    for n in range(2, c + 1):
        if is_pp[n]:
            for p in primes:
                k = n
                while k % p == 0:
                    k //= p
                if k == 1:
                    lam = mp.log(p)
                    break
            result.append((n, mp.log(n), lam / mp.sqrt(n)))
    return result, primes


@pytest.mark.parametrize("c", [2, 13, 100])
def test_fast_prime_power_enumerator_is_raw_identical(c):
    from connes_cvs.operator import prime_powers_up_to

    with mp.workdps(80):
        expected, expected_primes = _legacy_prime_powers(c)
        actual, actual_primes = prime_powers_up_to(c)
    assert actual_primes == expected_primes
    assert [row[0] for row in actual] == [row[0] for row in expected]
    for actual_row, expected_row in zip(actual, expected):
        assert actual_row[1]._mpf_ == expected_row[1]._mpf_
        assert actual_row[2]._mpf_ == expected_row[2]._mpf_


def test_stable_b_taylor_is_high_precision():
    from connes_cvs.kernels import stable_B

    with mp.workdps(220):
        L = mp.log(13)
        beta = mp.mpf("9.999e-6") / L
        got = mp.im(stable_B(beta, L))
        reference = (beta * L - mp.sin(beta * L)) / (L * beta * beta)
        assert abs(got - reference) < mp.mpf("1e-210")


def test_stable_b_matches_triangular_weight_integral_away_from_zero():
    from connes_cvs.kernels import stable_B

    with mp.workdps(100):
        L = mp.log(13)
        beta = 1 / L  # beta*L = 1 exercises the non-series formula.
        reference = mp.quad(
            lambda y: (1 - y / L) * mp.exp(mp.j * beta * y),
            [0, L],
        )
        assert abs(stable_B(beta, L) - reference) < mp.mpf("1e-95")


@pytest.mark.parametrize("dps", [30, 80, 200])
@pytest.mark.parametrize("scaled_beta", ["0.000009999", "0.00001", "0.000010001", "0.1", "0.999", "1.001"])
def test_stable_b_retains_working_precision_across_formula_transitions(
    dps, scaled_beta
):
    from connes_cvs.kernels import stable_B

    with mp.workdps(dps + 50):
        L_reference = mp.log(13)
        beta_reference = mp.mpf(scaled_beta) / L_reference
        reference = mp.quad(
            lambda y: (1 - y / L_reference) * mp.exp(mp.j * beta_reference * y),
            [0, L_reference],
        )
    with mp.workdps(dps):
        L = mp.log(13)
        beta = mp.mpf(scaled_beta) / L
        got = stable_B(beta, L)
        # Inputs are rounded at the tested precision, so compare against the
        # same rounded inputs recomputed with guard digits.
        with mp.extradps(50):
            L_guard = mp.mpf(L)
            beta_guard = mp.mpf(beta)
            guarded_reference = mp.quad(
                lambda y: (1 - y / L_guard) * mp.exp(mp.j * beta_guard * y),
                [0, L_guard],
            )
        del reference
        scale = max(abs(guarded_reference), mp.mpf(1))
        assert abs(got - guarded_reference) / scale < mp.power(10, -(dps - 8))


def test_archimedean_cache_isolated_across_x_sign_without_manual_clear():
    from connes_cvs import operator

    with mp.workdps(30):
        L = mp.log(13)
        operator._hplus_cache_clear()
        positive = operator.psi_arch(1, L, 20, 30)
        negative_without_clear = operator.psi_arch(-1, L, 20, 30)
        operator._hplus_cache_clear()
        negative_fresh = operator.psi_arch(-1, L, 20, 30)
        expected_negative = -positive
    assert negative_without_clear._mpf_ == negative_fresh._mpf_
    assert negative_without_clear._mpf_ == expected_negative._mpf_


def test_hplus_cache_isolated_across_mpmath_precision(monkeypatch):
    from connes_cvs import operator

    monkeypatch.setattr(operator, "h_plus", operator._h_plus_mpmath)
    operator._hplus_cache_clear()
    with mp.workdps(30):
        low = operator._h_plus_cached(mp.mpf(1), 30)
    with mp.workdps(80):
        high_after_low = operator._h_plus_cached(mp.mpf(1), 80)
        operator._hplus_cache_clear()
        high_fresh = operator._h_plus_cached(mp.mpf(1), 80)
    assert high_after_low._mpf_ == high_fresh._mpf_
    assert high_after_low._mpf_ != low._mpf_


def test_hplus_cache_isolated_across_flint_precision_when_available():
    from connes_cvs import operator

    if not operator.HAS_FLINT:
        pytest.skip("requires python-flint")
    old_prec = operator.flint_ctx.prec
    try:
        operator._hplus_cache_clear()
        with mp.workdps(50):
            operator.flint_ctx.prec = 170
            operator._h_plus_cached(mp.mpf(1), 50)
            operator.flint_ctx.prec = 240
            high_after_low = operator._h_plus_cached(mp.mpf(1), 50)
            operator._hplus_cache_clear()
            high_fresh = operator._h_plus_cached(mp.mpf(1), 50)
        assert high_after_low._mpf_ == high_fresh._mpf_
    finally:
        operator.flint_ctx.prec = old_prec
        operator._hplus_cache_clear()


@pytest.mark.parametrize("function_name", ["stable_A", "stable_B"])
@pytest.mark.parametrize("bad_beta", [mp.nan, mp.inf, -mp.inf])
def test_public_scalar_kernels_reject_nonfinite_beta(function_name, bad_beta):
    from connes_cvs import kernels

    function = getattr(kernels, function_name)
    with pytest.raises(ValueError, match="beta must be finite"):
        function(bad_beta, mp.log(13))


@pytest.mark.parametrize(
    ("function_name", "args"),
    [
        ("stable_A", (0,)),
        ("stable_B", (0,)),
        ("S_hat_x", (0, 0)),
        ("dS_hat_x_dx", (0, 0)),
    ],
)
@pytest.mark.parametrize("bad_L", [mp.nan, mp.inf, -mp.inf, 0, -1])
def test_public_kernels_reject_invalid_cutoff(function_name, args, bad_L):
    from connes_cvs import kernels

    function = getattr(kernels, function_name)
    message = "L must be finite" if not mp.isfinite(bad_L) else "L must be greater than zero"
    with pytest.raises(ValueError, match=message):
        function(*args, bad_L)


@pytest.mark.parametrize("function_name", ["S_hat_x", "dS_hat_x_dx"])
@pytest.mark.parametrize("parameter", ["tau", "x"])
@pytest.mark.parametrize("bad_value", [mp.nan, mp.inf, -mp.inf])
def test_composite_kernels_reject_nonfinite_arguments(function_name, parameter, bad_value):
    from connes_cvs import kernels

    values = {"tau": mp.mpf(0), "x": mp.mpf(0)}
    values[parameter] = bad_value
    function = getattr(kernels, function_name)
    with pytest.raises(ValueError, match=rf"{parameter} must be finite"):
        function(values["tau"], values["x"], mp.log(13))


@pytest.mark.parametrize("function_name", ["stable_A", "stable_B"])
def test_public_scalar_kernels_reject_nonreal_input(function_name):
    from connes_cvs import kernels

    function = getattr(kernels, function_name)
    with pytest.raises(TypeError, match="beta must be a finite real number"):
        function(mp.mpc(1, 0), mp.log(13))


def test_stable_b_private_series_has_finite_termination_backstop():
    from connes_cvs.kernels import _stable_b_imag

    with pytest.raises(ArithmeticError, match="did not converge"):
        _stable_b_imag(mp.nan, mp.mpf(1), mp.nan, mp.nan)


@pytest.mark.parametrize("dps", [30, 50, 80])
@pytest.mark.parametrize("n", [0, 1, 3])
def test_psi_index_parity_is_raw_identical(dps, n):
    from connes_cvs import operator

    mp.mp.dps = dps
    if operator.HAS_FLINT:
        operator.flint_ctx.prec = int(dps * 3.5)
        operator.flint_ctx.threads = 1
    L = mp.log(13)
    prime_data, _ = operator.prime_powers_up_to(13)
    positive = operator._compute_psi_pair(n, L, 20, dps, prime_data)
    negative = operator._compute_psi_pair(-n, L, 20, dps, prime_data)
    assert negative[0]._mpf_ == (-positive[0])._mpf_
    assert negative[1]._mpf_ == positive[1]._mpf_


def test_upper_triangle_assembly_matches_legacy_raw_matrix():
    from connes_cvs.operator import _compute_psi_pair, prime_powers_up_to

    dps, c, N, T = 50, 13, 3, 20
    mp.mp.dps = dps
    L = mp.log(c)
    prime_data, _ = prime_powers_up_to(c)
    psi = {}
    psi_d = {}
    for n in range(-N, N + 1):
        psi[n], psi_d[n] = _compute_psi_pair(n, L, T, dps, prime_data)
    dimension = 2 * N + 1
    legacy = mp.matrix(dimension, dimension)
    direct = mp.matrix(dimension, dimension)
    for i in range(dimension):
        m = i - N
        for j in range(dimension):
            n = j - N
            legacy[i, j] = psi_d[n] if m == n else (psi[m] - psi[n]) / (m - n)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            avg = (legacy[i, j] + legacy[j, i]) / 2
            legacy[i, j] = legacy[j, i] = avg
    for i in range(dimension):
        m = i - N
        for j in range(i, dimension):
            n = j - N
            value = psi_d[n] if m == n else (psi[m] - psi[n]) / (m - n)
            direct[i, j] = value
            direct[j, i] = value
    for i in range(dimension):
        for j in range(dimension):
            assert direct[i, j]._mpf_ == legacy[i, j]._mpf_


def test_ground_state_rejects_symmetric_matrix_that_breaks_even_sector():
    from connes_cvs import compute_ground_state

    Q = mp.diag([1, 2, 3])
    with pytest.raises(ValueError, match="exactly centrosymmetric"):
        compute_ground_state(Q)


def test_built_matrix_is_raw_exactly_centrosymmetric_and_is_accepted():
    from connes_cvs import build_galerkin_matrix, compute_ground_state

    Q = build_galerkin_matrix(c=13, N=2, T=20, dps=30)
    dimension = Q.rows
    for i in range(dimension):
        for j in range(dimension):
            assert Q[i, j]._mpf_ == Q[dimension - 1 - i, dimension - 1 - j]._mpf_

    eigenvalue, eigenvector = compute_ground_state(Q)
    residual = Q * eigenvector - eigenvalue * eigenvector
    assert mp.norm(residual) <= mp.mpf("1e-28")
