"""
Stable kernel functions for the CvS Galerkin matrix.

These implement the Fourier-space kernels A(beta) and B(beta) used to
evaluate the archimedean piece of the Weil explicit formula without
catastrophic cancellation near beta = 0.

The key identity is:

    A(beta) = integral_0^L exp(i*beta*y) dy  (split into real/imag via trig)
    B(beta) = integral_0^L (1 - y/L) * exp(i*beta*y) dy  (derivative kernel)

Both are written in terms of sin(beta*L/2) to avoid the near-cancellation
in (exp(i*beta*L) - 1) / (i*beta) when beta*L is small.

References
----------
- Connes & van Suijlekom, arXiv:2511.23257, Proposition 4.1
- BUGFIX_eps_threshold.md in the research repository
"""

from __future__ import annotations

import mpmath as mp


def _finite_real(value: object, *, name: str) -> mp.mpf:
    """Convert one public kernel argument to a finite real ``mpf``.

    ``mpmath`` deliberately permits NaN and infinity.  Those values are not
    meaningful kernel parameters, and NaN used to make the adaptive series in
    :func:`_stable_b_imag` run forever.  Keep this validation at the public
    boundary so the fused operator kernel does not pay for it at every
    quadrature node.
    """
    try:
        converted = mp.mpf(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not mp.isfinite(converted):
        raise ValueError(f"{name} must be finite")
    return converted


def _positive_cutoff(value: object) -> mp.mpf:
    """Convert the logarithmic cutoff and require ``L > 0``."""
    converted = _finite_real(value, name="L")
    if converted <= 0:
        raise ValueError("L must be greater than zero")
    return converted


def _stable_b_imag(beta: mp.mpf, L: mp.mpf, bL: mp.mpf, sin_bL: mp.mpf) -> mp.mpf:
    """Return the imaginary part of ``B(beta)`` at the active precision.

    Near zero, direct subtraction in ``bL - sin(bL)`` loses digits.  The
    previous fixed four-term polynomial was adequate only to roughly 64
    decimal places near its threshold.  This series stops only when the next
    term rounds away at the caller's current mpmath precision, so the kernel
    remains meaningful for the package's high-precision workloads.
    """
    abs_bL = abs(bL)
    if abs_bL >= 1:
        return (bL - sin_bL) / (L * beta * beta)

    if abs_bL >= mp.mpf("1e-5"):
        # Re-evaluate the cancellation-prone numerator with guard digits.
        # Keeping this branch below |beta*L| = 1 makes fifteen extra decimal
        # digits enough to recover essentially the caller's full precision,
        # while avoiding a long Taylor loop across ordinary quadrature nodes.
        with mp.extradps(15):
            guarded_bL = beta * L
            guarded = (guarded_bL - mp.sin(guarded_bL)) / (L * beta * beta)
        return +guarded

    # L * (x - sin(x)) / x**2
    #   = L * sum_{k>=0} (-1)^k x**(2k+1) / (2k+3)!, x = beta*L.
    term = bL / 6
    total = term
    # A finite cap is a defensive backstop: equality-based convergence must
    # never become an unbounded loop if this private helper is called directly
    # with a non-finite value.  For |bL| < 1e-5 each term gains more than ten
    # decimal digits, so dps + 8 is deliberately far above the terms required.
    for k in range(max(32, mp.mp.dps + 8)):
        term *= -(bL * bL) / ((2 * k + 4) * (2 * k + 5))
        updated = total + term
        if updated == total:
            return L * total
        total = updated
    raise ArithmeticError("stable_B Taylor series did not converge")


def _stable_A_unchecked(beta: mp.mpf, L: mp.mpf) -> mp.mpc:
    """Evaluate ``A`` after the public boundary has validated its inputs."""
    if beta == 0:
        return mp.mpc(L, 0)
    bL = beta * L
    sin_half = mp.sin(bL / 2)
    sin_full = mp.sin(bL)
    return mp.mpc(sin_full / beta, 2 * sin_half * sin_half / beta)


def _stable_B_unchecked(beta: mp.mpf, L: mp.mpf) -> mp.mpc:
    """Evaluate ``B`` after the public boundary has validated its inputs."""
    if beta == 0:
        return mp.mpc(L / 2, 0)
    bL = beta * L
    sin_half = mp.sin(bL / 2)
    real_part = 2 * sin_half * sin_half / (L * beta * beta)
    imag_part = _stable_b_imag(beta, L, bL, mp.sin(bL))
    return mp.mpc(real_part, imag_part)


def stable_A(beta: mp.mpf, L: mp.mpf) -> mp.mpc:
    """
    Stable evaluation of A(beta) = integral_0^L exp(i*beta*y) dy.

    Returns A(beta) = sin(beta*L)/beta + 2i*sin^2(beta*L/2)/beta,
    which avoids cancellation in the naive (exp(i*beta*L) - 1)/(i*beta).

    Parameters
    ----------
    beta : mpmath.mpf
        Frequency parameter.
    L : mpmath.mpf
        Log-cutoff L = log(c).

    Returns
    -------
    mp.mpc
        Complex value of A(beta).
    """
    beta = _finite_real(beta, name="beta")
    L = _positive_cutoff(L)
    return _stable_A_unchecked(beta, L)


def stable_B(beta: mp.mpf, L: mp.mpf) -> mp.mpc:
    """
    Stable evaluation of
    B(beta) = integral_0^L (1 - y/L) * exp(i*beta*y) dy.

    The imaginary part uses a precision-adaptive Taylor series for
    ``|beta*L| < 1e-5`` to avoid loss of significance in
    ``(beta*L - sin(beta*L)) / (L * beta**2)``.

    Parameters
    ----------
    beta : mpmath.mpf
        Frequency parameter.
    L : mpmath.mpf
        Log-cutoff L = log(c).

    Returns
    -------
    mp.mpc
        Complex value of B(beta).
    """
    beta = _finite_real(beta, name="beta")
    L = _positive_cutoff(L)
    return _stable_B_unchecked(beta, L)


def S_hat_x(tau: mp.mpf, x: mp.mpf, L: mp.mpf) -> mp.mpc:
    """
    Compute S_hat(tau, x), the Fourier kernel for the archimedean integral.

    S_hat_x(tau, x) = sin(2*pi*x) * I_c(tau, x) - cos(2*pi*x) * I_s(tau, x)

    where I_c and I_s are symmetric/antisymmetric combinations of A(alpha +/- tau).

    Parameters
    ----------
    tau : mpmath.mpf
        Spectral parameter.
    x : mpmath.mpf
        Basis index (integer in practice).
    L : mpmath.mpf
        Log-cutoff L = log(c).

    Returns
    -------
    mp.mpc
        Complex value of S_hat(tau, x).
    """
    PI = mp.pi
    tau = _finite_real(tau, name="tau")
    x = _finite_real(x, name="x")
    L = _positive_cutoff(L)
    alpha = 2 * PI * x / L
    s2pi = mp.sin(2 * PI * x)
    c2pi = mp.cos(2 * PI * x)
    A_plus = _stable_A_unchecked(alpha - tau, L)
    A_minus = _stable_A_unchecked(-(alpha + tau), L)
    I_c = (A_plus + A_minus) / 2
    I_s = (A_plus - A_minus) / (2j)
    return s2pi * I_c - c2pi * I_s


def dS_hat_x_dx(tau: mp.mpf, x: mp.mpf, L: mp.mpf) -> mp.mpc:
    """
    Compute d/dx S_hat(tau, x), the derivative kernel for diagonal entries.

    Uses stable_B instead of stable_A for the derivative.

    Parameters
    ----------
    tau : mpmath.mpf
        Spectral parameter.
    x : mpmath.mpf
        Basis index.
    L : mpmath.mpf
        Log-cutoff L = log(c).

    Returns
    -------
    mp.mpc
        Complex value of dS_hat/dx(tau, x).
    """
    PI = mp.pi
    tau = _finite_real(tau, name="tau")
    x = _finite_real(x, name="x")
    L = _positive_cutoff(L)
    alpha = 2 * PI * x / L
    s2pi = mp.sin(2 * PI * x)
    c2pi = mp.cos(2 * PI * x)
    B_plus = _stable_B_unchecked(alpha - tau, L)
    B_minus = _stable_B_unchecked(-(alpha + tau), L)
    C = c2pi * (B_plus + B_minus) / 2 + s2pi * (B_plus - B_minus) / (2j)
    return 2 * PI * C
