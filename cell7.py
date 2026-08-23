# ============================================================
# CELL 7 — SYMBOLIC / ANALYTIC ARCHIMEDEAN DICTIONARY AUDIT
#
# Purpose:
#
#   Derive and numerically verify the chain
#
#       completed-zeta archimedean factor
#           ->
#       h_+(tau)
#           ->
#       Fourier response of the basis
#           ->
#       S_x(tau)
#           ->
#       psi_arch(x)
#
# This is deliberately a small diagnostic cell.
# It does NOT construct a large Galerkin matrix.
#
# ============================================================

import mpmath as mp

from cell import (
    compute_L,
    canonical_pairs,
    F_basis,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 80

c = 13
N = 4
T = 20

L = compute_L(c)

print("=" * 70)
print("CELL 7 — SYMBOLIC / ANALYTIC ARCHIMEDEAN DICTIONARY AUDIT")
print("=" * 70)
print()
print("Parameters:")
print(f"c = {c}")
print(f"N = {N}")
print(f"T = {T}")
print(f"dps = {mp.mp.dps}")
print(f"L = {mp.nstr(L, 60)}")
print(f"2*pi/L = {mp.nstr(2*mp.pi/L, 60)}")
print()


# ============================================================
# 1. ARCHIMEDEAN h_+(tau)
#
# Repository definition:
#
#   h_+(tau)
#       = Re psi(1/4 + i*tau/2) - log(pi)
#
# We compare this with the logarithmic derivative of
#
#   pi^(-s/2) Gamma(s/2)
#
# at
#
#   s = 1/2 + i*tau.
#
# ============================================================

def h_plus_repository(tau):
    tau = mp.mpf(tau)
    return (
        mp.re(mp.digamma(mp.mpf("0.25") + 0.5j * tau))
        - mp.log(mp.pi)
    )


def h_plus_from_completed_factor(tau):
    """
    Real part of

        d/ds log[ pi^(-s/2) Gamma(s/2) ]

    evaluated at s = 1/2 + i*tau.
    """
    tau = mp.mpf(tau)
    s = mp.mpf("0.5") + 1j * tau

    value = (
        -mp.log(mp.pi) / 2
        + mp.digamma(s / 2) / 2
    )

    return mp.re(value)


print("-" * 70)
print("1. ARCHIMEDEAN h_+(tau)")
print("-" * 70)
print()
print("We compare the repository h_+(tau) with")
print("Re d/ds log[pi^(-s/2) Gamma(s/2)]")
print("at s = 1/2 + i*tau.")
print()
print("IMPORTANT: this also exposes any factor-of-two convention.")
print()

for tau in [0, mp.mpf("0.5"), 1, 2, 5, 10]:
    hp = h_plus_repository(tau)
    hc = h_plus_from_completed_factor(tau)

    print(f"tau = {mp.nstr(tau, 20)}")
    print(f"  repository h_+       = {mp.nstr(hp, 50)}")
    print(f"  completed-factor Re  = {mp.nstr(hc, 50)}")
    print(f"  difference            = {mp.nstr(hp - hc, 30)}")
    print()


# ============================================================
# 2. THE FOURIER RESPONSE
#
# The canonical basis uses
#
#   e_k(t) = exp(2*pi*i*k*t/L)
#
# with t in [-L/2, L/2].
#
# For a single Fourier mode:
#
#   integral_{-L/2}^{L/2}
#       exp(2*pi*i*k*t/L) exp(-i*tau*t) dt
#
# =
#
#   exp(i*tau*L/2)
#   *
#   (exp(-i*tau*L)-1)
#   /
#   (i*(2*pi*k/L-tau)).
#
# Cell.py's F_basis is the real part of this quantity divided
# by sqrt(L), with the canonical +/-k combination.
#
# We independently evaluate the defining integral and compare.
# ============================================================

def direct_mode_integral(k, tau):
    """
    Direct numerical integral for the Fourier transform of

        exp(2*pi*i*k*t/L)

    over [-L/2, L/2].
    """
    k = int(k)
    tau = mp.mpf(tau)

    f = lambda t: (
        mp.exp(2j * mp.pi * k * t / L)
        * mp.exp(-1j * tau * t)
    )

    return mp.quad(f, [-L/2, L/2])


def analytic_mode_integral(k, tau):
    """
    Closed-form Fourier integral for one Fourier mode.
    """
    k = int(k)
    tau = mp.mpf(tau)

    a = 2 * mp.pi * k / L
    denom = a - tau

    if denom == 0:
        return mp.mpc(L)

    return (
        mp.exp(1j * tau * L / 2)
        * (mp.exp(-1j * tau * L) - 1)
        / (1j * denom)
    )


print("-" * 70)
print("2. FOURIER MODE INTEGRAL")
print("-" * 70)
print()

max_mode_error = mp.mpf("0")

for k in [0, 1, -1, 2]:
    for tau in [
        mp.mpf("0.3"),
        mp.mpf("1.0"),
        mp.mpf("2.0"),
        mp.mpf("2.4496332798546520107426711685"),
    ]:
        direct = direct_mode_integral(k, tau)
        analytic = analytic_mode_integral(k, tau)

        err = abs(direct - analytic)
        max_mode_error = max(max_mode_error, err)

        print(
            f"k={k:2d}, tau={mp.nstr(tau, 12)}  "
            f"|direct-analytic| = {mp.nstr(err, 20)}"
        )

print()
print(
    "Maximum Fourier-mode error =",
    mp.nstr(max_mode_error, 30),
)
print()


# ============================================================
# 3. CANONICAL BASIS RESPONSE
#
# Compare cell.py F_basis with an independent construction from
# the defining Fourier integrals.
#
# ============================================================

def direct_F_basis(k, tau):
    """
    Independent construction of the canonical real-even basis
    response.

    k = 0:
        e_0(t) = 1/sqrt(L)

    k > 0:
        canonical even basis corresponds to
            [e_{+k} + e_{-k}] / sqrt(2)
    """

    total = mp.mpc(0)

    for kk, ck in canonical_pairs(k):
        total += ck * direct_mode_integral(kk, tau)

    return mp.re(total / mp.sqrt(L))


print("-" * 70)
print("3. CANONICAL BASIS RESPONSE")
print("-" * 70)
print()

max_F_error = mp.mpf("0")

for k in range(0, N + 1):
    for tau in [
        mp.mpf("0.3"),
        mp.mpf("1.0"),
        mp.mpf("2.0"),
        mp.mpf("3.0"),
        mp.mpf("5.0"),
    ]:
        repo = F_basis(k, tau, L)
        direct = direct_F_basis(k, tau)

        err = abs(repo - direct)
        max_F_error = max(max_F_error, err)

        print(
            f"k={k:2d}, tau={mp.nstr(tau, 12)}"
        )
        print(
            f"  repository F = {mp.nstr(repo, 45)}"
        )
        print(
            f"  direct F     = {mp.nstr(direct, 45)}"
        )
        print(
            f"  |error|      = {mp.nstr(err, 20)}"
        )
        print()

print(
    "Maximum |F_repository - F_direct| =",
    mp.nstr(max_F_error, 30),
)
print()


# ============================================================
# 4. S_x(tau)
#
# For the real-even basis response, the relevant source kernel
# used by the current Cell-6 construction is
#
#   S_x(tau) = F_x(tau)^2
#
# for an individual canonical basis function.
#
# We compare the source construction with the direct Fourier
# interpretation.
#
# ============================================================

def S_from_F(k, tau):
    F = F_basis(k, tau, L)
    return F * F


def S_direct(k, tau):
    F = direct_F_basis(k, tau)
    return F * F


print("-" * 70)
print("4. SOURCE KERNEL S_x(tau)")
print("-" * 70)
print()

max_S_error = mp.mpf("0")

for k in range(0, N + 1):
    for tau in [
        mp.mpf("0.5"),
        mp.mpf("1.0"),
        mp.mpf("2.0"),
        mp.mpf("4.0"),
    ]:
        repo = S_from_F(k, tau)
        direct = S_direct(k, tau)

        err = abs(repo - direct)
        max_S_error = max(max_S_error, err)

        print(
            f"k={k:2d}, tau={mp.nstr(tau, 12)}  "
            f"|S_repo-S_direct| = {mp.nstr(err, 20)}"
        )

print()
print(
    "Maximum |S_repository - S_direct| =",
    mp.nstr(max_S_error, 30),
)
print()


# ============================================================
# 5. ARCHIMEDEAN QUADRATURE FROM THE DICTIONARY
#
# We now construct
#
#   A_k = C * integral h_+(tau) S_k(tau) d tau
#
# using the normalization appearing in operator.py.
#
# The purpose here is NOT yet to assert that this is the final
# CvS formula. It is to isolate the normalization and compare
# it against the repository's psi_arch contribution.
#
# We therefore print the raw integral first, followed by several
# candidate normalisations.
#
# ============================================================

def arch_raw(k, T):
    f = lambda tau: (
        h_plus_repository(tau)
        * S_from_F(k, tau)
    )

    # h_+ and S are even here.
    return 2 * mp.quad(f, [0, T])


print("-" * 70)
print("5. ARCHIMEDEAN RAW INTEGRAL")
print("-" * 70)
print()

for k in range(0, N + 1):
    raw = arch_raw(k, T)

    print(f"k = {k}")
    print(f"  integral h_+ S = {mp.nstr(raw, 50)}")
    print(f"  /(2*pi^2)      = {mp.nstr(raw/(2*mp.pi**2), 50)}")
    print(f"  /(pi^2)        = {mp.nstr(raw/(mp.pi**2), 50)}")
    print(f"  /(4*pi^2)      = {mp.nstr(raw/(4*mp.pi**2), 50)}")
    print()


# ============================================================
# 6. FACTOR-OF-TWO DIAGNOSTIC FOR h_+
#
# The most important question from the symbolic trace is whether
# the repository h_+ corresponds to
#
#   Re [d/ds log(pi^(-s/2) Gamma(s/2))]
#
# or twice that quantity.
#
# We therefore print the ratio directly.
#
# ============================================================

print("-" * 70)
print("6. h_+ NORMALISATION DIAGNOSTIC")
print("-" * 70)
print()

for tau in [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("5"),
]:
    hp = h_plus_repository(tau)
    hc = h_plus_from_completed_factor(tau)

    print(f"tau = {mp.nstr(tau, 15)}")
    print(f"  h_repository = {mp.nstr(hp, 45)}")
    print(f"  h_completed  = {mp.nstr(hc, 45)}")

    if hc != 0:
        print(
            f"  ratio        = {mp.nstr(hp/hc, 30)}"
        )
    print()


# ============================================================
# SUMMARY
# ============================================================

print("=" * 70)
print("CELL 7 SUMMARY")
print("=" * 70)
print()
print("Maximum Fourier-mode error:")
print(mp.nstr(max_mode_error, 30))
print()
print("Maximum canonical F_basis error:")
print(mp.nstr(max_F_error, 30))
print()
print("Maximum S-kernel error:")
print(mp.nstr(max_S_error, 30))
print()
print(
    "The key diagnostic is the relationship between the repository"
)
print(
    "h_+(tau) and the logarithmic derivative of"
)
print(
    "pi^(-s/2) Gamma(s/2) at s=1/2+i*tau."
)
print()
print("END CELL 7")
print("=" * 70)
