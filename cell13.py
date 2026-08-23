# ============================================================
# CELL 13 — CONTROLLED PARSEVAL / FOURIER-NORM AUDIT
#
# Purpose:
#
# Cell 12 found:
#
#     ||v||^2 = 1
#
# but naive numerical integration of
#
#     (1/2*pi) integral |H(tau)|^2 dtau
#
# appeared to converge toward approximately 0.686 rather
# than 1.
#
# This cell determines whether that discrepancy is:
#
#   (A) quadrature failure,
#   (B) finite-cutoff tail,
#   (C) a normalization error,
#   (D) an error in the closed transform.
#
# We use:
#
#   1. the existing closed transform;
#   2. increasingly fine piecewise quadrature;
#   3. symmetric frequency intervals;
#   4. an independent direct transform at selected points;
#   5. the exact coefficient-space norm.
#
# ============================================================

import mpmath as mp

from cell import (
    compute_L,
    canonical_to_full,
)

from connes_cvs import (
    build_galerkin_matrix,
    compute_ground_state,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 60

c = 13
N = 6
T = 40

L = compute_L(c)
omega = 2 * mp.pi / L


print("=" * 70)
print("CELL 13 — CONTROLLED PARSEVAL / FOURIER-NORM AUDIT")
print("=" * 70)
print()

print("Parameters:")
print(f"c = {c}")
print(f"N = {N}")
print(f"T = {T}")
print(f"dps = {mp.mp.dps}")
print(f"L = {mp.nstr(L, 50)}")
print(f"omega = 2*pi/L = {mp.nstr(omega, 50)}")
print()


# ============================================================
# 1. BUILD GALERKIN GROUND STATE
# ============================================================

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=mp.mp.dps,
)

lambda_min, v_star = compute_ground_state(Q)

target_norm = mp.fdot(v_star, v_star)


print("-" * 70)
print("1. EXACT FINITE-DIMENSIONAL NORM")
print("-" * 70)
print()

print("lambda_min =")
print(mp.nstr(lambda_min, 50))
print()

print("||v_star||^2 =")
print(mp.nstr(target_norm, 50))
print()

print("v^T Q v - lambda_min =")
print(
    mp.nstr(
        mp.fdot(v_star, Q * v_star) - lambda_min,
        30,
    )
)

print()


# ============================================================
# 2. FOURIER COEFFICIENTS
# ============================================================

def full_coefficients(v):
    u = canonical_to_full(v, N)

    return {
        k: mp.mpf(u[k + N])
        for k in range(-N, N + 1)
    }


u = full_coefficients(v_star)


print("-" * 70)
print("2. FULL FOURIER COEFFICIENT NORM")
print("-" * 70)
print()

full_norm = mp.mpf("0")

for k in range(-N, N + 1):
    full_norm += abs(u[k]) ** 2

print("sum |u_k|^2 =")
print(mp.nstr(full_norm, 50))
print()

print("|sum |u_k|^2 - ||v||^2| =")
print(
    mp.nstr(
        abs(full_norm - target_norm),
        30,
    )
)

print()


# ============================================================
# 3. CLOSED CENTRED FOURIER TRANSFORM
#
# This is the same expression validated in Cells 10-12.
# ============================================================

def g_repo(k, tau):
    tau = mp.mpf(tau)

    denom = omega * k - tau

    if denom == 0:
        return L

    return (
        mp.exp(-1j * tau * L) - 1
    ) / (
        1j * denom
    )


def H_closed(tau):

    tau = mp.mpf(tau)

    total = mp.mpc("0")

    for k in range(-N, N + 1):
        total += u[k] * g_repo(k, tau)

    total /= mp.sqrt(L)

    return mp.exp(1j * tau * L / 2) * total


def H_abs2(tau):
    H = H_closed(tau)
    return mp.re(H * mp.conj(H))


# ============================================================
# 4. DIRECT FOURIER TRANSFORM
#
# Independent check at selected points.
# ============================================================

def phi(k, t):

    if k == 0:
        return 1 / mp.sqrt(L)

    return (
        mp.sqrt(2 / L)
        * mp.cos(2 * mp.pi * k * t / L)
    )


def f_direct(t):

    total = mp.mpf("0")

    for k in range(N + 1):
        total += v_star[k] * phi(k, t)

    return total


def H_direct(tau):

    tau = mp.mpf(tau)

    return mp.quad(
        lambda x:
            f_direct(x + L / 2)
            * mp.exp(-1j * tau * x),
        [-L / 2, L / 2],
    )


print("-" * 70)
print("3. DIRECT / CLOSED FOURIER CHECK")
print("-" * 70)
print()

test_taus = [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    omega,
    mp.mpf("3"),
    mp.mpf("5"),
    mp.mpf("7.5"),
    mp.mpf("10"),
]

max_err = mp.mpf("0")

for tau in test_taus:

    hd = H_direct(tau)
    hc = H_closed(tau)

    err = abs(hd - hc)

    max_err = max(max_err, err)

    print("tau =", mp.nstr(tau, 20))
    print("  direct =", mp.nstr(hd, 35))
    print("  closed =", mp.nstr(hc, 35))
    print("  error  =", mp.nstr(err, 20))
    print()

print(
    "Maximum direct/closed error =",
    mp.nstr(max_err, 30),
)

print()


# ============================================================
# 5. NAIVE QUADRATURE
#
# Reproduce the Cell-12 method exactly, but only for R=20.
# This provides a baseline.
# ============================================================

print("-" * 70)
print("4. NAIVE CELL-12 QUADRATURE")
print("-" * 70)
print()

R = mp.mpf("20")

naive = mp.quad(
    lambda tau: H_abs2(tau),
    [-R, 0, R],
)

naive_norm = naive / (2 * mp.pi)

print("R =", R)
print("naive spectral norm =")
print(mp.nstr(naive_norm, 40))

print("naive error =")
print(
    mp.nstr(
        abs(naive_norm - target_norm),
        30,
    )
)

print()


# ============================================================
# 6. PIECEWISE QUADRATURE
#
# Split into short intervals.
#
# We deliberately use a modest interval width first.
# If the naive result is a quadrature artefact, this should
# change it substantially.
# ============================================================

def piecewise_integral(R, width):

    R = mp.mpf(R)
    width = mp.mpf(width)

    points = [-R]

    x = -R

    while x < R:
        x_next = min(x + width, R)
        points.append(x_next)
        x = x_next

    total = mp.mpf("0")

    for a, b in zip(points[:-1], points[1:]):

        total += mp.quad(
            lambda tau: H_abs2(tau),
            [a, b],
        )

    return total


print("-" * 70)
print("5. PIECEWISE QUADRATURE")
print("-" * 70)
print()

R_values = [
    mp.mpf("10"),
    mp.mpf("20"),
    mp.mpf("40"),
]

width_values = [
    mp.mpf("2"),
    mp.mpf("1"),
]

for width in width_values:

    print("Interval width =", width)
    print()

    for R in R_values:

        integral = piecewise_integral(
            R,
            width,
        )

        spectral_norm = integral / (2 * mp.pi)

        error = abs(
            spectral_norm - target_norm
        )

        print("R =", mp.nstr(R, 10))
        print(
            "  spectral norm =",
            mp.nstr(spectral_norm, 40),
        )
        print(
            "  error         =",
            mp.nstr(error, 25),
        )

    print()


# ============================================================
# 7. FREQUENCY-ALIGNED BREAKPOINTS
#
# The transform has special structure at
#
#     tau = k * omega
#
# for integer k.
#
# Include those points explicitly in the quadrature mesh.
# ============================================================

def aligned_integral(R):

    R = mp.mpf(R)

    points = {-R, mp.mpf("0"), R}

    for k in range(-N, N + 1):

        tau = mp.mpf(k) * omega

        if -R < tau < R:
            points.add(tau)

    # Add a regular mesh of width 1.
    x = -R

    while x < R:
        points.add(x)
        x = min(x + 1, R)

    points = sorted(points)

    total = mp.mpf("0")

    for a, b in zip(points[:-1], points[1:]):

        total += mp.quad(
            lambda tau: H_abs2(tau),
            [a, b],
        )

    return total


print("-" * 70)
print("6. FREQUENCY-ALIGNED QUADRATURE")
print("-" * 70)
print()

for R in R_values:

    integral = aligned_integral(R)

    spectral_norm = integral / (2 * mp.pi)

    error = abs(
        spectral_norm - target_norm
    )

    print("R =", mp.nstr(R, 10))
    print(
        "  spectral norm =",
        mp.nstr(spectral_norm, 40),
    )
    print(
        "  error         =",
        mp.nstr(error, 25),
    )
    print()


# ============================================================
# 8. DIRECT TIME-DOMAIN NORM
#
# This is another independent check:
#
#     integral |f(x)|^2 dx
#
# should equal ||v||^2.
# ============================================================

print("-" * 70)
print("7. DIRECT TIME-DOMAIN NORM")
print("-" * 70)
print()

time_norm = mp.quad(
    lambda x:
        abs(f_direct(x + L / 2)) ** 2,
    [-L / 2, L / 2],
)

print("time-domain norm =")
print(mp.nstr(time_norm, 50))

print()

print("|time norm - coefficient norm| =")
print(
    mp.nstr(
        abs(time_norm - target_norm),
        30,
    )
)

print()


# ============================================================
# 9. PARSEVAL CONVENTION CHECK
#
# There is a possible normalization trap here, so explicitly
# report the raw integral and the values obtained under the
# two most common Fourier conventions.
# ============================================================

print("-" * 70)
print("8. NORMALISATION CHECK")
print("-" * 70)
print()

R = mp.mpf("20")

integral20 = aligned_integral(R)

print("Raw integral over [-20,20] =")
print(mp.nstr(integral20, 40))

print()

print("Raw integral / (2*pi) =")
print(
    mp.nstr(
        integral20 / (2 * mp.pi),
        40,
    )
)

print()

print("Raw integral / (2*pi)^2 =")
print(
    mp.nstr(
        integral20 / (2 * mp.pi) ** 2,
        40,
    )
)

print()

print("Raw integral =")
print(mp.nstr(integral20, 40))

print()

print("Target coefficient norm =")
print(mp.nstr(target_norm, 40))

print()


# ============================================================
# 10. SUMMARY
# ============================================================

print("=" * 70)
print("CELL 13 SUMMARY")
print("=" * 70)
print()

print(
    "Coefficient-space norm =",
    mp.nstr(target_norm, 30),
)

print(
    "Full Fourier coefficient norm =",
    mp.nstr(full_norm, 30),
)

print(
    "Direct time-domain norm =",
    mp.nstr(time_norm, 30),
)

print(
    "Maximum direct/closed Fourier error =",
    mp.nstr(max_err, 20),
)

print()

print(
    "The purpose of Cell 13 is to determine whether the"
)

print(
    "apparent Cell-12 Parseval discrepancy survives controlled"
)

print(
    "piecewise and frequency-aligned quadrature."
)

print()

print("=" * 70)
print("END CELL 13")
print("=" * 70)
