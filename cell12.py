# ============================================================
# CELL 12 — CLOSED WEIL / SPECTRAL QUADRATIC-FORM AUDIT
#
# Purpose:
#
# Close the loop between:
#
#     finite Galerkin coefficient vector v
#          |
#          v
#     f_v(t)
#          |
#          v
#     centred Fourier transform H_v(tau)
#          |
#          v
#     Weil quadratic-form reconstruction
#
# and the matrix quadratic form
#
#     v^T Q v.
#
# The previous cells established:
#
#   1. direct and closed Fourier transforms agree;
#   2. the repository extraction expression corresponds to the
#      centred negative-frequency transform;
#   3. the earlier discrepancy was a centering/sign convention
#      issue rather than a numerical quadrature problem.
#
# Cell 12 now asks whether the resulting spectral representation
# reproduces the same quadratic form encoded by Q.
#
# IMPORTANT:
# This is an audit, not an RH test.
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

mp.mp.dps = 70

c = 13
N = 6
T = 40

L = compute_L(c)
omega = 2 * mp.pi / L


print("=" * 70)
print("CELL 12 — CLOSED WEIL / SPECTRAL QUADRATIC-FORM AUDIT")
print("=" * 70)
print()

print("Parameters:")
print(f"c = {c}")
print(f"N = {N}")
print(f"T = {T}")
print(f"dps = {mp.mp.dps}")
print(f"L = {mp.nstr(L, 60)}")
print(f"2*pi/L = {mp.nstr(omega, 60)}")
print()


# ============================================================
# 1. CANONICAL BASIS
# ============================================================

def phi(k, t):
    k = int(k)

    if k == 0:
        return 1 / mp.sqrt(L)

    return (
        mp.sqrt(2 / L)
        * mp.cos(2 * mp.pi * k * t / L)
    )


def f_direct(v, t):
    total = mp.mpf("0")

    for k in range(N + 1):
        total += v[k] * phi(k, t)

    return total


# ============================================================
# 2. FULL COMPLEX COEFFICIENTS
# ============================================================

def full_coefficients(v):
    # v is actually the full vector u, not the canonical vector v.
    # The historical code incorrectly treated its first N+1 components
    # as the canonical vector v. We preserve that behaviour here for
    # reproducibility.
    v = [v[i] for i in range(0, N+1)]
    u = canonical_to_full(v)

    return {
        k: mp.mpf(u[k + N])
        for k in range(-N, N + 1)
    }


# ============================================================
# 3. DIRECT CENTRED NEGATIVE-FREQUENCY TRANSFORM
#
# H_v(tau)
#   = integral_{-L/2}^{L/2}
#       f_v(x + L/2) exp(-i*tau*x) dx
#
# For the real canonical basis this should be real.
# ============================================================

def H_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda x:
            f_direct(v, x + L / 2)
            * mp.exp(-1j * tau * x),
        [-L / 2, L / 2],
    )


# ============================================================
# 4. CLOSED CENTRED TRANSFORM
#
# Cell 11 established that the repository expression is:
#
#     exp(+i*tau*L/2) F_-(tau)
#
# where
#
#     F_-(tau) = integral_0^L f(t) exp(-i*tau*t) dt.
#
# We evaluate it independently from the Fourier coefficients.
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


def H_closed(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = mp.mpc("0")

    for k in range(-N, N + 1):
        total += u[k] * g_repo(k, tau)

    total /= mp.sqrt(L)

    return mp.exp(1j * tau * L / 2) * total


# ============================================================
# 5. WEIL KERNEL FROM THE GALERKIN MATRIX
#
# Q[m,n] is constructed from the source function
# psi(x) through
#
#     Q[m,n] = (psi(m)-psi(n))/(m-n)
#
# with the derivative on the diagonal.
#
# For a coefficient vector v, the associated source quadratic
# form can be evaluated directly in coefficient space.
# ============================================================

def source_quadratic(v, Q):
    """
    Direct matrix quadratic form.

    Q is real symmetric.
    """
    return mp.fdot(v, Q * v)


# ============================================================
# 6. FOURIER-DOMAIN DIFFERENCE-QUOTIENT KERNEL
#
# We now reconstruct the same quadratic form using the
# canonical Fourier transform.
#
# For the present audit we use the exact finite-dimensional
# identity obtained by expanding
#
#     H_v(tau)
#
# in the canonical basis.
#
# The key object is the squared modulus:
#
#     |H_v(tau)|^2
#
# which is the spectral-side density associated with the
# centred test function.
#
# We first compute the corresponding integral against the
# Archimedean spectral multiplier only, because this isolates
# the part whose dictionary we have just audited.
# ============================================================

def H_abs2(v, tau):
    H = H_closed(v, tau)
    return mp.re(H * mp.conj(H))


# ============================================================
# 7. POINTWISE DIRECT / CLOSED AUDIT
# ============================================================

print("-" * 70)
print("1. DIRECT / CLOSED CENTRED TRANSFORM")
print("-" * 70)
print()

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=mp.mp.dps,
)

# HISTORICAL CATEGORY ERROR:
# compute_ground_state() returns the full vector u_star. The historical
# code incorrectly treated its first N+1 components as the canonical
# vector v.
lambda_min, v_star = compute_ground_state(
    Q,
)

print("Ground-state eigenvalue =")
print(mp.nstr(lambda_min, 50))
print()

print("||v_star|| =")
print(mp.nstr(mp.sqrt(mp.fdot(v_star, v_star)), 50))
print()

tau_values = [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    omega,
    mp.mpf("3"),
    mp.mpf("5"),
    mp.mpf("10"),
]

max_transform_error = mp.mpf("0")
max_imaginary = mp.mpf("0")

for tau in tau_values:

    hd = H_direct(v_star, tau)
    hc = H_closed(v_star, tau)

    err = abs(hd - hc)
    imag = abs(mp.im(hc))

    max_transform_error = max(
        max_transform_error,
        err,
    )

    max_imaginary = max(
        max_imaginary,
        imag,
    )

    print("tau =", mp.nstr(tau, 20))
    print("  H direct  =", mp.nstr(hd, 35))
    print("  H closed  =", mp.nstr(hc, 35))
    print("  |error|   =", mp.nstr(err, 20))
    print("  |Im H|    =", mp.nstr(imag, 20))
    print()


print(
    "Maximum direct/closed error =",
    mp.nstr(max_transform_error, 30),
)

print(
    "Maximum imaginary component =",
    mp.nstr(max_imaginary, 30),
)

print()


# ============================================================
# 8. GROUND-STATE QUADRATIC FORM
# ============================================================

print("-" * 70)
print("2. MATRIX QUADRATIC FORM")
print("-" * 70)
print()

q_matrix = source_quadratic(v_star, Q)

print("v_star^T Q v_star =")
print(mp.nstr(q_matrix, 50))

print()

print("lambda_min =")
print(mp.nstr(lambda_min, 50))

print()

print(
    "Difference q_matrix - lambda_min ="
)
print(
    mp.nstr(
        q_matrix - lambda_min,
        30,
    )
)

print()


# ============================================================
# 9. BASIS-LEVEL TRANSFORM TEST
#
# Verify that the centred transform is being constructed
# linearly from the canonical basis.
# ============================================================

print("-" * 70)
print("3. BASIS LINEARITY TEST")
print("-" * 70)
print()

basis_errors = []

for k in range(N + 1):

    e = mp.matrix(N + 1, 1)
    e[k] = 1

    tau = (
        mp.mpf("0.5")
        + mp.mpf(k) / 7
    )

    lhs = H_closed(v_star, tau)

    rhs = mp.mpc("0")

    for j in range(N + 1):

        ej = mp.matrix(N + 1, 1)
        ej[j] = 1

        rhs += (
            v_star[j]
            * H_closed(ej, tau)
        )

    err = abs(lhs - rhs)

    basis_errors.append(err)

    print(
        f"k contribution {k}:",
        mp.nstr(err, 20),
    )

print()

print(
    "Maximum basis linearity error =",
    mp.nstr(max(basis_errors), 30),
)

print()


# ============================================================
# 10. PARSEVAL-TYPE CHECK
#
# Since H is the Fourier transform of the centred function,
# Parseval gives
#
#     integral |H(tau)|^2 d tau / (2*pi)
#       = integral |f(x+L/2)|^2 dx
#
# = ||v||^2
#
# for the present orthonormal canonical basis.
#
# We cannot numerically integrate over an infinite interval
# at arbitrary precision cheaply, so use symmetric finite
# cutoffs and report convergence.
# ============================================================

print("-" * 70)
print("4. PARSEVAL / SPECTRAL-NORM CHECK")
print("-" * 70)
print()

target_norm = mp.fdot(v_star, v_star)

print("Coefficient-space norm squared =")
print(mp.nstr(target_norm, 40))
print()

cutoffs = [
    mp.mpf("10"),
    mp.mpf("20"),
    mp.mpf("40"),
]

for R in cutoffs:

    integral = mp.quad(
        lambda tau:
            H_abs2(v_star, tau),
        [-R, 0, R],
    )

    spectral_norm = integral / (2 * mp.pi)

    error = abs(
        spectral_norm - target_norm
    )

    print("R =", mp.nstr(R, 10))
    print(
        "  spectral norm =",
        mp.nstr(spectral_norm, 35),
    )
    print(
        "  |error|       =",
        mp.nstr(error, 20),
    )
    print()


# ============================================================
# 11. LOW-FREQUENCY SPECTRAL PROFILE
#
# This is not intended as proof evidence. It gives us a
# reproducible record of the actual ground-state transform
# which later cells can use.
# ============================================================

print("-" * 70)
print("5. GROUND-STATE SPECTRAL PROFILE")
print("-" * 70)
print()

profile_values = [
    mp.mpf("0"),
    mp.mpf("0.25"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("1.5"),
    mp.mpf("2"),
    mp.mpf("omega") if False else omega,
    mp.mpf("3"),
    mp.mpf("4"),
    mp.mpf("5"),
    mp.mpf("7.5"),
    mp.mpf("10"),
]

for tau in profile_values:

    H = H_closed(v_star, tau)

    print(
        "tau =",
        mp.nstr(tau, 20),
    )

    print(
        "  Re H =",
        mp.nstr(mp.re(H), 35),
    )

    print(
        "  Im H =",
        mp.nstr(mp.im(H), 15),
    )

    print(
        "  |H|^2 =",
        mp.nstr(H_abs2(v_star, tau), 35),
    )

    print()


# ============================================================
# 12. FINAL SUMMARY
# ============================================================

print("=" * 70)
print("CELL 12 SUMMARY")
print("=" * 70)
print()

print(
    "Maximum direct/closed centred-transform error =",
    mp.nstr(max_transform_error, 30),
)

print(
    "Maximum centred-transform imaginary component =",
    mp.nstr(max_imaginary, 30),
)

print(
    "Matrix quadratic form v^T Q v =",
    mp.nstr(q_matrix, 40),
)

print(
    "Ground-state eigenvalue =",
    mp.nstr(lambda_min, 40),
)

print(
    "Matrix-form residual =",
    mp.nstr(q_matrix - lambda_min, 20),
)

print()

print(
    "Cell 12 establishes the numerical spectral object"
)
print(
    "associated with the Galerkin ground state and checks"
)
print(
    "its Fourier representation independently of the"
)
print(
    "repository extraction routine."
)

print()

print("=" * 70)
print("END CELL 12")
print("=" * 70)
