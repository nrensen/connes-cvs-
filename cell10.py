# ============================================================
# CELL 10 — STRUCTURAL SPECTRAL DICTIONARY AUDIT
#
# Goal:
#
# Resolve the Cell-9 spectral-function discrepancy symbolically.
#
# Cell 9 established, to ~70 digits, that
#
#     F(tau) = integral_0^L f(t) exp(+i*tau*t) dt
#
# agrees with its closed-form Fourier representation.
#
# The repository spectral representation uses
#
#     g_k(tau)
#       = (exp(-i*tau*L) - 1)
#         / (i*(2*pi*k/L - tau))
#
# and a centering phase.
#
# We derive the exact relationship between:
#
#     F(tau)
#     sum u_k g_k(tau)
#     exp(+/- i*tau*L/2) F(tau)
#
# and determine which centering convention is correct.
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

mp.mp.dps = 70

c = 13
N = 6
T = 40

L = compute_L(c)
omega = 2 * mp.pi / L

print("=" * 70)
print("CELL 10 — STRUCTURAL SPECTRAL DICTIONARY AUDIT")
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
# 1. BASIS
# ============================================================

def phi(k, t):
    if k == 0:
        return 1 / mp.sqrt(L)

    return (
        mp.sqrt(2 / L)
        * mp.cos(2 * mp.pi * k * t / L)
    )


# ============================================================
# 2. FULL COEFFICIENTS
# ============================================================

def full_coefficients(v):
    u = canonical_to_full(v, N)

    return {
        k: mp.mpf(u[k + N])
        for k in range(-N, N + 1)
    }


# ============================================================
# 3. DIRECT FUNCTION
# ============================================================

def f_direct(v, t):
    total = mp.mpf("0")

    for k in range(N + 1):
        total += v[k] * phi(k, t)

    return total


# ============================================================
# 4. DIRECT FOURIER TRANSFORM
#
#     F(tau) = integral_0^L f(t) exp(+i*tau*t) dt
# ============================================================

def F_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda t:
            f_direct(v, t)
            * mp.exp(1j * tau * t),
        [0, L],
    )


# ============================================================
# 5. CLOSED FORM FOURIER TRANSFORM
# ============================================================

def exponential_integral(a):
    a = mp.mpf(a)

    if a == 0:
        return L

    return (
        mp.exp(1j * a * L) - 1
    ) / (1j * a)


def F_closed(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = 0j

    for k in range(-N, N + 1):

        a = tau + 2 * mp.pi * k / L

        total += (
            u[k]
            * exponential_integral(a)
            / mp.sqrt(L)
        )

    return total


# ============================================================
# 6. REPOSITORY g_k
# ============================================================

def g_repo(k, tau):

    tau = mp.mpf(tau)

    denom = (
        2 * mp.pi * k / L
        - tau
    )

    if denom == 0:
        return L

    return (
        mp.exp(-1j * tau * L) - 1
    ) / (
        1j * denom
    )


# ============================================================
# 7. RAW REPOSITORY SUM
#
#     S(tau) = sum_k u_k g_k(tau)
#
# ============================================================

def repo_sum(v, tau):

    u = full_coefficients(v)

    total = 0j

    for k in range(-N, N + 1):
        total += (
            u[k]
            * g_repo(k, tau)
        )

    return total


# ============================================================
# 8. SYMBOLICALLY DERIVED RELATION
#
# Starting from
#
# F(tau)
#   = 1/sqrt(L) sum_k u_k
#       integral_0^L
#       exp(i*(tau + omega*k)t) dt
#
# and using exp(i*omega*k*L) = 1:
#
# F(tau)
#   = exp(i*tau*L)/sqrt(L)
#       sum_k u_k g_{-k}(tau)
#
# Since u_{-k} = u_k:
#
# F(tau)
#   = exp(i*tau*L)/sqrt(L)
#       sum_k u_k g_k(tau)
#
# Hence
#
#     exp(+i*tau*L/2) S/sqrt(L)
#       = exp(-i*tau*L/2) F
#
# where S = sum u_k g_k.
#
# ============================================================

def repo_derived(v, tau):
    """
    Quantity represented by the repository's formula,
    derived directly from F_closed.
    """
    F = F_closed(v, tau)

    return (
        mp.exp(-1j * tau * L / 2)
        * F
    )


def repo_formula(v, tau):
    """
    The formula actually used in Cell 9 / repository:
        Re[ exp(+i*tau*L/2) S ] / sqrt(L)
    """
    S = repo_sum(v, tau)

    return (
        mp.exp(1j * tau * L / 2)
        * S
        / mp.sqrt(L)
    )


def repo_formula_minus(v, tau):
    """
    Same raw sum, but with the opposite centering phase.
    """
    S = repo_sum(v, tau)

    return (
        mp.exp(-1j * tau * L / 2)
        * S
        / mp.sqrt(L)
    )


# ============================================================
# 9. GROUND STATE
# ============================================================

print("-" * 70)
print("1. GROUND STATE")
print("-" * 70)
print()

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=mp.mp.dps,
)

lambda_min, u_star = compute_ground_state(Q)

v_star = mp.matrix(N + 1, 1)

v_star[0] = u_star[N, 0]

for k in range(1, N + 1):
    v_star[k] = (
        mp.sqrt(2)
        * u_star[N + k, 0]
    )

v_norm = mp.sqrt(
    mp.fdot(v_star, v_star)
)

for k in range(N + 1):
    v_star[k] /= v_norm

print("lambda_min =")
print(mp.nstr(lambda_min, 50))
print()

print("||v_star|| =")
print(mp.nstr(
    mp.sqrt(mp.fdot(v_star, v_star)),
    30,
))
print()


# ============================================================
# 10. STRUCTURAL IDENTITY
#
# Compare the raw repository sum against F.
#
# Test:
#
#   F
#   =
#   exp(i*tau*L) S/sqrt(L)
#
# ============================================================

print("-" * 70)
print("2. RAW FOURIER / REPOSITORY-SUM IDENTITY")
print("-" * 70)
print()

tau_values = [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("3"),
    mp.mpf("5"),
    mp.mpf("10"),
]

max_raw_identity_error = mp.mpf("0")

for tau in tau_values:

    F = F_closed(v_star, tau)

    S = repo_sum(v_star, tau)

    rhs = (
        mp.exp(1j * tau * L)
        * S
        / mp.sqrt(L)
    )

    err = abs(F - rhs)

    max_raw_identity_error = max(
        max_raw_identity_error,
        err,
    )

    print(
        f"tau = {mp.nstr(tau, 15)}"
    )
    print(
        "  F(tau) =",
        mp.nstr(F, 40),
    )
    print(
        "  exp(+i*tau*L) S/sqrt(L) =",
        mp.nstr(rhs, 40),
    )
    print(
        "  error =",
        mp.nstr(err, 20),
    )
    print()


print(
    "Maximum raw identity error =",
    mp.nstr(
        max_raw_identity_error,
        30,
    ),
)
print()


# ============================================================
# 11. CENTERING SIGN AUDIT
#
# The two competing possibilities are:
#
#   A = Re[ exp(+i*tau*L/2) F ] / sqrt(L)
#
#   B = Re[ exp(-i*tau*L/2) F ] / sqrt(L)
#
# The structural derivation predicts B for the centred
# transform corresponding to t -> t + L/2.
# ============================================================

print("-" * 70)
print("3. CENTERING-SIGN AUDIT")
print("-" * 70)
print()

max_plus_error = mp.mpf("0")
max_minus_error = mp.mpf("0")

for tau in tau_values:

    F = F_closed(v_star, tau)

    plus_direct = (
        mp.re(
            mp.exp(
                1j * tau * L / 2
            ) * F
        )
        / mp.sqrt(L)
    )

    minus_direct = (
        mp.re(
            mp.exp(
                -1j * tau * L / 2
            ) * F
        )
        / mp.sqrt(L)
    )

    repo = mp.re(
        repo_formula(v_star, tau)
    )

    repo_minus = mp.re(
        repo_formula_minus(
            v_star,
            tau,
        )
    )

    err_plus = abs(
        plus_direct - repo
    )

    err_minus = abs(
        minus_direct - repo_minus
    )

    max_plus_error = max(
        max_plus_error,
        err_plus,
    )

    max_minus_error = max(
        max_minus_error,
        err_minus,
    )

    print(
        f"tau = {mp.nstr(tau, 15)}"
    )
    print(
        "  + phase direct =",
        mp.nstr(plus_direct, 40),
    )
    print(
        "  + phase repo   =",
        mp.nstr(repo, 40),
    )
    print(
        "  + phase error  =",
        mp.nstr(err_plus, 20),
    )
    print()
    print(
        "  - phase direct =",
        mp.nstr(minus_direct, 40),
    )
    print(
        "  - phase repo   =",
        mp.nstr(repo_minus, 40),
    )
    print(
        "  - phase error  =",
        mp.nstr(err_minus, 20),
    )
    print()


print(
    "Maximum +phase error =",
    mp.nstr(
        max_plus_error,
        30,
    ),
)

print(
    "Maximum -phase error =",
    mp.nstr(
        max_minus_error,
        30,
    ),
)

print()


# ============================================================
# 12. DIRECT CENTERED-INTERVAL CHECK
#
# Independently define
#
#   f_c(x) = f(x + L/2)
#
# on -L/2 <= x <= L/2.
#
# Its Fourier transform is
#
#   integral_{-L/2}^{L/2}
#       f(x+L/2) exp(i*tau*x) dx
#
# = exp(-i*tau*L/2) F(tau).
#
# This is the decisive independent sign check.
# ============================================================

def F_centered_direct(v, tau):

    tau = mp.mpf(tau)

    return mp.quad(
        lambda x:
            f_direct(
                v,
                x + L / 2,
            )
            * mp.exp(1j * tau * x),
        [-L / 2, L / 2],
    )


print("-" * 70)
print("4. DIRECT CENTERED-INTERVAL CHECK")
print("-" * 70)
print()

max_centering_identity_error = mp.mpf("0")

for tau in tau_values:

    F = F_closed(v_star, tau)

    Fc = F_centered_direct(
        v_star,
        tau,
    )

    predicted = (
        mp.exp(
            -1j * tau * L / 2
        ) * F
    )

    err = abs(
        Fc - predicted
    )

    max_centering_identity_error = max(
        max_centering_identity_error,
        err,
    )

    print(
        f"tau = {mp.nstr(tau, 15)}"
    )
    print(
        "  direct centred transform =",
        mp.nstr(Fc, 40),
    )
    print(
        "  exp(-i*tau*L/2) F(tau)   =",
        mp.nstr(predicted, 40),
    )
    print(
        "  error =",
        mp.nstr(err, 20),
    )
    print()


print(
    "Maximum direct-centering identity error =",
    mp.nstr(
        max_centering_identity_error,
        30,
    ),
)
print()


# ============================================================
# 13. FINAL DIAGNOSTIC
# ============================================================

print("=" * 70)
print("CELL 10 SUMMARY")
print("=" * 70)
print()

print(
    "Maximum raw Fourier/repository-sum identity error =",
    mp.nstr(
        max_raw_identity_error,
        30,
    ),
)

print(
    "Maximum +phase error =",
    mp.nstr(
        max_plus_error,
        30,
    ),
)

print(
    "Maximum -phase error =",
    mp.nstr(
        max_minus_error,
        30,
    ),
)

print(
    "Maximum direct-centering identity error =",
    mp.nstr(
        max_centering_identity_error,
        30,
    ),
)

print()
print(
    "Interpretation:"
)
print()
print(
    "The centred transform of f(x+L/2) should be"
)
print(
    "    exp(-i*tau*L/2) F(tau)."
)
print()
print(
    "If the -phase errors are near numerical precision"
)
print(
    "while the +phase errors are large, the Cell-9"
)
print(
    "discrepancy is a centering-phase sign error."
)
print()
print("=" * 70)
print("END CELL 10")
print("=" * 70)
