# ============================================================
# CELL 17 — QUADRATIC ARCHIMEDEAN KERNEL AUDIT
#
# Purpose
# -------
# Establish the exact relationship between:
#
#   1. The Archimedean quadratic form obtained from
#      psi_arch -> divided differences.
#
#   2. The quadratic Volterra kernel
#
#         K_v(omega) =
#             2 int_0^omega T_v(t) T_v(omega-t) dt
#
#   3. The direct Archimedean integral built from K_v.
#
# The key algebraic identity under test is
#
#   D_v(omega) = pi K_v(omega),
#
# where
#
#   D_v(omega)
#       = sum_{m,n} u_m u_n
#           [sin(2*pi*m*omega) - sin(2*pi*n*omega)]
#           / (m-n)
#
# with the diagonal interpreted by continuity.
#
# The Archimedean functional used by the repository is
#
#   psi_arch(x)
#       = 1/(2*pi^2)
#         int_{-T}^T h_+(r) Re(S_hat_x(r)) dr.
#
# Consequently, after taking the quadratic divided difference
# and using D_v = pi K_v,
#
#   A_arch(v)
#       = 1/pi
#         int_0^T h_+(r)
#           int_0^L
#             K_v(1-y/L) cos(r*y) dy dr.
#
# This cell deliberately does NOT use G_complex().
#
# It is intended to settle the Cell-5 discrepancy by comparing
# two genuinely quadratic constructions.
# ============================================================

import time
import mpmath as mp

from cell import (
    compute_L,
    canonical_to_full,
    full_to_canonical,
    T_canonical,
    K_canonical,
    get_ground_state,
)

from connes_cvs.operator import (
    psi_arch,
    psi_arch_deriv,
    h_plus,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 60

c = 13
N = 8
T = 40

DISPLAY_DIGITS = 40

TEST_OMEGA = [
    mp.mpf("0.05"),
    mp.mpf("0.125"),
    mp.mpf("0.25"),
    mp.mpf("0.5"),
    mp.mpf("0.75"),
    mp.mpf("0.875"),
    mp.mpf("0.95"),
]

TEST_K_NUMERICAL = [
    mp.mpf("0.125"),
    mp.mpf("0.25"),
    mp.mpf("0.5"),
    mp.mpf("0.75"),
]


# ============================================================
# HELPERS
# ============================================================

def nstr(x):
    return mp.nstr(x, DISPLAY_DIGITS)


def elapsed(start):
    return time.perf_counter() - start


def full_coefficients(v):
    """
    Return the full symmetric Fourier coefficient vector u
    corresponding to canonical real-even coefficients v.
    """
    return canonical_to_full(
        v,
        len(v) - 1,
    )


# ============================================================
# 1. HEADER
# ============================================================

print("=" * 78)
print("CELL 17 — QUADRATIC ARCHIMEDEAN KERNEL AUDIT")
print("=" * 78)
print()

print("Parameters:")
print(f"  c   = {c}")
print(f"  N   = {N}")
print(f"  T   = {T}")
print(f"  dps = {mp.mp.dps}")
print(
    f"  L   = {nstr(compute_L(c))}"
)
print()


# ============================================================
# 2. GROUND STATE
# ============================================================

print("-" * 78)
print("1. GROUND STATE")
print("-" * 78)
print()

L = compute_L(c)

ground_start = time.perf_counter()

lambda_min, u_star, ground_meta = get_ground_state(
    c=c,
    N=N,
    T=T,
    dps=mp.mp.dps,
    verbose=True,
)

ground_elapsed = elapsed(ground_start)

u_star = mp.matrix(u_star)

v_star = full_to_canonical(
    u_star,
    N,
)

print()
print("Ground-state result:")
print(f"  lambda_min = {nstr(lambda_min)}")
print(
    f"  ||u_star|| = "
    f"{nstr(mp.sqrt(mp.fdot(u_star, u_star)))}"
)
print(
    f"  ||v_star|| = "
    f"{nstr(mp.sqrt(mp.fdot(v_star, v_star)))}"
)
print(
    f"  elapsed    = {ground_elapsed:.6f} s"
)
print()


# ============================================================
# 3. FULL-SYMMETRIC COEFFICIENTS
# ============================================================

print("-" * 78)
print("2. COORDINATE CONSISTENCY")
print("-" * 78)
print()

u_from_v = canonical_to_full(
    v_star,
    N,
)

roundtrip_error = mp.sqrt(
    mp.fdot(
        u_from_v - u_star,
        u_from_v - u_star,
    )
)

print(
    "canonical -> full -> canonical consistency:"
)
print(
    f"  ||u(v_star) - u_star|| = "
    f"{nstr(roundtrip_error)}"
)
print()


# ============================================================
# 4. CLOSED-FORM QUADRATIC KERNEL
# ============================================================

def K_fourier(v, omega):
    """
    Closed finite-Fourier evaluation of

        K_v(omega)
          = 2 int_0^omega
              T_v(t) T_v(omega-t) dt.

    This is algebraically equivalent to K_canonical(), but avoids
    a nested numerical quadrature.
    """
    omega = mp.mpf(omega)

    if omega <= 0:
        return mp.mpf("0")

    if omega >= 1:
        raise ValueError(
            "K_fourier expects 0 < omega < 1"
        )

    u = full_coefficients(v)

    total = mp.mpc("0")

    for i, m in enumerate(
        range(-N, N + 1)
    ):
        um = u[i]

        for j, n in enumerate(
            range(-N, N + 1)
        ):
            un = u[j]

            if m == n:
                integral = omega
            else:
                delta = m - n

                integral = (
                    mp.expm1(
                        2j * mp.pi
                        * delta
                        * omega
                    )
                    /
                    (
                        2j * mp.pi * delta
                    )
                )

            total += (
                2
                * um
                * un
                * mp.exp(
                    2j * mp.pi * n * omega
                )
                * integral
            )

    return mp.re(total)


# ============================================================
# 5. DIRECT ALGEBRAIC IDENTITY
# ============================================================

def D_divdiff(v, omega):
    """
    Direct quadratic divided-difference kernel D_v(omega).
    """
    omega = mp.mpf(omega)
    u = full_coefficients(v)

    total = mp.mpf("0")

    for i, m in enumerate(
        range(-N, N + 1)
    ):
        um = u[i]

        for j, n in enumerate(
            range(-N, N + 1)
        ):
            un = u[j]

            if m == n:
                kernel = (
                    2
                    * mp.pi
                    * omega
                    * mp.cos(
                        2 * mp.pi * n * omega
                    )
                )
            else:
                kernel = (
                    mp.sin(
                        2 * mp.pi * m * omega
                    )
                    -
                    mp.sin(
                        2 * mp.pi * n * omega
                    )
                ) / (m - n)

            total += (
                um
                * un
                * kernel
            )

    return total


print("-" * 78)
print("3. ALGEBRAIC IDENTITY: D_v(omega) = pi K_v(omega)")
print("-" * 78)
print()

max_identity_error = mp.mpf("0")

for omega in TEST_OMEGA:

    D = D_divdiff(
        v_star,
        omega,
    )

    K = K_fourier(
        v_star,
        omega,
    )

    error = abs(
        D - mp.pi * K
    )

    max_identity_error = max(
        max_identity_error,
        error,
    )

    print(
        f"omega = {nstr(omega)}"
    )
    print(
        f"  D_v       = {nstr(D)}"
    )
    print(
        f"  pi*K_v    = {nstr(mp.pi * K)}"
    )
    print(
        f"  abs error = {nstr(error)}"
    )
    print()

print(
    "Maximum |D_v - pi*K_v| = "
    f"{nstr(max_identity_error)}"
)
print()


# ============================================================
# 6. DIRECT NUMERICAL CONVOLUTION CHECK
# ============================================================

print("-" * 78)
print("4. NUMERICAL CONVOLUTION CROSS-CHECK")
print("-" * 78)
print()

max_K_error = mp.mpf("0")

for omega in TEST_K_NUMERICAL:

    start = time.perf_counter()

    K_num = K_canonical(
        v_star,
        omega,
    )

    K_closed = K_fourier(
        v_star,
        omega,
    )

    t = elapsed(start)

    error = abs(
        K_num - K_closed
    )

    max_K_error = max(
        max_K_error,
        error,
    )

    print(
        f"omega = {nstr(omega)}"
    )
    print(
        f"  K numerical = {nstr(K_num)}"
    )
    print(
        f"  K closed    = {nstr(K_closed)}"
    )
    print(
        f"  abs error   = {nstr(error)}"
    )
    print(
        f"  elapsed     = {t:.6f} s"
    )
    print()

print(
    "Maximum numerical/closed K error = "
    f"{nstr(max_K_error)}"
)
print()


# ============================================================
# 7. ARCHIMEDEAN DIVIDED-DIFFERENCE QUADRATIC FORM
# ============================================================

print("-" * 78)
print("5. ARCHIMEDEAN DIVIDED-DIFFERENCE FORM")
print("-" * 78)
print()

psi_arch_values = {}
psi_arch_derivatives = {}

psi_start = time.perf_counter()

for n in range(
    0,
    N + 1,
):

    psi_arch_values[n] = psi_arch(
        mp.mpf(n),
        L,
        T,
        mp.mp.dps,
    )

    psi_arch_derivatives[n] = (
        psi_arch_deriv(
            mp.mpf(n),
            L,
            T,
            mp.mp.dps,
        )
    )

    if n:
        psi_arch_values[-n] = (
            -psi_arch_values[n]
        )

        psi_arch_derivatives[-n] = (
            psi_arch_derivatives[n]
        )

psi_elapsed = elapsed(
    psi_start
)

print(
    f"psi_arch construction elapsed = "
    f"{psi_elapsed:.6f} s"
)
print()

Q_arch = mp.matrix(
    2 * N + 1,
    2 * N + 1,
)

for i in range(
    2 * N + 1
):
    m = i - N

    for j in range(
        i,
        2 * N + 1
    ):
        n = j - N

        if m == n:
            value = (
                psi_arch_derivatives[n]
            )
        else:
            value = (
                psi_arch_values[m]
                - psi_arch_values[n]
            ) / (m - n)

        Q_arch[i, j] = value

        if i != j:
            Q_arch[j, i] = value


A_divdiff = (
    u_star.T
    * Q_arch
    * u_star
)[0]

print(
    "A_arch via psi_arch divided differences ="
)
print(
    nstr(A_divdiff)
)
print()


# ============================================================
# 8. DIRECT K-BASED ARCHIMEDEAN INTEGRAL
# ============================================================

def K_y(v, y):
    """
    K_v(1-y/L), for 0 <= y <= L.
    """
    omega = 1 - y / L

    if omega <= 0:
        return mp.mpf("0")

    return K_fourier(
        v,
        omega,
    )


def direct_K_inner(r):
    """
    J(r) =
        int_0^L
          K_v(1-y/L) cos(r*y) dy.
    """
    r = mp.mpf(r)

    return mp.quad(
        lambda y:
            K_y(v_star, y)
            * mp.cos(r * y),
        [0, L],
    )


def direct_K_integrand(r):
    return (
        h_plus(
            mp.mpf(r),
            mp.mp.dps,
        )
        * direct_K_inner(r)
    )


print("-" * 78)
print("6. DIRECT K-BASED ARCHIMEDEAN INTEGRAL")
print("-" * 78)
print()

print(
    "Computing"
)
print(
    "  A_K = (1/pi) int_0^T h(r)"
)
print(
    "        * int_0^L K_v(1-y/L) cos(r*y) dy dr"
)
print()
print(
    "This is the intentionally expensive part of cell17."
)
print()

direct_start = time.perf_counter()

A_K = (
    mp.quad(
        direct_K_integrand,
        [0, T],
    )
    / mp.pi
)

direct_elapsed = elapsed(
    direct_start
)

print(
    "A_arch via direct K integral ="
)
print(
    nstr(A_K)
)
print()

print(
    f"direct K integral elapsed = "
    f"{direct_elapsed:.6f} s"
)
print()


# ============================================================
# 9. THREE-WAY COMPARISON
# ============================================================

print("-" * 78)
print("7. THREE-WAY COMPARISON")
print("-" * 78)
print()

print("A_divdiff =")
print(nstr(A_divdiff))
print()

print("A_K       =")
print(nstr(A_K))
print()

difference = abs(
    A_divdiff - A_K
)

relative_difference = (
    difference
    /
    max(
        abs(A_divdiff),
        abs(A_K),
        mp.mpf("1"),
    )
)

print("|A_divdiff - A_K| =")
print(nstr(difference))
print()

print("relative difference =")
print(nstr(relative_difference))
print()


# ============================================================
# 10. HISTORICAL CELL-5 LINEAR FUNCTIONAL
# ============================================================

from cell import G_complex

print("-" * 78)
print("8. HISTORICAL CELL-5 LINEAR QUANTITY (CONTRAST ONLY)")
print("-" * 78)
print()

linear_start = time.perf_counter()

A_linear = (
    mp.quad(
        lambda r:
            h_plus(
                r,
                mp.mp.dps,
            )
            * mp.re(
                G_complex(
                    v_star,
                    r,
                    L,
                )
            ),
        [0, T],
    )
    / mp.pi
)

linear_elapsed = elapsed(
    linear_start
)

print("A_linear =")
print(nstr(A_linear))
print()

print(
    f"linear integral elapsed = "
    f"{linear_elapsed:.6f} s"
)
print()

print(
    "This value is deliberately NOT treated as a candidate"
)
print(
    "for the quadratic Archimedean functional."
)
print()


# ============================================================
# 11. FINAL DIAGNOSTIC
# ============================================================

print("-" * 78)
print("9. FINAL DIAGNOSTIC")
print("-" * 78)
print()

print("Kernel identity:")
print(
    f"  max |D - pi*K| = "
    f"{nstr(max_identity_error)}"
)
print()

print("Numerical convolution:")
print(
    f"  max |K_numeric - K_closed| = "
    f"{nstr(max_K_error)}"
)
print()

print("Archimedean quadratic form:")
print(
    f"  |A_divdiff - A_K| = "
    f"{nstr(difference)}"
)
print(
    f"  relative difference = "
    f"{nstr(relative_difference)}"
)
print()

print("Historical linear quantity:")
print(
    f"  A_linear = {nstr(A_linear)}"
)
print()

print("=" * 78)
print("CELL 17 COMPLETE")
print("=" * 78)
print()
