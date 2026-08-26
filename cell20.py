# ============================================================
# cell20.py — CORRECTED CELL-5 ARCHIMEDEAN AUDIT
#
# Purpose
# -------
# This cell is derived directly from Cell 5.
#
# The sole mathematical correction is to replace the historical
# Cell-5 explicit Archimedean quantity
#
#     (1/pi) int_0^T h_+(r) Re(G_complex(v,r)) dr
#
# with the mathematically required quadratic functional
#
#     (1/pi) int_0^T h_+(r)
#         int_0^L K_v(1-y/L) cos(r*y) dy dr.
#
# This is the same correction implemented literally in the
# abandoned cell5_corrected.py, but uses the closed finite-
# Fourier representation K_fourier() established and validated
# in Cell 17, avoiding the nested numerical convolution.
#
# Everything else is retained from Cell 5 as far as practical.
#
# In particular:
#
#   - same c, N, T
#   - same ground-state construction
#   - same canonical/full-space conversion
#   - same prime matrix
#   - same pole matrix
#   - same repository Galerkin matrix
#
# The purpose is to determine whether the corrected explicit
# Archimedean calculation now agrees with the repository
# Archimedean quadratic form.
# ============================================================

import time

import mpmath as mp

from connes_cvs import (
    build_galerkin_matrix,
    compute_ground_state,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 80

c = 13
N = 8
T = 60
dps = 80

L = mp.log(c)
beta = L / (4 * mp.pi)

DISPLAY_DIGITS = 50


def nstr(x):
    return mp.nstr(x, DISPLAY_DIGITS)


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# 1. BUILD GROUND STATE
# ============================================================

print("=" * 72)
print("CELL 20 — CORRECTED CELL-5 ARCHIMEDEAN AUDIT")
print("=" * 72)

print()
print("Parameters:")
print("  c   =", c)
print("  N   =", N)
print("  T   =", T)
print("  dps =", dps)
print("  L   =", nstr(L))

print()
print("-" * 72)
print("1. GROUND STATE")
print("-" * 72)

ground_start = time.perf_counter()

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=dps,
)

lam_min, eigvec = compute_ground_state(Q)

ground_elapsed = elapsed(ground_start)

print()
print("lambda_min =")
print(nstr(lam_min))

print()
print(
    f"ground-state calculation elapsed = "
    f"{ground_elapsed:.6f} s"
)


# ============================================================
# 2. CANONICAL GROUND-STATE VECTOR
# ============================================================

coefficients = [
    mp.mpf(eigvec[i, 0])
    for i in range(eigvec.rows)
]

norm = mp.sqrt(
    sum(
        x * x
        for x in coefficients
    )
)

coefficients = [
    x / norm
    for x in coefficients
]


# Canonical real-even coefficients:
#
#   v_0 = c_0
#   v_k = sqrt(2) c_k, k >= 1
#

v_star = mp.matrix(N + 1, 1)

v_star[0] = coefficients[N]

for k in range(1, N + 1):
    v_star[k] = (
        mp.sqrt(2)
        * coefficients[N + k]
    )


print()
print("Canonical ground-state norm =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                v_star,
                v_star,
            )
        )
    )
)


# ============================================================
# 3. FULL-SYMMETRIC COEFFICIENT VECTOR
# ============================================================
#
# The full vector is
#
#     u_{0} = v_0
#
#     u_{+k} = u_{-k} = v_k / sqrt(2)
#
# for k > 0.
#
# This is the representation used by the repository quadratic
# forms and by Cell 17's K_fourier().
# ============================================================

u_star = mp.matrix(
    2 * N + 1,
    1,
)

for m in range(-N, N + 1):

    if m == 0:

        u_star[m + N] = v_star[0]

    else:

        u_star[m + N] = (
            v_star[abs(m)]
            / mp.sqrt(2)
        )


print()
print("Full-space dimension:")
print(
    "  u_star =",
    u_star.rows,
    "x",
    u_star.cols,
)

print()
print("||u_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                u_star,
                u_star,
            )
        )
    )
)


# ============================================================
# 4. PRIME-POWER LIST
# ============================================================

def prime_power_terms_cell20(c):
    """
    Return (q, Lambda(q)) for all prime powers q <= c.
    """

    c_int = int(
        mp.floor(c)
    )

    terms = []

    for p0 in range(
        2,
        c_int + 1,
    ):

        is_prime = True

        for d0 in range(
            2,
            int(mp.sqrt(p0)) + 1,
        ):

            if p0 % d0 == 0:

                is_prime = False
                break

        if not is_prime:
            continue

        q0 = p0

        while q0 <= c_int:

            terms.append(
                (
                    mp.mpf(q0),
                    mp.log(p0),
                )
            )

            q0 *= p0

    return terms


prime_terms = (
    prime_power_terms_cell20(c)
)


# ============================================================
# 5. PRIME MATRIX
# ============================================================

def Q_prime_power_cell20(
    q,
    Lambda_q,
):

    size = 2 * N + 1

    Qq = mp.matrix(
        size,
        size,
    )

    a = (
        1
        - mp.log(q) / L
    )

    prefactor = (
        -1 / mp.pi
        * Lambda_q
        / mp.sqrt(q)
    )

    def psi(x):

        return (
            prefactor
            * mp.sin(
                2 * mp.pi
                * x
                * a
            )
        )

    def psi_derivative(x):

        return (
            prefactor
            * 2
            * mp.pi
            * a
            * mp.cos(
                2 * mp.pi
                * x
                * a
            )
        )

    for i, m in enumerate(
        range(-N, N + 1)
    ):

        for j, n in enumerate(
            range(-N, N + 1)
        ):

            if m != n:

                Qq[i, j] = (
                    psi(m)
                    - psi(n)
                ) / mp.mpf(
                    m - n
                )

            else:

                Qq[i, j] = (
                    psi_derivative(m)
                )

    return Qq


Q_prime = mp.matrix(
    2 * N + 1,
    2 * N + 1,
)

for q, Lambda_q in prime_terms:

    Q_prime += (
        Q_prime_power_cell20(
            q,
            Lambda_q,
        )
    )


prime_matrix_form = mp.fdot(
    u_star,
    Q_prime * u_star,
)


print()
print("-" * 72)
print("2. PRIME MATRIX")
print("-" * 72)

print()
print("<u_star, Q_prime u_star> =")
print(
    nstr(prime_matrix_form)
)


# ============================================================
# 6. POLE MATRIX
# ============================================================

def psi_pole_cell20(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * mp.sin(
            2 * mp.pi
            * x
            * (
                1
                - y / L
            )
        )
    )

    return (
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L],
        )
    )


def psi_pole_derivative_cell20(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * (
            2 * mp.pi
            * (
                1
                - y / L
            )
        )
        * mp.cos(
            2 * mp.pi
            * x
            * (
                1
                - y / L
            )
        )
    )

    return (
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L],
        )
    )


size = 2 * N + 1

Q_pole = mp.matrix(
    size,
    size,
)

for i, m in enumerate(
    range(-N, N + 1)
):

    for j, n in enumerate(
        range(-N, N + 1)
    ):

        if m != n:

            Q_pole[i, j] = (
                psi_pole_cell20(m)
                - psi_pole_cell20(n)
            ) / mp.mpf(
                m - n
            )

        else:

            Q_pole[i, j] = (
                psi_pole_derivative_cell20(m)
            )


pole_matrix_form = mp.fdot(
    u_star,
    Q_pole * u_star,
)


print()
print("-" * 72)
print("3. POLE MATRIX")
print("-" * 72)

print()
print("<u_star, Q_pole u_star> =")
print(
    nstr(pole_matrix_form)
)


# ============================================================
# 7. POLE SANITY CHECK
# ============================================================
#
# Retain the historical Cell-5 check:
#
#     <u,Q_pole u> = Re(2 G_complex(v,i/2))
#
# G_complex is replaced by sum_v_G from the current cell.py
# interface, which takes L explicitly.
# ============================================================

from cell import sum_v_G


pole_explicit = mp.re(
    2 * sum_v_G(
        v_star,
        1j / 2,
        L,
    )
)

print()
print("-" * 72)
print("4. POLE SANITY CHECK")
print("-" * 72)

print()
print("<u_star, Q_pole u_star> =")
print(
    nstr(pole_matrix_form)
)

print()
print("2 Re sum_v_G(v_star, i/2) =")
print(
    nstr(pole_explicit)
)

print()
print("difference =")
print(
    nstr(
        pole_matrix_form
        - pole_explicit
    )
)


# ============================================================
# 8. h_+(r)
# ============================================================

def h_plus_cell20(r):

    r = mp.mpf(r)

    return (
        mp.re(
            mp.digamma(
                mp.mpf("0.25")
                + 1j * r / 2
            )
        )
        - mp.log(mp.pi)
    )


# ============================================================
# 9. CLOSED-FORM QUADRATIC KERNEL
# ============================================================
#
# This is the efficient version of the K_v calculation
# established in Cell 17.
#
# K_v(omega)
#
#   = 2 int_0^omega
#       T_v(t) T_v(omega-t) dt
#
# but evaluated directly from the finite Fourier expansion.
#
# Cell 17 independently verified this against K_canonical().
# ============================================================

def K_fourier(v, omega):

    omega = mp.mpf(omega)

    if omega <= 0:
        return mp.mpf("0")

    if omega >= 1:
        raise ValueError(
            "K_fourier expects "
            "0 < omega < 1"
        )

    # Construct full symmetric Fourier coefficients.
    #
    # u_0 = v_0
    #
    # u_{+k} = u_{-k} = v_k / sqrt(2)

    u = mp.matrix(
        2 * N + 1,
        1,
    )

    for m in range(
        -N,
        N + 1,
    ):

        if m == 0:

            u[m + N] = v[0]

        else:

            u[m + N] = (
                v[abs(m)]
                / mp.sqrt(2)
            )

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
                        2j
                        * mp.pi
                        * delta
                        * omega
                    )
                    / (
                        2j
                        * mp.pi
                        * delta
                    )
                )

            total += (
                2
                * um
                * un
                * mp.exp(
                    2j
                    * mp.pi
                    * n
                    * omega
                )
                * integral
            )

    return mp.re(total)


# ============================================================
# 10. DIRECT K-BASED ARCHIMEDEAN FUNCTIONAL
# ============================================================
#
# This is the efficient implementation of the exact correction
# introduced by cell5_corrected.py.
#
# It is algebraically equivalent to:
#
#   (1/pi^2) int h_+(r) D_v(r) dr
#
# because Cell 17 established
#
#   D_v = pi K_v.
#
# ============================================================

def K_y(v, y):

    omega = (
        1
        - y / L
    )

    if omega <= 0:

        return mp.mpf("0")

    return K_fourier(
        v,
        omega,
    )


def corrected_arch_integrand(
    v,
    r,
):

    r = mp.mpf(r)

    inner = mp.quad(
        lambda y:
            K_y(v, y)
            * mp.cos(r * y),
        [0, L],
    )

    return (
        h_plus_cell20(r)
        * inner
    )


def corrected_arch_explicit(
    v,
    T,
):

    T = mp.mpf(T)

    return (
        mp.quad(
            lambda r:
                corrected_arch_integrand(
                    v,
                    r,
                ),
            [0, T],
        )
        / mp.pi
    )


# ============================================================
# 11. CORRECTED EXPLICIT ARCHIMEDEAN VALUE
# ============================================================

print()
print("-" * 72)
print("5. CORRECTED EXPLICIT ARCHIMEDEAN FUNCTIONAL")
print("-" * 72)

print()
print(
    "Computing:"
)

print(
    "  A_arch = (1/pi) int_0^T h_+(r)"
)

print(
    "          * int_0^L"
)

print(
    "              K_v(1-y/L)"
    " cos(r*y) dy dr"
)

print()

arch_start = time.perf_counter()

explicit_arch = (
    corrected_arch_explicit(
        v_star,
        T,
    )
)

arch_elapsed = elapsed(
    arch_start
)

print()
print(
    "Corrected explicit Archimedean ="
)

print(
    nstr(explicit_arch)
)

print()
print(
    f"corrected Archimedean elapsed = "
    f"{arch_elapsed:.6f} s"
)


# ============================================================
# 12. REPOSITORY ARCHIMEDEAN MATRIX
# ============================================================
#
# The repository matrix contains the complete finite-T
# quadratic form:
#
#     Q_T
#
# Remove the prime and pole contributions to isolate:
#
#     Q_arch = Q_T - Q_prime - Q_pole
# ============================================================

print()
print("-" * 72)
print("6. REPOSITORY ARCHIMEDEAN MATRIX")
print("-" * 72)

repo_start = time.perf_counter()

Q_arch = (
    Q
    - Q_prime
    - Q_pole
)

repo_arch = mp.fdot(
    u_star,
    Q_arch * u_star,
)

repo_arch_elapsed = elapsed(
    repo_start
)

print()
print(
    "Repository Archimedean ="
)

print(
    nstr(repo_arch)
)

print()
print(
    f"matrix subtraction/form elapsed = "
    f"{repo_arch_elapsed:.6f} s"
)


# ============================================================
# 13. COMPARISON
# ============================================================

difference = (
    repo_arch
    - explicit_arch
)

relative_difference = (
    abs(difference)
    / max(
        abs(repo_arch),
        abs(explicit_arch),
        mp.mpf("1"),
    )
)

print()
print("-" * 72)
print("7. CORRECTED EXPLICIT VS REPOSITORY")
print("-" * 72)

print()
print("Repository Archimedean =")
print(
    nstr(repo_arch)
)

print()
print("Corrected explicit =")
print(
    nstr(explicit_arch)
)

print()
print("difference =")
print(
    nstr(difference)
)

print()
print("absolute difference =")
print(
    nstr(abs(difference))
)

print()
print("relative difference =")
print(
    nstr(relative_difference)
)


# ============================================================
# 14. HISTORICAL LINEAR QUANTITY — CONTRAST ONLY
# ============================================================
#
# This is deliberately retained as a labelled contrast.
# It is NOT treated as a candidate Archimedean value.
# ============================================================

print()
print("-" * 72)
print("8. HISTORICAL CELL-5 LINEAR QUANTITY")
print("-" * 72)

linear_start = time.perf_counter()

historical_linear = (
    mp.quad(
        lambda r:
            h_plus_cell20(r)
            * mp.re(
                sum_v_G(
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

print()
print("Historical linear quantity =")
print(
    nstr(historical_linear)
)

print()
print(
    f"historical linear integral elapsed = "
    f"{linear_elapsed:.6f} s"
)

print()
print(
    "This quantity is retained solely for historical"
)

print(
    "comparison. It is NOT the corrected"
)

print(
    "Archimedean functional."
)


# ============================================================
# 15. FINAL DIAGNOSTIC
# ============================================================

print()
print("=" * 72)
print("CELL 20 — FINAL DIAGNOSTIC")
print("=" * 72)

print()
print("The corrected explicit Archimedean calculation uses")
print("the quadratic K_v functional identified in Cell 17.")
print()

print(
    "Repository Archimedean =",
    nstr(repo_arch),
)

print(
    "Corrected explicit      =",
    nstr(explicit_arch),
)

print(
    "Absolute difference     =",
    nstr(abs(difference)),
)

print(
    "Relative difference     =",
    nstr(relative_difference),
)

print()
print("Timing:")
print(
    f"  ground state          = "
    f"{ground_elapsed:.6f} s"
)

print(
    f"  corrected Archimedean = "
    f"{arch_elapsed:.6f} s"
)

print(
    f"  historical linear     = "
    f"{linear_elapsed:.6f} s"
)

print()
print("=" * 72)
print("CELL 20 COMPLETE")
print("=" * 72)
