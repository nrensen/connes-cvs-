# ============================================================
# CELL 21 — MODERN CELL-5 REIMPLEMENTATION
#
# Exploratory run:
#
#     working dps     = 40
#     forensic ground = (c=13, N=8, T=400, dps=150)
#
# The forensic ground state is retrieved from the persistent
# cache.  Cell 21 itself uses T=60 for the Galerkin calculation
# and for the explicit Archimedean integral.
#
# The principal Archimedean calculation uses the genuinely
# quadratic K_v construction.  It does NOT use sum_v_G.
# ============================================================

import time

import mpmath as mp

from connes_cvs import build_galerkin_matrix

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    canonical_to_full,
    full_to_canonical,
    compute_L,
    prime_power_terms,
    sum_v_G,
)


# ============================================================
# PARAMETERS
# ============================================================

WORKING_DPS = 20

# Cell-21 calculation cutoff.
T = 60

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

DISPLAY_DIGITS = 30


def nstr(x):
    return mp.nstr(x, DISPLAY_DIGITS)


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# HEADER
# ============================================================

print("=" * 72)
print("CELL 21 — MODERN CELL-5 REIMPLEMENTATION")
print("=" * 72)

print()
print("Parameters:")
print(f"  c              = {c}")
print(f"  N              = {N}")
print(f"  T              = {T}")
print(f"  working_dps    = {WORKING_DPS}")

print()
print("Forensic ground state:")
print(f"  c              = {FORENSIC_GROUND_STATE['c']}")
print(f"  N              = {FORENSIC_GROUND_STATE['N']}")
print(f"  T              = {FORENSIC_GROUND_STATE['T']}")
print(f"  generation_dps = {FORENSIC_GROUND_STATE['dps']}")

print()
print("L =")
print(nstr(L))


# ============================================================
# 1. FORENSIC GROUND STATE
#
# IMPORTANT:
#
# The ground state is deliberately obtained from the canonical
# FORENSIC_GROUND_STATE configuration.
#
# Its T=400 eigenvector is then evaluated against the T=60
# Galerkin matrix below.
# ============================================================

print()
print("-" * 72)
print("1. FORENSIC GROUND STATE")
print("-" * 72)

ground_start = time.perf_counter()

lambda_forensic, u_star, ground_meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
    verbose=True,
)

ground_elapsed = elapsed(ground_start)

u_star = mp.matrix(u_star)

print()
print("Forensic ground-state result:")
print("  lambda =")
print(nstr(lambda_forensic))

print()
print("  ||u_star|| =")
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

print()
print(
    f"ground-state retrieval elapsed = "
    f"{ground_elapsed:.6f} s"
)


# ============================================================
# 2. CANONICAL GROUND-STATE VECTOR
# ============================================================

print()
print("-" * 72)
print("2. CANONICAL GROUND-STATE VECTOR")
print("-" * 72)

v_star = full_to_canonical(
    u_star,
    N,
)

print()
print("||v_star|| =")
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
# 3. BUILD THE CELL-21 GALERKIN MATRIX
#
# This is Q_T with T=60.
#
# We need this because lambda_forensic belongs to the T=400
# ground-state problem and therefore cannot itself be used as
# the T=60 total quadratic form.
# ============================================================

print()
print("-" * 72)
print("3. CELL-21 GALERKIN MATRIX")
print("-" * 72)

matrix_start = time.perf_counter()

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=WORKING_DPS,
)

matrix_elapsed = elapsed(
    matrix_start
)

print()
print(
    f"Q(T={T}) construction elapsed = "
    f"{matrix_elapsed:.6f} s"
)


# ============================================================
# 4. TOTAL QUADRATIC FORM
#
#     total = u_star^T Q_T u_star
#
# This is the quadratic form evaluated on the forensic
# ground-state vector.
# ============================================================

print()
print("-" * 72)
print("4. TOTAL QUADRATIC FORM")
print("-" * 72)

total_start = time.perf_counter()

total_form = mp.fdot(
    u_star,
    Q * u_star,
)

total_elapsed = elapsed(
    total_start
)

print()
print("total =")
print(nstr(total_form))

print()
print(
    f"total quadratic-form elapsed = "
    f"{total_elapsed:.6f} s"
)


# ============================================================
# 5. PRIME MATRIX
#
# This follows the audited Cell-20 construction.
# ============================================================

print()
print("-" * 72)
print("5. PRIME MATRIX")
print("-" * 72)


def Q_prime_power(q, Lambda_q):

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

    values = {
        m: psi(m)
        for m in range(-N, N + 1)
    }

    derivatives = {
        m: psi_derivative(m)
        for m in range(-N, N + 1)
    }

    for i, m in enumerate(
        range(-N, N + 1)
    ):
        for j, n in enumerate(
            range(-N, N + 1)
        ):

            if m == n:
                Qq[i, j] = derivatives[m]

            else:
                Qq[i, j] = (
                    values[m]
                    - values[n]
                ) / mp.mpf(m - n)

    return Qq


prime_start = time.perf_counter()

Q_prime = mp.matrix(
    2 * N + 1,
    2 * N + 1,
)

for q, Lambda_q in prime_power_terms(c):
    Q_prime += Q_prime_power(
        q,
        Lambda_q,
    )

prime_form = mp.fdot(
    u_star,
    Q_prime * u_star,
)

prime_elapsed = elapsed(
    prime_start
)

print()
print("prime =")
print(nstr(prime_form))

print()
print(
    f"prime elapsed = "
    f"{prime_elapsed:.6f} s"
)


# ============================================================
# 6. POLE MATRIX
#
# This follows the audited Cell-20 construction.
#
# The basis responses are calculated once per basis index rather
# than once per matrix element.
# ============================================================

print()
print("-" * 72)
print("6. POLE MATRIX")
print("-" * 72)


def psi_pole(x):

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
        mp.quad(
            integrand,
            [0, L],
        )
        / mp.pi
    )


def psi_pole_derivative(x):

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
        mp.quad(
            integrand,
            [0, L],
        )
        / mp.pi
    )


pole_start = time.perf_counter()

pole_values = {
    m: psi_pole(m)
    for m in range(-N, N + 1)
}

pole_derivatives = {
    m: psi_pole_derivative(m)
    for m in range(-N, N + 1)
}

Q_pole = mp.matrix(
    2 * N + 1,
    2 * N + 1,
)

for i, m in enumerate(
    range(-N, N + 1)
):
    for j, n in enumerate(
        range(-N, N + 1)
    ):

        if m == n:
            Q_pole[i, j] = (
                pole_derivatives[m]
            )

        else:
            Q_pole[i, j] = (
                pole_values[m]
                - pole_values[n]
            ) / mp.mpf(m - n)

pole_form = mp.fdot(
    u_star,
    Q_pole * u_star,
)

pole_elapsed = elapsed(
    pole_start
)

print()
print("pole =")
print(nstr(pole_form))

print()
print(
    f"pole elapsed = "
    f"{pole_elapsed:.6f} s"
)


# ============================================================
# 7. ARCHIMEDEAN MATRIX FORM
#
#     Q_arch = Q_T - Q_prime - Q_pole
#
# and
#
#     arch = u_star^T Q_arch u_star
#
# Equivalently:
#
#     arch = total - prime - pole
#
# We use the scalar subtraction to avoid another matrix
# multiplication.
# ============================================================

print()
print("-" * 72)
print("7. ARCHIMEDEAN MATRIX FORM")
print("-" * 72)

arch_matrix = (
    total_form
    - prime_form
    - pole_form
)

print()
print("arch = total - prime - pole =")
print(nstr(arch_matrix))


# ============================================================
# 8. h_+(r)
# ============================================================

def h_plus(r):

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
#
# This is the Cell-17 K_fourier construction, reproduced here
# rather than importing Cell 17's executable machinery.
#
# K_v(omega)
#     = 2 int_0^omega
#           T_v(t) T_v(omega-t) dt
#
# evaluated analytically from the finite Fourier expansion.
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

    u = canonical_to_full(
        v,
        N,
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
#
#     A_arch =
#
#       (1/pi) int_0^T h_+(r)
#           int_0^L
#               K_v(1-y/L) cos(r*y) dy dr
#
# This is the corrected Cell-5 quantity.
#
# IMPORTANT:
#
#     sum_v_G is NOT used here.
# ============================================================

def K_y(y):

    omega = (
        1
        - y / L
    )

    if omega <= 0:
        return mp.mpf("0")

    return K_fourier(
        v_star,
        omega,
    )


def corrected_arch_integrand(r):

    r = mp.mpf(r)

    inner = mp.quad(
        lambda y:
            K_y(y)
            * mp.cos(r * y),
        [0, L],
    )

    return (
        h_plus(r)
        * inner
    )


print()
print("-" * 72)
print("8. EXPLICIT ARCHIMEDEAN")
print("-" * 72)

print()
print(
    "Computing corrected quadratic "
    "K-Fourier Archimedean functional..."
)

arch_start = time.perf_counter()

explicit_arch = (
    mp.quad(
        corrected_arch_integrand,
        [0, T],
    )
    / mp.pi
)

arch_elapsed = elapsed(
    arch_start
)

print()
print("K-fourier =")
print(nstr(explicit_arch))

print()
print(
    f"explicit Archimedean elapsed = "
    f"{arch_elapsed:.6f} s"
)


# ============================================================
# 11. AGREEMENT
# ============================================================

difference = (
    arch_matrix
    - explicit_arch
)

relative_difference = (
    abs(difference)
    / max(
        abs(arch_matrix),
        abs(explicit_arch),
        mp.mpf("1"),
    )
)

print()
print("-" * 72)
print("9. AGREEMENT")
print("-" * 72)

print()
print("matrix arch =")
print(nstr(arch_matrix))

print()
print("explicit arch =")
print(nstr(explicit_arch))

print()
print("difference =")
print(nstr(difference))

print()
print("relative =")
print(nstr(relative_difference))


# ============================================================
# 12. HISTORICAL CELL-5 LINEAR QUANTITY
#
# Deliberately retained ONLY as a forensic comparison.
#
#     (1/pi) int h_+(r)
#         Re(sum_v_G(v,r,L)) dr
#
# This is linear in v and is NOT the Archimedean quadratic
# functional.
# ============================================================

print()
print("-" * 72)
print("10. HISTORICAL CELL-5 LINEAR QUANTITY")
print("-" * 72)

linear_start = time.perf_counter()

historical_linear = (
    mp.quad(
        lambda r:
            h_plus(r)
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
print("linear result =")
print(nstr(historical_linear))

print()
print(
    "[WRONG; retained for forensic comparison]"
)


# ============================================================
# 13. TIMING / FINAL SUMMARY
# ============================================================

print()
print("=" * 72)
print("CELL 21 COMPLETE")
print("=" * 72)

print()
print("Quadratic forms:")
print(
    "  total          =",
    nstr(total_form),
)
print(
    "  prime          =",
    nstr(prime_form),
)
print(
    "  pole           =",
    nstr(pole_form),
)
print(
    "  arch           =",
    nstr(arch_matrix),
)

print()
print("Explicit Archimedean:")
print(
    "  K-fourier      =",
    nstr(explicit_arch),
)

print()
print("Agreement:")
print(
    "  matrix arch    =",
    nstr(arch_matrix),
)
print(
    "  explicit arch  =",
    nstr(explicit_arch),
)
print(
    "  difference     =",
    nstr(difference),
)
print(
    "  relative       =",
    nstr(relative_difference),
)

print()
print("Historical Cell-5 expression:")
print(
    "  linear result  =",
    nstr(historical_linear),
)
print(
    "  [WRONG; retained for forensic comparison]"
)

print()
print("Timing:")
print(
    f"  ground state   = {ground_elapsed:.6f} s"
)
print(
    f"  Q(T={T})       = {matrix_elapsed:.6f} s"
)
print(
    f"  total form     = {total_elapsed:.6f} s"
)
print(
    f"  prime          = {prime_elapsed:.6f} s"
)
print(
    f"  pole           = {pole_elapsed:.6f} s"
)
print(
    f"  Archimedean    = {arch_elapsed:.6f} s"
)
print(
    f"  linear         = {linear_elapsed:.6f} s"
)

print()
print("=" * 72)
print("CELL 21 COMPLETE")
print("=" * 72)
