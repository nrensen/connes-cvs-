# ============================================================
# CELL 24 — FINITE-T ARCHIMEDEAN CONVERGENCE MAP
#
# Purpose
# -------
# Investigate the finite-T behaviour of the Archimedean
# quadratic functional using the validated Cell-23 analytic
# reduction.
#
# IMPORTANT:
#
# The forensic ground state is held FIXED.
#
# The ground state is the canonical forensic ground state:
#
#     c = 13
#     N = 8
#     T = 400
#     generation_dps = 150
#
# The T varied in this cell is ONLY the upper limit of the
# Archimedean r-integral:
#
#     A_arch(T)
#       = 1/pi * int_0^T h_+(r) J_v(r) dr
#
# No Galerkin matrix is rebuilt as T changes.
#
# This isolates finite-T behaviour of the Archimedean functional
# for one fixed ground state.
#
# No numerical y-integration is performed.
# ============================================================

import time

import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    canonical_to_full,
    full_to_canonical,
    compute_L,
)

# ============================================================
# PARAMETERS
# ============================================================

WORKING_DPS = 40

# Initial finite-T survey.
#
# The relatively dense low-T grid is intended to show the
# shape of convergence rather than merely provide isolated
# endpoint values.
T_VALUES = (
    list(range(5, 105, 5))
    + [120, 150, 200, 250, 300, 350, 400]
)

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

DISPLAY_DIGITS = WORKING_DPS


def nstr(x):
    return mp.nstr(
        x,
        DISPLAY_DIGITS,
    )


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# HEADER
# ============================================================

print("=" * 78)
print("CELL 24 — FINITE-T ARCHIMEDEAN CONVERGENCE MAP")
print("=" * 78)

print()
print("Parameters:")
print(f"  c              = {c}")
print(f"  N              = {N}")
print(f"  working_dps    = {WORKING_DPS}")

print()
print("T values:")
print(f"  {T_VALUES}")

print()
print("Forensic ground state:")
print(
    f"  c              = "
    f"{FORENSIC_GROUND_STATE['c']}"
)
print(
    f"  N              = "
    f"{FORENSIC_GROUND_STATE['N']}"
)
print(
    f"  T              = "
    f"{FORENSIC_GROUND_STATE['T']}"
)
print(
    f"  generation_dps = "
    f"{FORENSIC_GROUND_STATE['dps']}"
)

print()
print("L =")
print(nstr(L))


# ============================================================
# 1. FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 78)
print("1. FORENSIC GROUND STATE")
print("-" * 78)

ground_wall_start = time.perf_counter()
ground_cpu_start = time.process_time()

lambda_forensic, u_star, ground_meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
    verbose=True,
)

ground_wall_elapsed = elapsed(
    ground_wall_start
)
ground_cpu_elapsed = (
    time.process_time()
    - ground_cpu_start
)

u_star = mp.matrix(u_star)

print()
print("lambda_forensic =")
print(nstr(lambda_forensic))

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

print()
print(
    f"ground-state retrieval wall = "
    f"{ground_wall_elapsed:.6f} s"
)
print(
    f"ground-state retrieval CPU  = "
    f"{ground_cpu_elapsed:.6f} s"
)


# ============================================================
# 2. CANONICAL / FULL REPRESENTATIONS
# ============================================================

print()
print("-" * 78)
print("2. CANONICAL / FULL REPRESENTATIONS")
print("-" * 78)

v_star = full_to_canonical(
    u_star,
    N,
)

u = canonical_to_full(
    v_star,
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

print()
print("||u - u_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                u - u_star,
                u - u_star,
            )
        )
    )
)


# ============================================================
# 3. FOURIER MODE DATA
# ============================================================

modes = list(
    range(-N, N + 1)
)

a = {
    m: (
        2
        * mp.pi
        * m
        / L
    )
    for m in modes
}

coeff = {
    m: u[i]
    for i, m in enumerate(modes)
}


# ============================================================
# 4. STABLE ELEMENTARY FUNCTIONS
# ============================================================

def sinc(z):
    """
    sin(z) / z, with the removable value 1 at z=0.
    """

    z = mp.mpf(z)

    if z == 0:
        return mp.mpf("1")

    return mp.sin(z) / z


def one_minus_cos(x):
    """
    1 - cos(x), evaluated through 2 sin^2(x/2).
    """

    x = mp.mpf(x)

    return (
        2
        * mp.sin(x / 2) ** 2
    )


# ============================================================
# 5. FOURIER MODE FUNCTIONS
# ============================================================

def S_mode(m, r):
    """
    S_m(r) =
        int_0^L sin(a_m y) cos(r y) dy
    """

    am = a[m]

    if m == 0:
        return mp.mpf("0")

    k_plus = am + r
    k_minus = am - r

    plus = (
        one_minus_cos(
            k_plus * L
        )
        / k_plus
        if k_plus != 0
        else mp.mpf("0")
    )

    minus = (
        one_minus_cos(
            k_minus * L
        )
        / k_minus
        if k_minus != 0
        else mp.mpf("0")
    )

    return (
        plus + minus
    ) / 2


def W(k):
    """
    W(k) =
        int_0^L (1-y/L) cos(k y) dy

    expressed in stable sinc form.
    """

    k = mp.mpf(k)

    return (
        L
        / 2
        * sinc(
            k * L / 2
        ) ** 2
    )


def C_mode(m, r):
    """
    C_m(r) =
        int_0^L
            (1-y/L)
            cos(a_m y)
            cos(r y)
        dy
    """

    am = a[m]

    return (
        W(am - r)
        + W(am + r)
    ) / 2


# ============================================================
# 6. ANALYTIC J_v(r)
# ============================================================

def analytic_J(r):
    """
    Analytic quadratic Fourier kernel J_v(r).

    This is the Cell-23 implementation, retained here rather
    than imported from Cell 23 so that Cell 24 remains an
    independently executable research cell.
    """

    r = mp.mpf(r)

    # S_m(r) is shared by many pair terms.
    S = {
        m: S_mode(m, r)
        for m in modes
    }

    total = mp.mpf("0")

    # --------------------------------------------------------
    # Diagonal
    # --------------------------------------------------------

    for m in modes:

        um = coeff[m]

        total += (
            2
            * um
            * um
            * C_mode(m, r)
        )

    # --------------------------------------------------------
    # Off-diagonal
    #
    # Only m < n is required because the two orientations
    # contribute identically.
    # --------------------------------------------------------

    for i, m in enumerate(modes):

        um = coeff[m]

        for n in modes[i + 1:]:

            un = coeff[n]

            total += (
                2
                * um
                * un
                / mp.pi
                * (
                    S[n]
                    - S[m]
                )
                / mp.mpf(m - n)
            )

    return total


# ============================================================
# 7. ARCHIMEDEAN SOURCE
# ============================================================

def h_plus(r):
    """
    h_+(r) =
        Re psi(1/4 + i r / 2) - log(pi)
    """

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
# 8. FINITE-T ARCHIMEDEAN FUNCTIONAL
# ============================================================

def archimedean_to_T(T):
    """
    Compute

        A_arch(T)
          = 1/pi * int_0^T h_+(r) J_v(r) dr

    using the Cell-23 analytic reduction.
    """

    T = mp.mpf(T)

    return (
        mp.quad(
            lambda r:
                h_plus(r)
                * analytic_J(r),
            [0, T],
        )
        / mp.pi
    )


# ============================================================
# 9. FINITE-T SCAN
# ============================================================

print()
print("-" * 78)
print("3. FINITE-T ARCHIMEDEAN CONVERGENCE")
print("-" * 78)

print()
print(
    "The ground state is held fixed while only the upper "
    "limit T of the Archimedean r-integral is varied."
)

print()
print(
    "No Galerkin matrix is rebuilt as T changes."
)

print()
print(
    "Columns:"
)
print(
    "  T"
    "              A_arch(T)"
    "                          delta"
    "                 wall_s"
    "          cpu_s"
)

results = []

previous_value = None

total_wall_start = time.perf_counter()
total_cpu_start = time.process_time()

for T in T_VALUES:

    wall_start = time.perf_counter()
    cpu_start = time.process_time()

    value = archimedean_to_T(T)

    wall_elapsed = (
        time.perf_counter()
        - wall_start
    )

    cpu_elapsed = (
        time.process_time()
        - cpu_start
    )

    if previous_value is None:
        delta = mp.mpf("0")
    else:
        delta = value - previous_value

    results.append(
        (
            T,
            value,
            delta,
            wall_elapsed,
            cpu_elapsed,
        )
    )

    print(
        f"{T:6d}  "
        f"{nstr(value):>{DISPLAY_DIGITS + 4}}  "
        f"{nstr(delta):>{DISPLAY_DIGITS + 4}}  "
        f"{wall_elapsed:12.6f}  "
        f"{cpu_elapsed:12.6f}"
    )

    previous_value = value


total_wall_elapsed = (
    time.perf_counter()
    - total_wall_start
)

total_cpu_elapsed = (
    time.process_time()
    - total_cpu_start
)


# ============================================================
# 10. SUMMARY
# ============================================================

print()
print("-" * 78)
print("4. FINITE-T SUMMARY")
print("-" * 78)

print()
print("T                  A_arch(T)")

for T, value, delta, wall, cpu in results:

    print(
        f"{T:6d}   {nstr(value)}"
    )

print()
print(
    "Final T =",
    T_VALUES[-1],
)

print()
print(
    "A_arch(final) ="
)
print(
    nstr(
        results[-1][1]
    )
)

print()
print(
    f"finite-T scan wall = "
    f"{total_wall_elapsed:.6f} s"
)

print(
    f"finite-T scan CPU  = "
    f"{total_cpu_elapsed:.6f} s"
)


# ============================================================
# 11. FINAL RESULT
# ============================================================

print()
print("=" * 78)
print("CELL 24 COMPLETE")
print("=" * 78)
