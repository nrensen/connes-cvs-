# ============================================================
# CELL_6 — ARCHIMEDEAN SOURCE / DICTIONARY AUDIT
#
# This cell uses the common definitions in cell.py.
#
# Goal:
#
#   Independently reconstruct the Archimedean matrix from the
#   source function
#
#       psi_R,T(x)
#         = 1/(2*pi^2) int_{-T}^T
#               h_+(r) S(r,x,L) dr
#
# where
#
#       S(r,x,L)
#         = int_0^L
#             sin(2*pi*x*(1-y/L)) cos(r y) dy.
#
# We then form the divided-difference matrix
#
#       Q_arch,ij =
#           (psi(i)-psi(j))/(i-j),  i != j
#
#       Q_arch,ii = psi'(i)
#
# and compare:
#
#   1. Source-derived Archimedean matrix
#   2. Repository Archimedean matrix
#   3. Explicit Weil-side Archimedean quadratic form
#
# No calculations are performed by cell.py.
# ============================================================

import mpmath as mp

from cell import (
    DEFAULT_DPS,
    compute_L,
    compute_beta,
    compute_delta,
    build_galerkin_matrix,
    compute_ground_state,
    normalise_ground_state,
    canonical_to_full,
    canonical_norm,
    prime_power_terms,
    sum_v_G,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = DEFAULT_DPS

c = 13
N = 8

# Keep this modest initially.  The source calculation is
# expensive because it contains nested quadrature.
T = 40

dps = DEFAULT_DPS

L = compute_L(c)
beta = compute_beta(L)


print("\n" + "=" * 60)
print("CELL_6 — ARCHIMEDEAN SOURCE / DICTIONARY AUDIT")
print("=" * 60)

print("\nParameters:")
print("N =", N)
print("c =", c)
print("T =", T)
print("L =", mp.nstr(L, 50))
print("dps =", dps)


# ============================================================
# 1. BUILD GROUND STATE
# ============================================================

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=dps,
)

lam_min, eigvec = compute_ground_state(Q)

v_star = normalise_ground_state(
    eigvec,
    N,
)

u_star = canonical_to_full(
    v_star,
    N,
)

print("\n" + "-" * 60)
print("1. GROUND STATE")
print("-" * 60)

print("\nlambda_min =")
print(mp.nstr(lam_min, 50))

print("\nCanonical dimension =")
print(len(v_star))

print("\nFull dimension =")
print(u_star.rows)

print("\n||v_star|| =")
print(mp.nstr(
    canonical_norm(v_star),
    50
))

print("\n||u_star|| =")
print(mp.nstr(
    mp.sqrt(mp.fdot(u_star, u_star)),
    50
))


# ============================================================
# 2. PRIME MATRIX
#
# Reconstruct the prime matrix independently from the
# divided-difference definition.
# ============================================================

def psi_prime_cell6(x, q, Lambda_q):

    x = mp.mpf(x)

    a = (
        1
        - mp.log(q) / L
    )

    prefactor = (
        -1
        / mp.pi
        * Lambda_q
        / mp.sqrt(q)
    )

    return (
        prefactor
        * mp.sin(
            2 * mp.pi * x * a
        )
    )


def psi_prime_derivative_cell6(x, q, Lambda_q):

    x = mp.mpf(x)

    a = (
        1
        - mp.log(q) / L
    )

    prefactor = (
        -1
        / mp.pi
        * Lambda_q
        / mp.sqrt(q)
    )

    return (
        prefactor
        * 2
        * mp.pi
        * a
        * mp.cos(
            2 * mp.pi * x * a
        )
    )


def build_prime_matrix_cell6():

    size = 2 * N + 1

    Q_prime = mp.matrix(
        size,
        size
    )

    terms = prime_power_terms(c)

    for q, Lambda_q in terms:

        values = {}

        derivatives = {}

        for x in range(-N, N + 1):

            values[x] = psi_prime_cell6(
                x,
                q,
                Lambda_q
            )

            derivatives[x] = (
                psi_prime_derivative_cell6(
                    x,
                    q,
                    Lambda_q
                )
            )

        for i, m in enumerate(
            range(-N, N + 1)
        ):

            for j, n in enumerate(
                range(-N, N + 1)
            ):

                if m != n:

                    Q_prime[i, j] += (
                        values[m]
                        - values[n]
                    ) / mp.mpf(m - n)

                else:

                    Q_prime[i, j] += (
                        derivatives[m]
                    )

    return Q_prime


Q_prime = build_prime_matrix_cell6()

prime_form = mp.fdot(
    u_star,
    Q_prime * u_star
)

print("\n" + "-" * 60)
print("2. PRIME MATRIX")
print("-" * 60)

print("\nPrime quadratic form =")
print(mp.nstr(
    prime_form,
    60
))


# ============================================================
# 3. POLE MATRIX
# ============================================================

def psi_pole_cell6(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * mp.sin(
            2
            * mp.pi
            * x
            * (1 - y / L)
        )
    )

    return (
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L]
        )
    )


def psi_pole_derivative_cell6(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * (
            2
            * mp.pi
            * (1 - y / L)
        )
        * mp.cos(
            2
            * mp.pi
            * x
            * (1 - y / L)
        )
    )

    return (
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L]
        )
    )


def build_pole_matrix_cell6():

    size = 2 * N + 1

    Q_pole = mp.matrix(
        size,
        size
    )

    values = {}
    derivatives = {}

    for x in range(-N, N + 1):

        values[x] = psi_pole_cell6(x)

        derivatives[x] = (
            psi_pole_derivative_cell6(x)
        )

    for i, m in enumerate(
        range(-N, N + 1)
    ):

        for j, n in enumerate(
            range(-N, N + 1)
        ):

            if m != n:

                Q_pole[i, j] = (
                    values[m]
                    - values[n]
                ) / mp.mpf(m - n)

            else:

                Q_pole[i, j] = (
                    derivatives[m]
                )

    return Q_pole


Q_pole = build_pole_matrix_cell6()

pole_form = mp.fdot(
    u_star,
    Q_pole * u_star
)

pole_explicit = mp.re(
    2 * sum_v_G(
        v_star,
        1j / 2,
        L,
    )
)

print("\n" + "-" * 60)
print("3. POLE DICTIONARY CHECK")
print("-" * 60)

print("\nPole matrix form =")
print(mp.nstr(
    pole_form,
    60
))

print("\n2 G(i/2) =")
print(mp.nstr(
    pole_explicit,
    60
))

print("\nDifference =")
print(mp.nstr(
    pole_form - pole_explicit,
    60
))


# ============================================================
# 4. ARCHIMEDEAN SOURCE FUNCTIONS
#
# h_+(r)
#
# and
#
# S(r,x,L).
# ============================================================

def h_plus_cell6(r):

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


def S_arch_cell6(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    integrand = lambda y: (
        mp.sin(
            2
            * mp.pi
            * x
            * (1 - y / L)
        )
        * mp.cos(r * y)
    )

    return mp.quad(
        integrand,
        [0, L]
    )


# ============================================================
# 5. SOURCE FUNCTION psi_R,T(x)
#
# The integrand is even in r, so
#
#   1/(2*pi^2) int_-T^T
#
# becomes
#
#   1/pi^2 int_0^T.
#
# ============================================================

def psi_arch_T_cell6(x, T_local):

    x = mp.mpf(x)
    T_local = mp.mpf(T_local)

    integrand = lambda r: (
        h_plus_cell6(r)
        * S_arch_cell6(r, x)
    )

    # Split the integration range into moderate intervals.
    #
    # This is deliberately explicit rather than relying on a
    # single enormous interval for the oscillatory integral.

    step = mp.mpf(5)

    points = [mp.mpf(0)]

    r = step

    while r < T_local:

        points.append(r)

        r += step

    points.append(T_local)

    integral = mp.quad(
        integrand,
        points
    )

    return (
        integral
        / (mp.pi ** 2)
    )


# ============================================================
# 6. SOURCE VALUES AT ALL INTEGER FOURIER INDICES
#
# We calculate these ONCE.
#
# This is important: later matrix construction uses the
# cached values rather than repeatedly evaluating the nested
# Archimedean integral.
# ============================================================

print("\n" + "-" * 60)
print("4. ARCHIMEDEAN SOURCE VALUES")
print("-" * 60)

print(
    "\nComputing psi_R,T(x) for"
    f" x = -{N}, ..., {N}"
)

psi_arch_values = {}

for x in range(-N, N + 1):

    value = psi_arch_T_cell6(
        x,
        T
    )

    psi_arch_values[x] = value

    print(
        "x =",
        f"{x:>3}",
        " psi_arch,T(x) =",
        mp.nstr(value, 40)
    )


# ============================================================
# 7. SOURCE-DERIVED ARCHIMEDEAN MATRIX
#
# IMPORTANT:
#
# The off-diagonal divided differences require only the source
# values above.
#
# The diagonal requires psi'(x).
#
# Rather than numerically differentiating the expensive source
# integral, calculate the derivative by differentiating the
# sine factor analytically.
# ============================================================

def S_arch_derivative_cell6(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.pi
        * (
            1 - y / L
        )
        * mp.cos(
            2
            * mp.pi
            * x
            * (1 - y / L)
        )
        * mp.cos(r * y)
    )

    return mp.quad(
        integrand,
        [0, L]
    )


def psi_arch_derivative_T_cell6(
    x,
    T_local,
):

    x = mp.mpf(x)
    T_local = mp.mpf(T_local)

    integrand = lambda r: (
        h_plus_cell6(r)
        * S_arch_derivative_cell6(
            r,
            x
        )
    )

    step = mp.mpf(5)

    points = [mp.mpf(0)]

    r = step

    while r < T_local:

        points.append(r)

        r += step

    points.append(T_local)

    integral = mp.quad(
        integrand,
        points
    )

    return (
        integral
        / (mp.pi ** 2)
    )


print("\nComputing Archimedean source derivatives...")

psi_arch_derivatives = {}

for x in range(-N, N + 1):

    value = psi_arch_derivative_T_cell6(
        x,
        T
    )

    psi_arch_derivatives[x] = value

    print(
        "x =",
        f"{x:>3}",
        " psi'_arch,T(x) =",
        mp.nstr(value, 40)
    )


def build_source_arch_matrix_cell6():

    size = 2 * N + 1

    Q_source = mp.matrix(
        size,
        size
    )

    for i, m in enumerate(
        range(-N, N + 1)
    ):

        for j, n in enumerate(
            range(-N, N + 1)
        ):

            if m != n:

                Q_source[i, j] = (
                    psi_arch_values[m]
                    - psi_arch_values[n]
                ) / mp.mpf(m - n)

            else:

                Q_source[i, j] = (
                    psi_arch_derivatives[m]
                )

    return Q_source


Q_source_arch = (
    build_source_arch_matrix_cell6()
)

source_arch_form = mp.fdot(
    u_star,
    Q_source_arch * u_star
)


# ============================================================
# 8. REPOSITORY ARCHIMEDEAN MATRIX
#
# This is the repository matrix with the independently
# reconstructed prime and pole matrices removed.
# ============================================================

Q_repository = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=dps,
)

Q_repository_arch = (
    Q_repository
    - Q_prime
    - Q_pole
)

repository_arch_form = mp.fdot(
    u_star,
    Q_repository_arch * u_star
)


# ============================================================
# 9. EXPLICIT WEIL-SIDE ARCHIMEDEAN QUADRATIC FORM
#
# A_T =
#
#   1/pi int_0^T h_+(r) g(r) dr
#
# where g(r) = Re G(r).
# ============================================================

def explicit_arch_form_cell6(T_local):

    T_local = mp.mpf(T_local)

    integrand = lambda r: (
        h_plus_cell6(r)
        * mp.re(
            sum_v_G(
                v_star,
                r,
                L,
            )
        )
    )

    step = mp.mpf(5)

    points = [mp.mpf(0)]

    r = step

    while r < T_local:

        points.append(r)

        r += step

    points.append(T_local)

    integral = mp.quad(
        integrand,
        points
    )

    return integral / mp.pi


explicit_arch_form = (
    explicit_arch_form_cell6(T)
)


# ============================================================
# 10. DICTIONARY COMPARISON
# ============================================================

print("\n" + "-" * 60)
print("5. ARCHIMEDEAN DICTIONARY COMPARISON")
print("-" * 60)

print("\nT =")
print(T)

print("\nSource-derived matrix form =")
print(mp.nstr(
    source_arch_form,
    60
))

print("\nRepository Archimedean matrix form =")
print(mp.nstr(
    repository_arch_form,
    60
))

print("\nExplicit Weil Archimedean form =")
print(mp.nstr(
    explicit_arch_form,
    60
))

print("\nSource - Repository =")
print(mp.nstr(
    source_arch_form
    - repository_arch_form,
    60
))

print("\nSource - Explicit =")
print(mp.nstr(
    source_arch_form
    - explicit_arch_form,
    60
))

print("\nRepository - Explicit =")
print(mp.nstr(
    repository_arch_form
    - explicit_arch_form,
    60
))


# ============================================================
# 11. FULL MATRIX COMPARISON
#
# Compare the source-derived Archimedean matrix with the
# repository Archimedean matrix element by element.
#
# This is particularly useful because a quadratic-form
# discrepancy can otherwise conceal a sign, factor, or
# transpose/convention error.
# ============================================================

print("\n" + "-" * 60)
print("6. MATRIX-LEVEL COMPARISON")
print("-" * 60)

max_abs_source_repository = mp.mpf(0)
max_i = None
max_j = None

for i in range(2 * N + 1):

    for j in range(2 * N + 1):

        diff = abs(
            Q_source_arch[i, j]
            - Q_repository_arch[i, j]
        )

        if diff > max_abs_source_repository:

            max_abs_source_repository = diff
            max_i = i
            max_j = j


print("\nMaximum absolute matrix difference =")
print(mp.nstr(
    max_abs_source_repository,
    60
))

if max_i is not None:

    m_max = max_i - N
    n_max = max_j - N

    print("\nLocation of maximum difference:")
    print(
        "m =",
        m_max,
        " n =",
        n_max
    )

    print("\nSource value =")
    print(mp.nstr(
        Q_source_arch[max_i, max_j],
        50
    ))

    print("\nRepository value =")
    print(mp.nstr(
        Q_repository_arch[max_i, max_j],
        50
    ))


# ============================================================
# 12. SYMMETRY CHECKS
#
# The Archimedean divided-difference matrix should be
# symmetric in this real formulation.
# ============================================================

def matrix_symmetry_error(M):

    size = M.rows

    maximum = mp.mpf(0)

    for i in range(size):

        for j in range(size):

            err = abs(
                M[i, j]
                - M[j, i]
            )

            if err > maximum:

                maximum = err

    return maximum


print("\n" + "-" * 60)
print("7. SYMMETRY CHECKS")
print("-" * 60)

print("\nSource matrix symmetry error =")
print(mp.nstr(
    matrix_symmetry_error(
        Q_source_arch
    ),
    40
))

print("\nRepository Archimedean symmetry error =")
print(mp.nstr(
    matrix_symmetry_error(
        Q_repository_arch
    ),
    40
))


# ============================================================
# 13. REQUIRED ARCHIMEDEAN VALUE
#
# Since
#
#     Q = Q_prime + Q_pole + Q_arch
#
# the Archimedean contribution required by the complete
# Galerkin quadratic form is
#
#     Q_required =
#         Q_form - prime_form - pole_form.
#
# ============================================================

Q_form = mp.fdot(
    u_star,
    Q_repository * u_star
)

required_arch = (
    Q_form
    - prime_form
    - pole_form
)


print("\n" + "-" * 60)
print("8. REQUIRED ARCHIMEDEAN VALUE")
print("-" * 60)

print("\nComplete Galerkin Q =")
print(mp.nstr(
    Q_form,
    60
))

print("\nPrime =")
print(mp.nstr(
    prime_form,
    60
))

print("\nPole =")
print(mp.nstr(
    pole_form,
    60
))

print("\nRequired Archimedean =")
print(mp.nstr(
    required_arch,
    60
))

print("\nSource-derived Archimedean - required =")
print(mp.nstr(
    source_arch_form
    - required_arch,
    60
))

print("\nRepository Archimedean - required =")
print(mp.nstr(
    repository_arch_form
    - required_arch,
    60
))

print("\nExplicit Archimedean - required =")
print(mp.nstr(
    explicit_arch_form
    - required_arch,
    60
))


# ============================================================
# 14. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("CELL_6 SUMMARY")
print("=" * 60)

print("\nT =")
print(T)

print("\nComplete Galerkin Q =")
print(mp.nstr(
    Q_form,
    60
))

print("\nPrime contribution =")
print(mp.nstr(
    prime_form,
    60
))

print("\nPole contribution =")
print(mp.nstr(
    pole_form,
    60
))

print("\nRequired Archimedean =")
print(mp.nstr(
    required_arch,
    60
))

print("\nSource-derived Archimedean =")
print(mp.nstr(
    source_arch_form,
    60
))

print("\nRepository Archimedean =")
print(mp.nstr(
    repository_arch_form,
    60
))

print("\nExplicit Weil Archimedean =")
print(mp.nstr(
    explicit_arch_form,
    60
))

print("\nSource - Repository =")
print(mp.nstr(
    source_arch_form
    - repository_arch_form,
    60
))

print("\nSource - Explicit =")
print(mp.nstr(
    source_arch_form
    - explicit_arch_form,
    60
))

print("\nRepository - Explicit =")
print(mp.nstr(
    repository_arch_form
    - explicit_arch_form,
    60
))

print("\n" + "=" * 60)
print("END CELL_6")
print("=" * 60)
