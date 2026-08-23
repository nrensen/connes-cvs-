import mpmath as mp

from connes_cvs import build_galerkin_matrix, compute_ground_state
from connes_cvs.operator import _compute_psi_pair


# ============================================================
# CELL5_CORRECTED2A
#
# Full corrected Archimedean audit.
#
# Purpose:
#
#   1. Reproduce the repository Galerkin matrix.
#   2. Isolate its prime, pole and Archimedean pieces.
#   3. Build the Archimedean matrix independently from
#      _compute_psi_pair().
#   4. Reconstruct the ground-state Archimedean functional
#      directly from the source.
#   5. Separate the zero-mode contribution from the
#      nonzero-mode contribution.
#
# IMPORTANT:
#
# We do NOT use the old nested
#
#       r integral
#          -> D_v(r)
#             -> basis loop
#                -> S(r,x) quadrature
#
# construction.
#
# Instead the ground-state source is assembled once for each
# required x-value and then integrated as a scalar function.
#
# This is intended to give us an independent check of the
# matrix/dictionary correspondence without the catastrophic
# runtime of cell5_corrected.py.
# ============================================================


# ============================================================
# Parameters
# ============================================================

mp.mp.dps = 80

c = 13
N = 8
T = 60
dps = 80

L = mp.log(c)
alpha = 2 * mp.pi / L


print("=" * 70)
print("CELL5_CORRECTED2A — FULL CORRECTED ARCHIMEDEAN AUDIT")
print("=" * 70)

print()
print("Parameters:")
print("c =", c)
print("N =", N)
print("T =", T)
print("dps =", dps)
print("L =", mp.nstr(L, 60))
print("2*pi/L =", mp.nstr(alpha, 60))


# ============================================================
# 1. Repository Galerkin matrix
# ============================================================

print()
print("-" * 70)
print("1. BUILD REPOSITORY GALERKIN MATRIX")
print("-" * 70)

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=dps,
)

lam_min, eigvec = compute_ground_state(Q)

print()
print("lambda_min =")
print(mp.nstr(lam_min, 50))


# ============================================================
# 2. Ground-state vector
# ============================================================

coefficients = [
    mp.mpf(eigvec[i, 0])
    for i in range(eigvec.rows)
]

norm = mp.sqrt(sum(x * x for x in coefficients))

coefficients = [
    x / norm
    for x in coefficients
]

#
# Convert full complex-index representation into the
# canonical real-even coefficients:
#
#   v_0 = c_0
#   v_k = sqrt(2) c_k
#

v_ground = [coefficients[N]]

for k in range(1, N + 1):
    v_ground.append(
        mp.sqrt(2) * coefficients[N + k]
    )


#
# Full-space coefficient vector.
#

u_ground = mp.matrix(2 * N + 1, 1)

for m in range(-N, N + 1):

    if m == 0:

        u_ground[m + N] = v_ground[0]

    else:

        u_ground[m + N] = (
            v_ground[abs(m)] / mp.sqrt(2)
        )


print()
print("||v_ground|| =")
print(
    mp.nstr(
        mp.sqrt(mp.fdot(v_ground, v_ground)),
        40
    )
)

print()
print("||u_ground|| =")
print(
    mp.nstr(
        mp.sqrt(mp.fdot(u_ground, u_ground)),
        40
    )
)


# ============================================================
# 3. Prime-power list
# ============================================================

def prime_power_terms(c):

    c_int = int(mp.floor(c))

    terms = []

    for p in range(2, c_int + 1):

        is_prime = True

        for d in range(
            2,
            int(mp.sqrt(p)) + 1
        ):

            if p % d == 0:

                is_prime = False
                break

        if not is_prime:
            continue

        q = p

        while q <= c_int:

            terms.append(
                (
                    mp.mpf(q),
                    mp.log(p)
                )
            )

            q *= p

    return terms


prime_terms = prime_power_terms(c)


# ============================================================
# 4. Explicit prime matrix
# ============================================================

def Q_prime_power(q, Lambda_q):

    size = 2 * N + 1

    Qq = mp.matrix(size, size)

    a = 1 - mp.log(q) / L

    prefactor = (
        -1
        / mp.pi
        * Lambda_q
        / mp.sqrt(q)
    )

    def psi(x):

        return (
            prefactor
            * mp.sin(
                2 * mp.pi * x * a
            )
        )

    def psi_prime(x):

        return (
            prefactor
            * 2 * mp.pi * a
            * mp.cos(
                2 * mp.pi * x * a
            )
        )

    for i, m in enumerate(
        range(-N, N + 1)
    ):

        for j, n in enumerate(
            range(-N, N + 1)
        ):

            if m == n:

                Qq[i, j] = psi_prime(m)

            else:

                Qq[i, j] = (
                    psi(m) - psi(n)
                ) / mp.mpf(m - n)

    return Qq


Q_prime = mp.matrix(
    2 * N + 1,
    2 * N + 1
)

for q, Lambda_q in prime_terms:

    Q_prime += Q_prime_power(
        q,
        Lambda_q
    )


# ============================================================
# 5. Pole matrix
# ============================================================

def psi_pole(x):

    x = mp.mpf(x)

    f = lambda y: (
        2
        * mp.cosh(y / 2)
        * mp.sin(
            2 * mp.pi
            * x
            * (1 - y / L)
        )
    )

    return (
        mp.quad(
            f,
            [0, L]
        )
        / mp.pi
    )


def psi_pole_prime(x):

    x = mp.mpf(x)

    f = lambda y: (
        2
        * mp.cosh(y / 2)
        * (
            2 * mp.pi
            * (1 - y / L)
        )
        * mp.cos(
            2 * mp.pi
            * x
            * (1 - y / L)
        )
    )

    return (
        mp.quad(
            f,
            [0, L]
        )
        / mp.pi
    )


Q_pole = mp.matrix(
    2 * N + 1,
    2 * N + 1
)

for i, m in enumerate(
    range(-N, N + 1)
):

    for j, n in enumerate(
        range(-N, N + 1)
    ):

        if m == n:

            Q_pole[i, j] = (
                psi_pole_prime(m)
            )

        else:

            Q_pole[i, j] = (
                psi_pole(m)
                - psi_pole(n)
            ) / mp.mpf(m - n)


# ============================================================
# 6. Repository Archimedean matrix
# ============================================================

Q_arch_repo = (
    Q
    - Q_prime
    - Q_pole
)


# ============================================================
# 7. Explicit Archimedean psi values
#
# Empty prime-data list isolates the Archimedean source.
# ============================================================

print()
print("-" * 70)
print("7. BUILD EXPLICIT ARCHIMEDEAN PSI DATA")
print("-" * 70)

psi_arch = {}
psi_arch_prime = {}

for n in range(-N, N + 1):

    p, pd = _compute_psi_pair(
        n,
        L,
        T,
        dps,
        []
    )

    psi_arch[n] = mp.mpf(p)
    psi_arch_prime[n] = mp.mpf(pd)


# ============================================================
# 8. Explicit Archimedean matrix
# ============================================================

Q_arch_explicit = mp.matrix(
    2 * N + 1,
    2 * N + 1
)

for i, m in enumerate(
    range(-N, N + 1)
):

    for j, n in enumerate(
        range(-N, N + 1)
    ):

        if m == n:

            Q_arch_explicit[i, j] = (
                psi_arch_prime[n]
            )

        else:

            Q_arch_explicit[i, j] = (
                psi_arch[m]
                - psi_arch[n]
            ) / mp.mpf(m - n)


# ============================================================
# 9. Matrix comparison
# ============================================================

print()
print("-" * 70)
print("9. MATRIX COMPARISON")
print("-" * 70)

max_err = mp.mpf(0)
max_m = None
max_n = None

for i in range(2 * N + 1):

    for j in range(2 * N + 1):

        err = abs(
            Q_arch_repo[i, j]
            - Q_arch_explicit[i, j]
        )

        if err > max_err:

            max_err = err
            max_m = i - N
            max_n = j - N


print()
print("Maximum |Q_arch_repo - Q_arch_explicit| =")
print(mp.nstr(max_err, 50))

print()
print("Location:")
print("m =", max_m, "n =", max_n)


print()
print("Selected entries:")

for m, n in [
    (0, 0),
    (0, 1),
    (1, 1),
    (1, 2),
    (N, N),
]:

    i = m + N
    j = n + N

    print()
    print(f"({m},{n})")

    print(
        "repo     =",
        mp.nstr(
            Q_arch_repo[i, j],
            40
        )
    )

    print(
        "explicit =",
        mp.nstr(
            Q_arch_explicit[i, j],
            40
        )
    )

    print(
        "error    =",
        mp.nstr(
            abs(
                Q_arch_repo[i, j]
                - Q_arch_explicit[i, j]
            ),
            25
        )
    )


# ============================================================
# 10. Ground-state quadratic forms
# ============================================================

print()
print("-" * 70)
print("10. GROUND-STATE QUADRATIC FORMS")
print("-" * 70)

A_repo = (
    u_ground.T
    * Q_arch_repo
    * u_ground
)[0]

A_explicit = (
    u_ground.T
    * Q_arch_explicit
    * u_ground
)[0]


print()
print("u^T Q_arch_repo u =")
print(mp.nstr(A_repo, 50))

print()
print("u^T Q_arch_explicit u =")
print(mp.nstr(A_explicit, 50))

print()
print("|difference| =")
print(
    mp.nstr(
        abs(A_repo - A_explicit),
        50
    )
)


# ============================================================
# 11. Canonical basis audit
# ============================================================

def canonical_to_full(k):

    u = mp.matrix(
        2 * N + 1,
        1
    )

    if k == 0:

        u[N] = 1

    else:

        u[N + k] = (
            1 / mp.sqrt(2)
        )

        u[N - k] = (
            1 / mp.sqrt(2)
        )

    return u


print()
print("-" * 70)
print("11. CANONICAL BASIS AUDIT")
print("-" * 70)

for k in range(
    0,
    min(N, 4) + 1
):

    u = canonical_to_full(k)

    a_repo = (
        u.T
        * Q_arch_repo
        * u
    )[0]

    a_exp = (
        u.T
        * Q_arch_explicit
        * u
    )[0]

    print()
    print("k =", k)

    print(
        "repo     =",
        mp.nstr(a_repo, 35)
    )

    print(
        "explicit =",
        mp.nstr(a_exp, 35)
    )

    print(
        "error    =",
        mp.nstr(
            abs(a_repo - a_exp),
            25
        )
    )


# ============================================================
# 12. DIRECT SOURCE CONSTRUCTION
#
# This is the important new part.
#
# We construct the Archimedean source corresponding to the
# ground-state test function directly.
#
# The source used here is the h_+ / digamma representation
# audited independently in Cell 6.
#
# For a fixed r, the ground-state source is assembled as a
# linear combination of the already-required basis responses.
#
# There is no inner basis loop inside the outer quadrature.
# ============================================================


def h_plus(t):

    t = mp.mpf(t)

    return mp.mpf("0.5") * (
        mp.digamma(
            (1 + 1j * t) / 4
        )
        + mp.digamma(
            (1 - 1j * t) / 4
        )
    ).real


#
# Source S(r,x).
#
# We deliberately keep this as a separate function so that
# its result can be inspected independently.
#

def S_direct(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    a = 2 * mp.pi * x / L

    f = lambda t: (
        h_plus(
            r + t
        )
        * mp.sin(a * t)
    )

    #
    # The integrand is even in t for the relevant combination,
    # so use the positive half and double it.
    #
    return 2 * mp.quad(
        f,
        [0, T]
    )


# ============================================================
# 13. Direct source values
#
# These are the quantities which enter the divided-difference
# representation.
# ============================================================

print()
print("-" * 70)
print("13. DIRECT SOURCE VALUES")
print("-" * 70)

S_direct_values = {}

for x in range(-N, N + 1):

    S_direct_values[x] = S_direct(
        0,
        x
    )

    print(
        "x =",
        x,
        " S_direct(0,x) =",
        mp.nstr(
            S_direct_values[x],
            35
        )
    )


# ============================================================
# 14. Direct source-derived Archimedean matrix
#
# For off-diagonal entries:
#
#       (S(m)-S(n))/(m-n)
#
# For the diagonal we differentiate S with respect to x.
#
# This section is intentionally evaluated only at the discrete
# basis points. It is therefore cheap compared with the old
# nested cell5 calculation.
# ============================================================

print()
print("-" * 70)
print("14. DIRECT SOURCE-DERIVED MATRIX")
print("-" * 70)


Q_source_direct = mp.matrix(
    2 * N + 1,
    2 * N + 1
)


for i, m in enumerate(
    range(-N, N + 1)
):

    for j, n in enumerate(
        range(-N, N + 1)
    ):

        if m != n:

            Q_source_direct[i, j] = (
                S_direct(
                    0,
                    m
                )
                -
                S_direct(
                    0,
                    n
                )
            ) / mp.mpf(m - n)

        else:

            #
            # Numerical derivative only at the diagonal.
            #
            h = mp.sqrt(mp.eps)

            Q_source_direct[i, j] = (
                S_direct(
                    0,
                    mp.mpf(m) + h
                )
                -
                S_direct(
                    0,
                    mp.mpf(m) - h
                )
            ) / (2 * h)


# ============================================================
# 15. Direct-source matrix comparison
# ============================================================

print()
print("-" * 70)
print("15. DIRECT-SOURCE MATRIX COMPARISON")
print("-" * 70)

max_source_err = mp.mpf(0)
source_m = None
source_n = None

for i in range(2 * N + 1):

    for j in range(2 * N + 1):

        err = abs(
            Q_source_direct[i, j]
            - Q_arch_explicit[i, j]
        )

        if err > max_source_err:

            max_source_err = err
            source_m = i - N
            source_n = j - N


print()
print(
    "Maximum |Q_source_direct - Q_arch_explicit| ="
)

print(
    mp.nstr(
        max_source_err,
        50
    )
)

print()
print("Location:")
print(
    "m =",
    source_m,
    "n =",
    source_n
)


# ============================================================
# 16. Direct-source ground-state quadratic form
# ============================================================

print()
print("-" * 70)
print("16. DIRECT-SOURCE GROUND-STATE FORM")
print("-" * 70)

A_source = (
    u_ground.T
    * Q_source_direct
    * u_ground
)[0]


print()
print("u^T Q_source_direct u =")
print(
    mp.nstr(
        A_source,
        50
    )
)

print()
print(
    "|source - explicit| ="
)

print(
    mp.nstr(
        abs(
            A_source
            - A_explicit
        ),
        50
    )
)


# ============================================================
# 17. Zero-mode decomposition
#
# This is especially important because the previous audit
# showed the largest discrepancy at (0,0).
#
# Split the quadratic form into:
#
#       zero-zero
#       zero/nonzero cross terms
#       nonzero/nonzero
# ============================================================

print()
print("-" * 70)
print("17. ZERO-MODE DECOMPOSITION")
print("-" * 70)


zero = mp.matrix(
    2 * N + 1,
    1
)

nonzero = mp.matrix(
    2 * N + 1,
    1
)

for i in range(2 * N + 1):

    if i == N:

        zero[i] = u_ground[i]
        nonzero[i] = 0

    else:

        zero[i] = 0
        nonzero[i] = u_ground[i]


A00 = (
    zero.T
    * Q_arch_explicit
    * zero
)[0]

A0n = (
    zero.T
    * Q_arch_explicit
    * nonzero
)[0]

An0 = (
    nonzero.T
    * Q_arch_explicit
    * zero
)[0]

Ann = (
    nonzero.T
    * Q_arch_explicit
    * nonzero
)[0]


print()
print("zero-zero =")
print(mp.nstr(A00, 50))

print()
print("zero-nonzero =")
print(mp.nstr(A0n, 50))

print()
print("nonzero-zero =")
print(mp.nstr(An0, 50))

print()
print("nonzero-nonzero =")
print(mp.nstr(Ann, 50))

print()
print("sum =")
print(
    mp.nstr(
        A00 + A0n + An0 + Ann,
        50
    )
)

print()
print("full explicit quadratic form =")
print(
    mp.nstr(
        A_explicit,
        50
    )
)


# ============================================================
# 18. Parity / symmetry checks
# ============================================================

print()
print("-" * 70)
print("18. SYMMETRY / PARITY")
print("-" * 70)

sym_err = mp.mpf(0)

for i in range(2 * N + 1):

    for j in range(2 * N + 1):

        sym_err = max(
            sym_err,
            abs(
                Q_arch_explicit[i, j]
                -
                Q_arch_explicit[j, i]
            )
        )


parity_err = mp.mpf(0)

for n in range(1, N + 1):

    parity_err = max(
        parity_err,
        abs(
            psi_arch[-n]
            + psi_arch[n]
        )
    )

    parity_err = max(
        parity_err,
        abs(
            psi_arch_prime[-n]
            - psi_arch_prime[n]
        )
    )


print()
print("Maximum matrix symmetry error =")
print(
    mp.nstr(
        sym_err,
        30
    )
)

print()
print("Maximum psi parity error =")
print(
    mp.nstr(
        parity_err,
        30
    )
)


# ============================================================
# 19. Final summary
# ============================================================

print()
print("=" * 70)
print("CELL5_CORRECTED2A SUMMARY")
print("=" * 70)

print()
print("Maximum repository/explicit matrix error:")
print(
    mp.nstr(
        max_err,
        50
    )
)

print()
print("Ground-state repository/explicit difference:")
print(
    mp.nstr(
        abs(
            A_repo
            - A_explicit
        ),
        50
    )
)

print()
print("Maximum direct-source/explicit matrix error:")
print(
    mp.nstr(
        max_source_err,
        50
    )
)

print()
print("Direct-source ground-state form:")
print(
    mp.nstr(
        A_source,
        50
    )
)

print()
print("Direct-source minus explicit:")
print(
    mp.nstr(
        abs(
            A_source
            - A_explicit
        ),
        50
    )
)

print()
print("Zero-zero contribution:")
print(
    mp.nstr(
        A00,
        50
    )
)

print()
print("Zero/nonzero contribution:")
print(
    mp.nstr(
        A0n + An0,
        50
    )
)

print()
print("Nonzero/nonzero contribution:")
print(
    mp.nstr(
        Ann,
        50
    )
)

print()
print("=" * 70)
print("END CELL5_CORRECTED2A")
print("=" * 70)
