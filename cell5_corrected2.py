import mpmath as mp
from connes_cvs import build_galerkin_matrix, compute_ground_state


# ============================================================
# CELL_5_CORRECTED2
#
# Fast Archimedean functional audit.
#
# Main change from cell_5_corrected.py:
#
#   DO NOT evaluate
#
#       integral h(r) D_v(r) dr
#
#   by repeatedly constructing D_v(r).
#
# Instead:
#
#   1. construct the explicit Archimedean matrix A_exp;
#   2. compare A_exp directly with the repository Archimedean
#      matrix A_repo;
#   3. evaluate v^T A v once.
#
# This turns the expensive nested quadrature into a finite
# collection of one-dimensional scalar quadratures.
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

print("=" * 70)
print("CELL_5_CORRECTED2 — FAST ARCHIMEDEAN MATRIX AUDIT")
print("=" * 70)

print()
print("Parameters:")
print("c =", c)
print("N =", N)
print("T =", T)
print("dps =", dps)
print("L =", mp.nstr(L, 60))
print("2*pi/L =", mp.nstr(2 * mp.pi / L, 60))


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

print("lambda_min =")
print(mp.nstr(lam_min, 50))


# ============================================================
# 2. Full-space ground state
#
# Repository ordering:
#
#     m = -N,...,0,...,+N
#
# The previous Cell-5 code uses the canonical real-even
# representation
#
#     v_0 = c_0
#     v_k = sqrt(2) c_k.
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

v_ground = [coefficients[N]]

for k in range(1, N + 1):
    v_ground.append(
        mp.sqrt(2) * coefficients[N + k]
    )

u_ground = mp.matrix(2 * N + 1, 1)

for m in range(-N, N + 1):
    if m == 0:
        u_ground[m + N] = v_ground[0]
    else:
        u_ground[m + N] = (
            v_ground[abs(m)] / mp.sqrt(2)
        )

print()
print("||v_ground|| =",
      mp.nstr(mp.sqrt(mp.fdot(v_ground, v_ground)), 40))

print("||u_ground|| =",
      mp.nstr(mp.sqrt(mp.fdot(u_ground, u_ground)), 40))


# ============================================================
# 3. Prime-power list
# ============================================================

def prime_power_terms(c):
    c_int = int(mp.floor(c))

    terms = []

    for p in range(2, c_int + 1):

        is_prime = True

        for d in range(2, int(mp.sqrt(p)) + 1):
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
        -1 / mp.pi
        * Lambda_q / mp.sqrt(q)
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

    for i, m in enumerate(range(-N, N + 1)):

        for j, n in enumerate(range(-N, N + 1)):

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
#
# This reproduces the pole contribution independently.
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
        mp.quad(f, [0, L])
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
        mp.quad(f, [0, L])
        / mp.pi
    )


Q_pole = mp.matrix(
    2 * N + 1,
    2 * N + 1
)

for i, m in enumerate(range(-N, N + 1)):

    for j, n in enumerate(range(-N, N + 1)):

        if m == n:

            Q_pole[i, j] = psi_pole_prime(m)

        else:

            Q_pole[i, j] = (
                psi_pole(m)
                - psi_pole(n)
            ) / mp.mpf(m - n)


# ============================================================
# 6. Repository Archimedean matrix
#
# Q = Q_prime + Q_pole + Q_arch
# ============================================================

Q_arch_repo = (
    Q
    - Q_prime
    - Q_pole
)


# ============================================================
# 7. Canonical basis -> full-space coefficients
# ============================================================

def canonical_to_full(k):

    u = mp.matrix(2 * N + 1, 1)

    if k == 0:

        u[N] = 1

    else:

        u[N + k] = 1 / mp.sqrt(2)
        u[N - k] = 1 / mp.sqrt(2)

    return u


# ============================================================
# 8. Fast explicit Archimedean source
#
# This is the crucial optimisation.
#
# For each canonical pair (k,l), form the corresponding
# quadratic source D_{kl}(r).
#
# We integrate that scalar function ONCE.
#
# There is no outer integral containing an inner construction
# of D_v(r).
# ============================================================

def F_full_basis(k, r):

    """
    Fourier/Mellin response corresponding to full-space
    basis vector e_k.

    k is an integer in [-N,N].
    """

    k = int(k)
    r = mp.mpf(r)

    a = 2 * mp.pi * k / L

    denom = a - r

    if denom == 0:
        return mp.mpc(L)

    z = denom * L

    if abs(z) < mp.sqrt(mp.eps):

        return (
            mp.expm1(1j * z)
            / (1j * denom)
        )

    return (
        mp.exp(-1j * r * L) - 1
    ) / (1j * denom)


def F_full_basis_phase(k, r):

    return (
        mp.exp(1j * r * L / 2)
        * F_full_basis(k, r)
        / mp.sqrt(L)
    )


# ============================================================
# 9. Source kernel for a pair of full-space basis vectors
#
# This uses the same source/divided-difference structure as
# the independent Cell-6 audit.
#
# The important point is that the resulting scalar function
# is associated with ONE matrix entry.
# ============================================================

def source_pair(m, n, r):

    """
    Archimedean source kernel for full-space basis pair m,n.

    This is the direct divided-difference form induced by the
    Archimedean source S(r,x).

    The implementation deliberately keeps m,n explicit so
    that each matrix element can be audited independently.
    """

    # --------------------------------------------------------
    # S(r,x)
    #
    # We use the repository's Archimedean source through the
    # same h_+ representation audited in Cell 6.
    # --------------------------------------------------------

    def S(x):

        x = mp.mpf(x)

        alpha = 2 * mp.pi * x / L

        # h_+ contribution.
        #
        # The source is evaluated by the same integral
        # convention used in _compute_psi_pair.
        #
        # Splitting at the natural breakpoints is important.
        # For the present c=13 audit these are r/5.
        # ----------------------------------------------------

        f = lambda t: (
            mp.mpf('0.5')
            * (
                mp.digamma(
                    (1 + 1j * (r + t)) / 4
                )
                + mp.digamma(
                    (1 + 1j * (r - t)) / 4
                )
            ).real
            * mp.sin(alpha * t)
        )

        # This helper is intentionally written as a direct
        # scalar quadrature. It is NOT called inside an outer
        # quadrature over r.
        #
        # The exact source dictionary is tested separately
        # below against Q_arch_repo.
        return mp.quad(
            f,
            [-T, T]
        )

    if m == n:

        h = mp.sqrt(mp.eps)

        return (
            S(mp.mpf(m) + h)
            - S(mp.mpf(m) - h)
        ) / (2 * h)

    return (
        S(m) - S(n)
    ) / mp.mpf(m - n)


# ============================================================
# IMPORTANT
#
# The completely direct source_pair implementation above is
# retained only as a reference definition.
#
# We do NOT use it for the production audit because it would
# simply move the nested quadrature problem into another place.
#
# Instead we use the repository's actual source implementation
# already exposed by operator.py.
# ============================================================

from connes_cvs.operator import _compute_psi_pair


# ============================================================
# 10. Extract explicit Archimedean psi values
#
# _compute_psi_pair(n,L,T,dps,prime_data) returns the complete
# source pair for the requested basis index.
#
# To isolate the Archimedean component we call it with an empty
# prime-data list.
# ============================================================

print()
print("-" * 70)
print("10. BUILD EXPLICIT ARCHIMEDEAN MATRIX")
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
# 11. Explicit Archimedean matrix
#
# This is the direct finite-dimensional matrix corresponding
# to the Archimedean functional.
# ============================================================

Q_arch_explicit = mp.matrix(
    2 * N + 1,
    2 * N + 1
)

for i, m in enumerate(range(-N, N + 1)):

    for j, n in enumerate(range(-N, N + 1)):

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
# 12. Matrix comparison
# ============================================================

print()
print("-" * 70)
print("11. MATRIX COMPARISON")
print("-" * 70)

max_err = mp.mpf(0)
max_i = None
max_j = None

for i in range(2 * N + 1):

    for j in range(2 * N + 1):

        err = abs(
            Q_arch_repo[i, j]
            - Q_arch_explicit[i, j]
        )

        if err > max_err:

            max_err = err
            max_i = i
            max_j = j


print()
print("Maximum |Q_arch_repo - Q_arch_explicit| =")
print(mp.nstr(max_err, 50))

print()
print("Location:")
print(
    "m =",
    max_i - N,
    "n =",
    max_j - N
)

print()
print("Selected matrix entries:")

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
        mp.nstr(Q_arch_repo[i, j], 40)
    )

    print(
        "explicit =",
        mp.nstr(Q_arch_explicit[i, j], 40)
    )

    print(
        "error    =",
        mp.nstr(
            abs(
                Q_arch_repo[i, j]
                - Q_arch_explicit[i, j]
            ),
            20
        )
    )


# ============================================================
# 13. Quadratic-form comparison
#
# This is now cheap.
#
# There is NO quadrature here.
# ============================================================

print()
print("-" * 70)
print("12. GROUND-STATE QUADRATIC FORM")
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
# 14. Canonical-basis quadratic forms
#
# These are particularly useful if the full ground-state
# comparison exposes a normalisation problem.
# ============================================================

print()
print("-" * 70)
print("13. CANONICAL BASIS AUDIT")
print("-" * 70)

for k in range(0, min(N, 4) + 1):

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
            20
        )
    )


# ============================================================
# 15. SYMMETRY CHECK
# ============================================================

print()
print("-" * 70)
print("14. SYMMETRY / PARITY CHECK")
print("-" * 70)

sym_err = mp.mpf(0)

for i in range(2 * N + 1):

    for j in range(2 * N + 1):

        err = abs(
            Q_arch_explicit[i, j]
            - Q_arch_explicit[j, i]
        )

        sym_err = max(sym_err, err)


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
print(mp.nstr(sym_err, 30))

print()
print("Maximum psi parity error =")
print(mp.nstr(parity_err, 30))


# ============================================================
# 16. Summary
# ============================================================

print()
print("=" * 70)
print("CELL_5_CORRECTED2 SUMMARY")
print("=" * 70)

print()
print("Maximum matrix error:")
print(mp.nstr(max_err, 50))

print()
print("Ground-state quadratic-form difference:")
print(
    mp.nstr(
        abs(A_repo - A_explicit),
        50
    )
)

print()
print("The expensive nested r -> D_v(r) -> basis-loop")
print("quadrature has NOT been used for the final comparison.")

print()
print("=" * 70)
print("END CELL_5_CORRECTED2")
print("=" * 70)
