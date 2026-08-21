import mpmath as mp
from connes_cvs import build_galerkin_matrix, compute_ground_state

# ============================================================
# Parameters
# ============================================================

mp.mp.dps = 80

c = 13
N = 8
T = 60
dps = 80

L = mp.log(c)
beta = L / (4 * mp.pi)

# ============================================================
# Build ground state
# ============================================================

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=dps,
)

lam_min, eigvec = compute_ground_state(Q)

# Repository/extract_zeros normalisation
coefficients = [
    mp.mpf(eigvec[i, 0])
    for i in range(eigvec.rows)
]

norm = mp.sqrt(sum(x*x for x in coefficients))
coefficients = [x / norm for x in coefficients]

# Canonical real-even coefficients:
#
#   v_0 = c_0
#   v_k = sqrt(2) c_k,  k >= 1
#
v_ground = [coefficients[N]]

for k in range(1, N + 1):
    v_ground.append(mp.sqrt(2) * coefficients[N + k])


# ============================================================
# Canonical basis response F_k(tau)
# ============================================================

def F_basis(k, tau):
    """
    F response to canonical basis vector e_k.
    """

    tau = mp.mpf(tau)
    exp_tL = mp.exp(-1j * tau * L)

    if k == 0:
        pairs = [(0, mp.mpf(1))]
    else:
        ck = 1 / mp.sqrt(2)
        pairs = [
            ( k, ck),
            (-k, ck),
        ]

    total = mp.mpc(0)

    for kk, ck in pairs:

        denom = 2 * mp.pi * kk / L - tau

        if denom == 0:
            term = mp.mpc(L)

        elif abs(denom * L) < mp.sqrt(mp.eps):
            term = (
                mp.expm1(1j * denom * L)
                / (1j * denom)
            )

        else:
            term = (
                exp_tL - 1
            ) / (1j * denom)

        total += ck * term

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * total
        / mp.sqrt(L)
    )


# ============================================================
# Canonical basis response F'(tau)
# ============================================================

def Fprime_basis(k, tau):
    """
    Analytic derivative of F_basis(k, tau).
    """

    tau = mp.mpf(tau)
    exp_tL = mp.exp(-1j * tau * L)

    if k == 0:
        pairs = [(0, mp.mpf(1))]
    else:
        ck = 1 / mp.sqrt(2)
        pairs = [
            ( k, ck),
            (-k, ck),
        ]

    H = mp.mpc(0)
    Hp = mp.mpc(0)

    for kk, ck in pairs:

        a = 2 * mp.pi * kk / L
        denom = a - tau

        # All zeta-zero ordinates used here are safely away
        # from the removable singularities.
        g = (
            exp_tL - 1
        ) / (1j * denom)

        gp = (
            -L * exp_tL / denom
            - 1j * (exp_tL - 1) / denom**2
        )

        H += ck * g
        Hp += ck * gp

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * (
            1j * L / 2 * H
            + Hp
        )
        / mp.sqrt(L)
    )


# ============================================================
# Pole functional
# ============================================================

def pole_basis(k):
    """
    P(e_k) for the canonical basis.
    """

    if k == 0:
        return 1 / beta**2

    return mp.sqrt(2) / (k**2 + beta**2)


# ============================================================
# Build the basic linear-functional rows
# ============================================================

def zero_row(j):
    """
    Linear functional

        v -> F(gamma_j)

    expressed in canonical coordinates.
    """

    gamma = mp.im(mp.zetazero(j))

    return mp.matrix([
        F_basis(k, gamma)
        for k in range(N + 1)
    ])


pole_row = mp.matrix([
    pole_basis(k)
    for k in range(N + 1)
])

derivative_row = mp.matrix([
    Fprime_basis(k, mp.im(mp.zetazero(1)))
    for k in range(N + 1)
])


# ============================================================
# General constrained problem
#
# Maximise |F'(gamma_1)| subject to
#
#   P(v) = 0
#   F(gamma_j) = 0, j=1,...,m
#   ||v|| = 1
#
# The maximum is
#
#   D_m = || projection of d onto ker(C_m) ||
# ============================================================

def constrained_D(m):
    """
    Return the maximum possible |F'(gamma_1)| under

        P(v) = 0
        F(gamma_j) = 0, j=1,...,m
        ||v|| = 1

    together with the constraint matrix and projected derivative.
    """

    rows = [pole_row]

    for j in range(1, m + 1):
        rows.append(zero_row(j))

    nr = len(rows)

    C = mp.matrix(nr, N + 1)

    for i in range(nr):
        for k in range(N + 1):
            C[i, k] = rows[i][k]

    G = C * C.T

    Cd = C * derivative_row

    y = mp.lu_solve(G, Cd)

    d_perp = derivative_row - C.T * y

    D = mp.sqrt(mp.fdot(d_perp, d_perp))

    return D, C, G, d_perp


# ============================================================
# Basic sanity checks
# ============================================================

gamma = [
    mp.im(mp.zetazero(j))
    for j in range(1, N + 1)
]

print("============================================================")
print("GROUND STATE")
print("============================================================")

print("lambda =", mp.nstr(lam_min, 40))

print("\nCanonical ground-state norm =")
print(mp.nstr(
    mp.sqrt(mp.fdot(mp.matrix(v_ground),
                    mp.matrix(v_ground))),
    40
))

print("\nFirst few zeta zeros:")
for j in range(1, min(5, N) + 1):
    print(
        j,
        mp.nstr(gamma[j-1], 35)
    )


# ============================================================
# Verify the ground state against the linear functionals
# ============================================================

print("\n============================================================")
print("GROUND-STATE FUNCTIONAL CHECKS")
print("============================================================")

for j in range(1, min(3, N) + 1):

    row = zero_row(j)

    value = mp.fdot(
        row,
        mp.matrix(v_ground)
    )

    print(
        f"F(gamma{j}) =",
        mp.nstr(value, 30)
    )

P_ground = mp.fdot(
    pole_row,
    mp.matrix(v_ground)
)

print(
    "P(v)       =",
    mp.nstr(P_ground, 30)
)

d_ground = mp.fdot(
    derivative_row,
    mp.matrix(v_ground)
)

print(
    "F'(gamma1) =",
    mp.nstr(d_ground, 30)
)


# ============================================================
# D_m sequence
# ============================================================

print("\n============================================================")
print("CONSTRAINED DERIVATIVE BOUND")
print("============================================================")
print()
print("m    constraints       D_m                         det(G_m)")
print("------------------------------------------------------------")

results = {}

for m in range(0, N):

    Dm, Cm, Gm, dpm = constrained_D(m)

    results[m] = (Dm, Cm, Gm, dpm)

    print(
        f"{m:<4}"
        f"{m+1:<18}"
        f"{mp.nstr(Dm, 30):<32}"
        f"{mp.nstr(mp.det(Gm), 20)}"
    )


# ============================================================
# Incremental geometry
#
# For each new zero constraint, determine how much of the
# currently available derivative direction it removes.
# ============================================================

print("\n============================================================")
print("INCREMENTAL CONSTRAINT GEOMETRY")
print("============================================================")

for m in range(1, N):

    D_prev, C_prev, G_prev, d_prev = results[m-1]

    # New constraint row = F(gamma_m)
    new_row = zero_row(m)

    # Project the new row into ker(C_prev)
    G_prev_inv_rhs = mp.lu_solve(
        G_prev,
        C_prev * new_row
    )

    new_perp = (
        new_row
        - C_prev.T * G_prev_inv_rhs
    )

    new_perp_norm = mp.sqrt(
        mp.fdot(new_perp, new_perp)
    )

    # Angle between the new constraint direction and
    # the currently available derivative direction.
    dot = mp.fdot(
        new_perp,
        d_prev
    )

    cos_theta = (
        dot
        / (new_perp_norm * D_prev)
    )

    print(f"\nAdding F(gamma_{m}) = 0")

    print(
        "||new constraint_perp|| =",
        mp.nstr(new_perp_norm, 40)
    )

    print(
        "||d_perp(previous)||    =",
        mp.nstr(D_prev, 40)
    )

    print(
        "cos(theta)              =",
        mp.nstr(cos_theta, 40)
    )

    print(
        "|cos(theta)|^2          =",
        mp.nstr(cos_theta**2, 40)
    )

    print(
        "D_previous              =",
        mp.nstr(D_prev, 40)
    )

    print(
        "D_new                   =",
        mp.nstr(results[m][0], 40)
    )

    print(
        "D_new / D_previous     =",
        mp.nstr(results[m][0] / D_prev, 40)
    )


# ============================================================
# Explicit maximiser for a selected m
# ============================================================

m_test = min(2, N - 1)

D_test, C_test, G_test, dperp_test = results[m_test]

v_star = dperp_test / D_test

print("\n============================================================")
print(f"MAXIMISING VECTOR FOR m = {m_test}")
print("============================================================")

print("D_m =")
print(mp.nstr(D_test, 60))

print("\n||v_star|| =")
print(mp.nstr(
    mp.sqrt(mp.fdot(v_star, v_star)),
    50
))

print("\nConstraint residuals:")

print(
    "P(v) =",
    mp.nstr(
        mp.fdot(pole_row, v_star),
        40
    )
)

for j in range(1, m_test + 1):

    row = zero_row(j)

    print(
        f"F(gamma{j}) =",
        mp.nstr(
            mp.fdot(row, v_star),
            40
        )
    )

print("\nDerivative:")

print(
    "F'(gamma1) =",
    mp.nstr(
        mp.fdot(derivative_row, v_star),
        50
    )
)


# ============================================================
# Local behaviour around the first few constrained zeros
# ============================================================

def F_vector(v, tau):
    return sum(
        v[k] * F_basis(k, tau)
        for k in range(N + 1)
    )


print("\n============================================================")
print("LOCAL F(tau) FOR v_star")
print("============================================================")

for j in range(1, m_test + 1):

    gj = gamma[j-1]

    print(f"\nAround gamma_{j} = {mp.nstr(gj, 30)}")

    for delta in ["-0.01", "0", "0.01"]:

        x = gj + mp.mpf(delta)

        print(
            delta,
            mp.nstr(
                F_vector(v_star, x),
                30
            )
        )

# ------------------------------------------------------------
# Build the contribution of each canonical v_k to F(tau)
#
# We construct a repository-style coefficient vector for
# the pure canonical basis vector e_k:
#
#   k=0: c_0 = 1
#   k>0: c_{+k}=c_{-k}=1/sqrt(2)
# ------------------------------------------------------------

def F_from_canonical_basis(k, tau):
    tau = mp.mpf(tau)

    total = mp.mpc(0)
    exp_tL = mp.exp(-1j * tau * L)

    if k == 0:
        coeff_pairs = [(0, mp.mpf(1))]
    else:
        ck = 1 / mp.sqrt(2)
        coeff_pairs = [(k, ck), (-k, ck)]

    for kk, ck in coeff_pairs:
        denom = 2 * mp.pi * kk / L - tau

        if denom == 0:
            term = mp.mpc(L)
        elif abs(denom * L) < mp.sqrt(mp.eps):
            term = mp.expm1(1j * denom * L) / (1j * denom)
        else:
            term = (exp_tL - 1) / (1j * denom)

        total += ck * term

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * total
        / mp.sqrt(L)
    )

def F_from_canonical_vector(v, tau):
    total = mp.mpf('0')

    for k in range(N + 1):
        total += v[k] * F_from_canonical_basis(k, tau)

    return total


# ============================================================
# Derivative response of each canonical basis vector
# ============================================================
def Fprime_from_canonical_basis(k, tau):
    tau = mp.mpf(tau)

    H = mp.mpc(0)
    Hp = mp.mpc(0)

    exp_tL = mp.exp(-1j * tau * L)

    if k == 0:
        coeff_pairs = [(0, mp.mpf(1))]
    else:
        ck = 1 / mp.sqrt(2)
        coeff_pairs = [(k, ck), (-k, ck)]

    for kk, ck in coeff_pairs:
        a = 2 * mp.pi * kk / L
        denom = a - tau

        # For the zeta zeros used here, we are safely away
        # from the removable singularities.
        g = (exp_tL - 1) / (1j * denom)

        gp = (
            -L * exp_tL / denom
            - 1j * (exp_tL - 1) / denom**2
        )

        H += ck * g
        Hp += ck * gp

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * (1j * L / 2 * H + Hp)
        / mp.sqrt(L)
    )

def T_canonical(v, t):
    """
    Trigonometric polynomial T_v(t) in the full symmetric
    coefficient representation.
    """

    total = mp.mpc(0)

    # m = 0
    total += v[0]

    # positive / negative pairs
    for k in range(1, N + 1):
        uk = v[k] / mp.sqrt(2)

        total += (
            uk * mp.exp(2j * mp.pi * k * t)
            + uk * mp.exp(-2j * mp.pi * k * t)
        )

    return total


def K_canonical(v, omega):
    """
    Volterra sine-chord kernel

        K_v(omega)
          = 2 int_0^omega T_v(t) T_v(omega-t) dt.
    """

    omega = mp.mpf(omega)

    if omega <= 0:
        return mp.mpf(0)

    if omega >= 1:
        # The Fourier-weight definition only needs
        # omega in [0,1] for the prime terms.
        raise ValueError("K_canonical expects 0 <= omega <= 1")

    integrand = lambda t: (
        T_canonical(v, t)
        * T_canonical(v, omega - t)
    )

    return 2 * mp.quad(integrand, [0, omega])


def ghat(v, xi):
    """
    Fourier weight

        ghat(xi) = pi K(1 - |xi|/Delta)

    for |xi| <= Delta.
    """

    xi = mp.mpf(xi)

    Delta = L / (2 * mp.pi)

    if abs(xi) > Delta:
        return mp.mpf(0)

    omega = 1 - abs(xi) / Delta

    return mp.pi * K_canonical(v, omega)

def G_complex_from_canonical_basis(k, z):
    z = mp.mpc(z)

    total = mp.mpc(0)

    exp_zL = mp.exp(-1j * z * L)

    if k == 0:
        coeff_pairs = [(0, mp.mpf(1))]
    else:
        ck = 1 / mp.sqrt(2)
        coeff_pairs = [(k, ck), (-k, ck)]

    for kk, ck in coeff_pairs:

        a = 2 * mp.pi * kk / L
        denom = a - z

        if denom == 0:
            term = mp.mpc(L)

        elif abs(denom * L) < mp.sqrt(mp.eps):
            term = mp.expm1(1j * denom * L) / (1j * denom)

        else:
            term = (exp_zL - 1) / (1j * denom)

        total += ck * term

    return (
        mp.exp(1j * z * L / 2)
        * total
        / mp.sqrt(L)
    )


def G_complex(v, z):
    total = mp.mpc(0)

    for k in range(N + 1):
        total += (
            v[k]
            * G_complex_from_canonical_basis(k, z)
        )

    return total



# ============================================================
# CELL_5 — ARCHIMEDEAN FUNCTIONAL AUDIT
#
# Goal:
#
# Independently audit the archimedean contribution by comparing
#
#   A) the repository finite-T Galerkin matrix, with its prime
#      and pole pieces removed,
#
# against
#
# B) the explicit Archimedean quadratic form
#
#     (1 / 2*pi^2) int_{-T}^{T}
#         h_+(r) D_v(r) dr
#
# which, by evenness, is
#
#     (1 / pi^2) int_0^T
#         h_+(r) D_v(r) dr.
#
# We also compare individual canonical matrix elements.
#
# DO NOT modify Cell_1 ... Cell_4.
# ============================================================

print("\n" + "=" * 60)
print("CELL_5 — ARCHIMEDEAN FUNCTIONAL AUDIT")
print("=" * 60)

mp.mp.dps = 80

print("\nParameters:")
print("N =", N)
print("c =", c)
print("L =", mp.nstr(L, 50))


# ============================================================
# 1. FULL-SPACE / CANONICAL EMBEDDING
# ============================================================

u_star = mp.matrix(2 * N + 1, 1)

for m in range(-N, N + 1):

    if m == 0:
        u_star[m + N] = v_star[0]
    else:
        u_star[m + N] = v_star[abs(m)] / mp.sqrt(2)


print("\nFull-space dimension:")
print("u_star =", u_star.rows, "x", u_star.cols)

print("\n||u_star|| =")
print(mp.nstr(
    mp.sqrt(mp.fdot(u_star, u_star)),
    50
))


# ============================================================
# 2. PRIME-POWER LIST
# ============================================================

def prime_power_terms_cell5(c):
    """
    Return (q, Lambda(q)) for all prime powers q <= c.
    """

    c_int = int(mp.floor(c))

    terms = []

    for p0 in range(2, c_int + 1):

        is_prime = True

        for d0 in range(2, int(mp.sqrt(p0)) + 1):

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
                    mp.log(p0)
                )
            )

            q0 *= p0

    return terms


prime_terms_5 = prime_power_terms_cell5(c)


# ============================================================
# 3. PRIME MATRIX
#
# Construct the complete full-space Q_prime directly from
# divided differences.
# ============================================================

def Q_prime_power_cell5(q, Lambda_q):

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
            * mp.sin(2 * mp.pi * x * a)
        )

    def psi_derivative(x):

        return (
            prefactor
            * 2 * mp.pi * a
            * mp.cos(2 * mp.pi * x * a)
        )

    for i, m in enumerate(range(-N, N + 1)):

        for j, n in enumerate(range(-N, N + 1)):

            if m != n:

                Qq[i, j] = (
                    psi(m) - psi(n)
                ) / mp.mpf(m - n)

            else:

                Qq[i, j] = psi_derivative(m)

    return Qq


Q_prime_5 = mp.matrix(
    2 * N + 1,
    2 * N + 1
)

for q, Lambda_q in prime_terms_5:

    Q_prime_5 += Q_prime_power_cell5(
        q,
        Lambda_q
    )


# ============================================================
# 4. POLE MATRIX
#
# psi_0(x) =
#
#   1/pi int_0^L
#       2 cosh(y/2)
#       sin(2*pi*x*(1-y/L))
#   dy
#
# Construct Q_psi0 directly from divided differences.
#
# ============================================================

def psi_pole_cell5(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * mp.sin(
            2 * mp.pi
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


def psi_pole_derivative_cell5(x):

    x = mp.mpf(x)

    integrand = lambda y: (
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
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L]
        )
    )


size = 2 * N + 1

Q_pole_5 = mp.matrix(size, size)

for i, m in enumerate(range(-N, N + 1)):

    for j, n in enumerate(range(-N, N + 1)):

        if m != n:

            Q_pole_5[i, j] = (
                psi_pole_cell5(m)
                - psi_pole_cell5(n)
            ) / mp.mpf(m - n)

        else:

            Q_pole_5[i, j] = (
                psi_pole_derivative_cell5(m)
            )


# ============================================================
# 5. BASIC POLE SANITY CHECK
#
# Compare the pole matrix quadratic form against 2*g(i/2).
# ============================================================

pole_matrix_form = mp.fdot(
    u_star,
    Q_pole_5 * u_star
)

pole_explicit = mp.re(
    2 * G_complex(
        v_star,
        1j / 2
    )
)

print("\n" + "-" * 60)
print("1. POLE MATRIX SANITY CHECK")
print("-" * 60)

print("\n<v_star, Q_pole v_star> =")
print(mp.nstr(
    pole_matrix_form,
    60
))

print("\n2 g(i/2) =")
print(mp.nstr(
    pole_explicit,
    60
))

print("\nDifference =")
print(mp.nstr(
    pole_matrix_form - pole_explicit,
    60
))


# ============================================================
# 6. EXPLICIT ARCHIMEDEAN FUNCTIONS
# ============================================================

def h_plus_cell5(r):
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
# 7. SOURCE-LEVEL ARCHIMEDEAN QUADRATIC INTEGRAND
#
# The repository defines
#
#   psi_arch,T(x)
#
#       = 1/(2*pi^2) int_{-T}^T
#           h_+(r) S(r,x) dr
#
# and the Galerkin matrix is formed by divided differences:
#
#   Q_arch(m,n)
#       = D[psi_arch](m,n).
#
# Therefore, by linearity,
#
#   <u,Q_arch u>
#
#       = 1/(2*pi^2) int_{-T}^T
#           h_+(r) D_v(r) dr
#
# where
#
#   D_v(r)
#       = sum_{m,n} u_m u_n D_S(m,n;r).
#
# Since D_v(r) and h_+(r) are even,
#
#   <u,Q_arch u>
#
#       = 1/pi^2 int_0^T
#           h_+(r) D_v(r) dr.
#
# This is the quantity that must be compared with the
# repository Archimedean quadratic form.
# ============================================================

def S_arch_cell5(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    integrand = lambda y: (
        mp.sin(
            2 * mp.pi
            * x
            * (1 - y / L)
        )
        * mp.cos(r * y)
    )

    return mp.quad(
        integrand,
        [0, L]
    )


def dS_direct_cell5(r, x):
    r = mp.mpf(r)
    x = mp.mpf(x)

    integrand = lambda y: (
        2 * mp.pi
        * (1 - y / L)
        * mp.cos(
            2 * mp.pi
            * x
            * (1 - y / L)
        )
        * mp.cos(r * y)
    )

    return mp.quad(
        integrand,
        [0, L]
    )


def source_divdiff_cell5(r, m, n):
    r = mp.mpf(r)

    if m != n:
        return (
            S_arch_cell5(r, m)
            - S_arch_cell5(r, n)
        ) / mp.mpf(m - n)

    return dS_direct_cell5(r, m)


def source_quadratic_integrand_cell5(r):
    r = mp.mpf(r)

    total = mp.mpf(0)

    for i, m in enumerate(
        range(-N, N + 1)
    ):
        for j, n in enumerate(
            range(-N, N + 1)
        ):
            total += (
                u_star[i]
                * u_star[j]
                * source_divdiff_cell5(
                    r,
                    m,
                    n,
                )
            )

    return total


def arch_explicit_cell5(T):
    T = mp.mpf(T)

    def integrand(r):
        return (
            h_plus_cell5(r)
            * source_quadratic_integrand_cell5(r)
        )

    step = mp.mpf(5)

    points = [mp.mpf(0)]

    r0 = step

    while r0 < T:
        points.append(r0)
        r0 += step

    points.append(T)

    integral = mp.quad(
        integrand,
        points
    )

    return integral / (mp.pi ** 2)


# ============================================================
# 8. REPOSITORY ARCHIMEDEAN MATRIX
#
# For each T:
#
#   Q_arch,T(repo)
#       = Q_T(repo) - Q_prime - Q_pole
#
# Then evaluate its quadratic form.
#
# ============================================================

def repository_arch_matrix(T):

    Q_T = build_galerkin_matrix(
        c=c,
        N=N,
        T=T,
        dps=80,
    )

    Q_arch = (
        Q_T
        - Q_prime_5
        - Q_pole_5
    )

    return Q_T, Q_arch


# ============================================================
# 9. FINITE-T COMPARISON
# ============================================================

print("\n" + "-" * 60)
print("2. ARCHIMEDEAN FINITE-T COMPARISON")
print("-" * 60)

print()
print(
    "T".ljust(8),
    "Repository A_T".ljust(34),
    "Explicit A_T".ljust(34),
    "Difference"
)

print("-" * 115)

arch_results_5 = {}

for T_arch in [20, 40, 80, 120, 200]:

    Q_T, Q_arch = repository_arch_matrix(
        T_arch
    )

    repo_A = mp.fdot(
        u_star,
        Q_arch * u_star
    )

    explicit_A = arch_explicit_cell5(
        T_arch
    )

    diff = repo_A - explicit_A

    arch_results_5[T_arch] = (
        repo_A,
        explicit_A,
        diff,
        Q_T,
        Q_arch
    )

    print(
        f"{T_arch:<8}",
        mp.nstr(repo_A, 28).ljust(34),
        mp.nstr(explicit_A, 28).ljust(34),
        mp.nstr(diff, 28)
    )


# ============================================================
# 10. INDEPENDENT ARCHIMEDEAN SOURCE
#
# We now evaluate the source function
#
#   psi_R,T(x)
#
# directly from the definition:
#
#   1/(2*pi^2) int_{-T}^T
#       h_+(r) S(r,x,L) dr
#
# with
#
#   S(r,x,L)
#      = int_0^L
#          sin(2*pi*x*(1-y/L))
#          cos(r y) dy.
#
# This gives an independent check of the source-level
# construction before the divided-difference matrix is formed.
# ============================================================

def psi_arch_T_cell5(x, T):

    x = mp.mpf(x)
    T = mp.mpf(T)

    integrand = lambda r: (
        h_plus_cell5(r)
        * S_arch_cell5(r, x)
    )

    # The integrand is even in r.
    #
    # Therefore
    #
    # 1/(2*pi^2) int_{-T}^T
    #     = 1/pi^2 int_0^T.
    #
    step = mp.mpf(5)

    points = [mp.mpf(0)]

    r0 = step

    while r0 < T:

        points.append(r0)

        r0 += step

    points.append(T)

    integral = mp.quad(
        integrand,
        points
    )

    return integral / (mp.pi ** 2)


# ============================================================
# 11. SOURCE-LEVEL TEST
#
# Compare psi_R,T(x) with the corresponding matrix source
# implied by the repository archimedean matrix only indirectly
# through the quadratic form.
#
# We evaluate several integer x values.
# ============================================================

print("\n" + "-" * 60)
print("3. ARCHIMEDEAN SOURCE-LEVEL SAMPLES")
print("-" * 60)

for T_test in [40, 120]:

    print(
        "\nT =",
        T_test
    )

    for x_test in [0, 1, 2, 3]:

        value = psi_arch_T_cell5(
            x_test,
            T_test
        )

        print(
            "x =",
            x_test,
            " psi_arch,T(x) =",
            mp.nstr(value, 40)
        )


# ============================================================
# 12. REQUIRED ARCHIMEDEAN VALUE
#
# Use the CURRENT Cell-5 Galerkin matrix and CURRENT v_star.
#
# Q_current = <u_star, Q u_star>
#
# Therefore
#
#   Q_arch,required =
#       Q_current - Q_prime - Q_pole.
# ============================================================

Q_current_form = mp.fdot(
    u_star,
    Q * u_star
)

prime_current_form = mp.fdot(
    u_star,
    Q_prime_5 * u_star
)

required_arch = (
    Q_current_form
    - prime_current_form
    - pole_matrix_form
)

print("\n" + "-" * 60)
print("4. REQUIRED ARCHIMEDEAN VALUE")
print("-" * 60)

print("\nCurrent Galerkin Q =")
print(mp.nstr(Q_current_form, 60))

print("\nPrime matrix form =")
print(mp.nstr(prime_current_form, 60))

print("\nPole matrix form =")
print(mp.nstr(pole_matrix_form, 60))

print("\nRequired Archimedean =")
print(mp.nstr(required_arch, 60))

# ============================================================
# 13. COMPARE REQUIRED VALUE WITH FINITE-T RESULTS
# ============================================================

print("\n" + "-" * 60)
print("5. REQUIRED VALUE VS FINITE-T")
print("-" * 60)

for T_test in [20, 40, 80, 120, 200]:

    repo_A, explicit_A, diff, _, _ = (
        arch_results_5[T_test]
    )

    print(
        "\nT =",
        T_test
    )

    print(
        "  Repository A_T =",
        mp.nstr(repo_A, 40)
    )

    print(
        "  Explicit A_T   =",
        mp.nstr(explicit_A, 40)
    )

    print(
        "  Required A     =",
        mp.nstr(required_arch, 40)
    )

    print(
        "  Explicit - Required =",
        mp.nstr(
            explicit_A - required_arch,
            40
        )
    )


# ============================================================
# 14. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("CELL_5 SUMMARY")
print("=" * 60)

print("\nRequired archimedean contribution:")
print(
    mp.nstr(
        required_arch,
        60
    )
)

print("\nExplicit A_200:")
print(
    mp.nstr(
        arch_results_5[200][1],
        60
    )
)

print("\nRepository A_200:")
print(
    mp.nstr(
        arch_results_5[200][0],
        60
    )
)

print("\nRepository - explicit at T=200:")
print(
    mp.nstr(
        arch_results_5[200][2],
        60
    )
)

print("\n" + "=" * 60)
print("END CELL_5")
print("=" * 60)
