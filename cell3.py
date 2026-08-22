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

# ============================================================
# CELL_3 — EXPLICIT WEIL-SIDE CHECK
#
# Independent evaluation of
#
#   <v,Q_inf v>
#
# from
#
#   prime + pole + archimedean
#
# rather than from the zeros.
# ============================================================

import mpmath as mp

print("\n" + "=" * 60)
print("CELL_3 — EXPLICIT WEIL-SIDE CHECK")
print("=" * 60)

mp.mp.dps = 80

# ------------------------------------------------------------
# Parameters / sanity checks
# ------------------------------------------------------------

print("\nParameters:")
print("N =", N)
print("c =", c)
print("L =", mp.nstr(L, 40))

norm_v = mp.sqrt(mp.fdot(v_star, v_star))

print("\n||v_star|| =")
print(mp.nstr(norm_v, 50))


# ============================================================
# 1. COMPLEX TEST FUNCTION g_v(z)
#
# For real z, this agrees with F_from_canonical_vector(v,z).
#
# Unlike F_from_canonical_basis(), this version does NOT take
# mp.re(...), because we need g(i/2) for the pole term.
# ============================================================

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


# ------------------------------------------------------------
# Sanity check:
# complex evaluator should reproduce the real evaluator.
# ------------------------------------------------------------

print("\nTest function consistency:")

for j in range(1, 4):
    gamma = mp.im(mp.zetazero(j))

    g1 = F_from_canonical_vector(v_star, gamma)
    g2 = G_complex(v_star, gamma)

    print(
        f"gamma_{j}:"
    )
    print(
        "  F_real    =",
        mp.nstr(g1, 40)
    )
    print(
        "  G_complex =",
        mp.nstr(g2, 40)
    )
    print(
        "  difference =",
        mp.nstr(g1 - mp.re(g2), 20)
    )


# ============================================================
# 2. GALERKIN QUADRATIC FORM
# ============================================================

# Q is the repository full-space matrix on indices
#
#       -N, ..., 0, ..., N
#
# whereas v_star is the canonical even-sector vector
#
#       (v_0, ..., v_N).
#
# Embed v_star isometrically into the full symmetric
# coefficient vector:
#
#       u_0 = v_0
#       u_{+k} = u_{-k} = v_k / sqrt(2).

u_star = mp.matrix(2 * N + 1, 1)

u_star[N] = v_star[0]

for k in range(1, N + 1):
    uk = v_star[k] / mp.sqrt(2)

    u_star[N + k] = uk
    u_star[N - k] = uk

print("\nFull-space dimension check:")
print("Q dimensions      =", Q.rows, "x", Q.cols)
print("u_star dimension  =", u_star.rows)

assert Q.rows == 2 * N + 1
assert Q.cols == 2 * N + 1
assert u_star.rows == 2 * N + 1

Qv = Q * u_star

Q_form = mp.fdot(
    u_star,
    Qv
)

print("\n" + "-" * 60)
print("1. GALERKIN QUADRATIC FORM")
print("-" * 60)

print(
    "<v_star, Q_inf v_star> ="
)
print(mp.nstr(Q_form, 60))

# ============================================================
# 3. PRIME CONTRIBUTION
#
#   -1/pi sum_{q=p^a <= c}
#       Lambda(q)/sqrt(q) *
#       ghat(log(q)/(2*pi))
#
# We obtain ghat directly from the canonical Fourier weight.
#
# For the symmetric coefficients:
#
#   T(t) = sum u_m exp(2*pi*i*m*t)
#
# and
#
#   K(w) = 2 int_0^w T(t) T(w-t) dt.
#
# We only need K at
#
#   w = 1 - log(q)/L.
#
# ============================================================

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


# ------------------------------------------------------------
# Prime powers q <= c
# ------------------------------------------------------------

def prime_power_terms(c):
    """
    Return (q, Lambda(q)) for all prime powers q <= c.
    """

    c_int = int(mp.floor(c))

    terms = []

    # Simple primality test is sufficient for this tiny c.
    for p0 in range(2, c_int + 1):

        is_prime = True

        for d0 in range(2, int(mp.sqrt(p0)) + 1):
            if p0 % d0 == 0:
                is_prime = False
                break

        if not is_prime:
            continue

        q = p0

        while q <= c_int:

            terms.append(
                (
                    mp.mpf(q),
                    mp.log(p0)
                )
            )

            q *= p0

    return terms


prime_terms = prime_power_terms(c)

print("\n" + "-" * 60)
print("2. PRIME CONTRIBUTION")
print("-" * 60)

print("\nPrime powers included:")

for q, Lambda_q in prime_terms:
    print(
        "q =",
        mp.nstr(q, 10),
        " Lambda(q) =",
        mp.nstr(Lambda_q, 20)
    )


prime_sum = mp.mpf(0)

print("\nIndividual prime-power terms:")

for q, Lambda_q in prime_terms:

    xi = mp.log(q) / (2 * mp.pi)

    gh = ghat(v_star, xi)

    term = (
        -1 / mp.pi
        * Lambda_q
        / mp.sqrt(q)
        * gh
    )

    prime_sum += term

    print(
        "q =",
        mp.nstr(q, 8),
        " xi =",
        mp.nstr(xi, 25),
        " ghat =",
        mp.nstr(gh, 30),
        " contribution =",
        mp.nstr(term, 30)
    )


print("\nPrime contribution =")
print(mp.nstr(prime_sum, 60))


# ============================================================
# 4. POLE CONTRIBUTION
#
#       2 g(i/2)
# ============================================================

print("\n" + "-" * 60)
print("3. POLE CONTRIBUTION")
print("-" * 60)

g_pole = G_complex(v_star, 1j / 2)

pole_sum = 2 * g_pole

print("g(i/2) =")
print(mp.nstr(g_pole, 60))

print("\n2 g(i/2) =")
print(mp.nstr(pole_sum, 60))


# ============================================================
# 5. ARCHIMEDEAN FUNCTION h_+(r)
#
#   h_+(r)
#      = Re psi_Gamma(1/4 + ir/2) - log(pi)
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


# ------------------------------------------------------------
# Sanity check h_+
# ------------------------------------------------------------

print("\nh_+(r) samples:")

for r in [0, 1, 5, 10, 20, 50]:
    print(
        "r =",
        r,
        " h_+ =",
        mp.nstr(h_plus(r), 30)
    )


# ============================================================
# 6. ARCHIMEDEAN CONTRIBUTION
#
#       1/(2*pi) int_R h_+(r) g(r) dr
#
# Since both h_+ and g are even:
#
#       1/pi int_0^T h_+(r) g(r) dr
#
# We initially use finite T and examine convergence.
# ============================================================

def arch_contribution(T):
    T = mp.mpf(T)

    integrand = lambda r: (
        h_plus(r)
        * mp.re(G_complex(v_star, r))
    )

    # Break the integral into pieces.  This is considerably
    # more stable for the oscillatory g(r).
    #
    # The exact breakpoints are not mathematically important;
    # they simply prevent mp.quad from attempting one enormous
    # oscillatory interval.

    step = mp.mpf(5)

    points = [mp.mpf(0)]

    x = step

    while x < T:
        points.append(x)
        x += step

    points.append(T)

    integral = mp.quad(
        integrand,
        points
    )

    return integral / mp.pi


print("\n" + "-" * 60)
print("4. ARCHIMEDEAN CONTRIBUTION")
print("-" * 60)

print("\nFinite-T convergence:")

arch_values = {}

for T_arch in [20, 40, 80, 120, 200, 400]:

    A = arch_contribution(T_arch)

    arch_values[T_arch] = A

    print(
        f"T = {T_arch:3d}   "
        f"A_T = {mp.nstr(A, 50)}"
    )


# Use the largest cutoff as the current numerical estimate.
arch_sum = arch_values[400]

print("\nCurrent archimedean estimate =")
print(mp.nstr(arch_sum, 60))


# ============================================================
# 7. COMPLETE EXPLICIT WEIL SUM
# ============================================================

weil_sum = (
    prime_sum
    + pole_sum
    + arch_sum
)

print("\n" + "-" * 60)
print("5. COMPLETE EXPLICIT WEIL FORM")
print("-" * 60)

print("Prime       =")
print(mp.nstr(prime_sum, 60))

print("\nPole        =")
print(mp.nstr(pole_sum, 60))

print("\nArchimedean =")
print(mp.nstr(arch_sum, 60))

print("\n--------------------------------")

print("\nWeil sum =")
print(mp.nstr(weil_sum, 60))


# ============================================================
# 8. COMPARISON WITH GALERKIN MATRIX
# ============================================================

difference = weil_sum - Q_form

relative_error = (
    abs(difference) / abs(Q_form)
)

print("\n" + "-" * 60)
print("6. DICTIONARY CHECK")
print("-" * 60)

print("\nGalerkin quadratic form =")
print(mp.nstr(Q_form, 60))

print("\nExplicit Weil form =")
print(mp.nstr(weil_sum, 60))

print("\nDifference =")
print(mp.nstr(difference, 60))

print("\nRelative difference =")
print(mp.nstr(relative_error, 30))


# ============================================================
# 9. FINAL CHECK
# ============================================================

print("\n" + "=" * 60)
print("END CELL_3")
print("=" * 60)


