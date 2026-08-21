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



# ============================================================
# CELL_4 — PRIME FUNCTIONAL AUDIT
#
# Independent audit of the prime contribution.
#
# We compare two constructions:
#
#   A) Direct divided-difference matrix:
#
#        Q_prime = Q_{psi_p}
#
#      using
#
#        psi_p(x)
#          = -1/pi sum_q Lambda(q)/sqrt(q)
#              sin(2*pi*x*(1-log(q)/L))
#
#   B) Explicit Guinand-Weil Fourier-side expression:
#
#        -1/pi sum_q Lambda(q)/sqrt(q)
#              ghat(log(q)/(2*pi))
#
# The two must agree if the finite dictionary is
# implemented consistently.
# ============================================================

print("\n" + "=" * 60)
print("CELL_4 — PRIME FUNCTIONAL AUDIT")
print("=" * 60)

mp.mp.dps = 80

print("\nParameters:")
print("N =", N)
print("c =", c)
print("L =", mp.nstr(L, 50))

# ------------------------------------------------------------
# Full-space coefficients u corresponding to v_star
#
# u_0 = v_0
# u_{+k} = u_{-k} = v_k / sqrt(2)
# ------------------------------------------------------------

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
# 1. PRIME-POWER LIST
# ============================================================

def prime_power_terms_cell4(c):
    """
    Return ordered list of (q, Lambda(q)) for prime powers q <= c.
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


prime_terms_4 = prime_power_terms_cell4(c)

print("\n" + "-" * 60)
print("1. PRIME POWERS")
print("-" * 60)

for q, lam in prime_terms_4:
    print(
        "q =",
        mp.nstr(q, 12),
        " Lambda(q) =",
        mp.nstr(lam, 30)
    )


# ============================================================
# 2. DIRECT PRIME SOURCE
#
# psi_p(x) =
#
#   -1/pi sum_q Lambda(q)/sqrt(q)
#       sin(2*pi*x*(1-log(q)/L))
#
# ============================================================

def psi_prime_q(x, q, Lambda_q):
    """
    Contribution of one prime power q to psi_p(x).
    """

    x = mp.mpf(x)

    a = 1 - mp.log(q) / L

    return (
        -1 / mp.pi
        * Lambda_q / mp.sqrt(q)
        * mp.sin(2 * mp.pi * x * a)
    )


def psi_prime(x):
    """
    Complete truncated prime source psi_p^(c)(x).
    """

    total = mp.mpf(0)

    for q, Lambda_q in prime_terms_4:

        total += psi_prime_q(
            x,
            q,
            Lambda_q
        )

    return total


# ============================================================
# 3. DIVIDED-DIFFERENCE MATRIX FOR ONE PRIME POWER
#
# Q_psi[m,n] =
#
#   (psi(m)-psi(n))/(m-n),  m != n
#
#   psi'(m),               m = n
#
# ============================================================

def Q_prime_power(q, Lambda_q):
    """
    Full-space divided-difference matrix for a single
    prime-power contribution.
    """

    size = 2 * N + 1

    Qq = mp.matrix(size, size)

    a = 1 - mp.log(q) / L

    prefactor = (
        -1 / mp.pi
        * Lambda_q / mp.sqrt(q)
    )

    def psi(x):
        return prefactor * mp.sin(
            2 * mp.pi * x * a
        )

    def psi_prime_derivative(x):
        return prefactor * (
            2 * mp.pi * a
            * mp.cos(2 * mp.pi * x * a)
        )

    for i, m in enumerate(range(-N, N + 1)):

        for j, n in enumerate(range(-N, N + 1)):

            if m != n:

                Qq[i, j] = (
                    psi(m) - psi(n)
                ) / mp.mpf(m - n)

            else:

                Qq[i, j] = psi_prime_derivative(m)

    return Qq


# ============================================================
# 4. DIRECT QUADRATIC FORM, PRIME POWER BY PRIME POWER
# ============================================================

print("\n" + "-" * 60)
print("2. DIRECT DIVIDED-DIFFERENCE PRIME TERMS")
print("-" * 60)

direct_prime_sum = mp.mpf(0)

direct_prime_terms = {}

for q, Lambda_q in prime_terms_4:

    Qq = Q_prime_power(q, Lambda_q)

    value = mp.fdot(
        u_star,
        Qq * u_star
    )

    direct_prime_terms[q] = value

    direct_prime_sum += value

    print(
        "q =",
        mp.nstr(q, 10),
        " direct quadratic form =",
        mp.nstr(value, 50)
    )

print("\nDirect prime contribution =")
print(mp.nstr(direct_prime_sum, 60))


# ============================================================
# 5. FOURIER-SIDE PRIME TERMS
#
# This deliberately uses the existing ghat() from Cell_3.
# ============================================================

print("\n" + "-" * 60)
print("3. FOURIER-SIDE PRIME TERMS")
print("-" * 60)

fourier_prime_sum = mp.mpf(0)

fourier_prime_terms = {}

for q, Lambda_q in prime_terms_4:

    xi = mp.log(q) / (2 * mp.pi)

    gh = ghat(
        v_star,
        xi
    )

    term = (
        -1 / mp.pi
        * Lambda_q
        / mp.sqrt(q)
        * gh
    )

    fourier_prime_terms[q] = term

    fourier_prime_sum += mp.re(term)

    print(
        "q =",
        mp.nstr(q, 10)
    )

    print(
        "  xi =",
        mp.nstr(xi, 30)
    )

    print(
        "  ghat =",
        mp.nstr(gh, 50)
    )

    print(
        "  Fourier quadratic form =",
        mp.nstr(term, 50)
    )

print("\nFourier-side prime contribution =")
print(mp.nstr(fourier_prime_sum, 60))


# ============================================================
# 6. TERM-BY-TERM COMPARISON
# ============================================================

print("\n" + "-" * 60)
print("4. TERM-BY-TERM COMPARISON")
print("-" * 60)

print()
print(
    "q".ljust(8),
    "DIRECT".ljust(32),
    "FOURIER".ljust(32),
    "DIFFERENCE"
)

print("-" * 110)

for q, Lambda_q in prime_terms_4:

    d = direct_prime_terms[q]
    f = fourier_prime_terms[q]

    diff = d - f

    print(
        mp.nstr(q, 6).ljust(8),
        mp.nstr(d, 25).ljust(32),
        mp.nstr(f, 25).ljust(32),
        mp.nstr(diff, 25)
    )


# ============================================================
# 7. TOTAL COMPARISON
# ============================================================

prime_difference = (
    direct_prime_sum
    - fourier_prime_sum
)

print("\n" + "-" * 60)
print("5. PRIME DICTIONARY CHECK")
print("-" * 60)

print("\nDirect divided-difference prime form =")
print(mp.nstr(direct_prime_sum, 60))

print("\nFourier-side prime form =")
print(mp.nstr(fourier_prime_sum, 60))

print("\nDifference =")
print(mp.nstr(prime_difference, 60))

print("\nRelative difference =")
print(mp.nstr(
    abs(prime_difference)
    / max(abs(direct_prime_sum), mp.mpf("1e-100")),
    40
))


# ============================================================
# 8. COMPARE DIRECT PRIME FORM WITH THE PRIME PART OF Q
#
# Q is the complete cutoff-free matrix from Cell_1.
#
# We cannot extract the repository's prime matrix directly
# from Q, so this section only records the direct prime result
# for comparison with the explicit formula.
# ============================================================

print("\n" + "-" * 60)
print("6. SUMMARY")
print("-" * 60)

print("\nDirect prime contribution =")
print(mp.nstr(direct_prime_sum, 60))

print("\nFourier-side prime contribution =")
print(mp.nstr(fourier_prime_sum, 60))

print("\nPrime-side discrepancy =")
print(mp.nstr(prime_difference, 60))


# ============================================================
# 9. FINAL CHECK
# ============================================================

print("\n" + "=" * 60)
print("END CELL_4")
print("=" * 60)

