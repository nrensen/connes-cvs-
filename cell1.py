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
# DIAGNOSTICS FOR THE m = 2 MAXIMISING VECTOR
# ============================================================

print("\n" + "=" * 60)
print("DIAGNOSTICS FOR v_star (m = 2)")
print("=" * 60)

# ------------------------------------------------------------
# 1. Canonical coefficients of v_star
# ------------------------------------------------------------

print("\nCanonical coefficients v_star:")
for k in range(N + 1):
    print(
        f"v[{k}] = {mp.nstr(v_star[k], 40)}"
    )

# ------------------------------------------------------------
# 2. F'(gamma_j) for the first five zeta zeros
# ------------------------------------------------------------

print("\nF'(gamma_j):")

for j in range(1, 6):
    gamma = mp.im(mp.zetazero(j))

    Fp = mp.mpf("0")

    for k in range(N + 1):
        Fp += (
            v_star[k]
            * Fprime_from_canonical_basis(k, gamma)
        )

    print(
        f"j={j}  gamma={mp.nstr(gamma, 30)}"
    )
    print(
        f"     F'(gamma_j) = {mp.nstr(Fp, 40)}"
    )

# ------------------------------------------------------------
# 3. F(gamma_j) for the first five zeta zeros
# ------------------------------------------------------------

print("\nF(gamma_j):")

for j in range(1, 6):
    gamma = mp.im(mp.zetazero(j))

    Fj = F_from_canonical_vector(v_star, gamma)

    print(
        f"j={j}  gamma={mp.nstr(gamma, 30)}"
    )
    print(
        f"     F(gamma_j) = {mp.nstr(Fj, 40)}"
    )

# ------------------------------------------------------------
# 4. Pole functional P(v_star)
# ------------------------------------------------------------

beta = L / (4 * mp.pi)

P_star = (
    v_star[0] / beta**2
    + mp.sqrt(2) * sum(
        v_star[k] / (k**2 + beta**2)
        for k in range(1, N + 1)
    )
)

print("\nPole functional:")
print(
    "P(v_star) =",
    mp.nstr(P_star, 50)
)

# ------------------------------------------------------------
# 5. Norm of v_star
# ------------------------------------------------------------

norm_star = mp.sqrt(mp.fdot(v_star, v_star))

print("\nNorm:")
print(
    "||v_star|| =",
    mp.nstr(norm_star, 50)
)

# ------------------------------------------------------------
# 6. First unconstrained zero: F(gamma_3)
# ------------------------------------------------------------

gamma3 = mp.im(mp.zetazero(3))

F_gamma3 = F_from_canonical_vector(v_star, gamma3)

print("\nFirst unconstrained zero:")
print(
    "gamma_3 =",
    mp.nstr(gamma3, 40)
)
print(
    "F(gamma_3) =",
    mp.nstr(F_gamma3, 50)
)


