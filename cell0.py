import mpmath as mp
from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros

mp.mp.dps = 80

c = 13
N = 8
T = 60

Q = build_galerkin_matrix(c=c, N=N, T=T, dps=80)
lam_min, eigvec = compute_ground_state(Q)

L = mp.log(c)

# Same normalisation used by extract_zeros
coefficients = [
    mp.mpf(eigvec[i, 0]) for i in range(eigvec.rows)
]

norm = mp.sqrt(sum(v*v for v in coefficients))
coefficients = [v/norm for v in coefficients]


def F_even(tau):
    tau = mp.mpf(tau)
    total = mp.mpc(0)

    exp_tL = mp.exp(-1j * tau * L)

    for k in range(-N, N + 1):
        ck = coefficients[k + N]

        if ck == 0:
            continue

        denom = 2 * mp.pi * k / L - tau

        if denom == 0:
            term = mp.mpc(L)
        elif abs(denom * L) < mp.sqrt(mp.eps):
            term = mp.expm1(1j * denom * L) / (1j * denom)
        else:
            term = (exp_tL - 1) / (1j * denom)

        total += ck * term

    total /= mp.sqrt(L)

    return mp.re(mp.exp(1j * tau * L / 2) * total)


gamma1 = mp.im(mp.zetazero(1))

print("lambda =", mp.nstr(lam_min, 30))
print("gamma1 =", mp.nstr(gamma1, 30))
print("F(gamma1) =", mp.nstr(F_even(gamma1), 30))

z = extract_zeros(
    eigvec,
    c=c,
    n_zeros=1,
    dps=80
)

print("repository gamma_detected =",
      mp.nstr(z[0]["gamma_detected"], 30))

print("repository error =",
      mp.nstr(z[0]["error"], 10))

def F_even_deriv_numeric(tau):
    return mp.diff(F_even, tau)

def F_even_deriv_analytic(tau):
    tau = mp.mpf(tau)

    H = mp.mpc(0)
    Hp = mp.mpc(0)

    exp_tL = mp.exp(-1j * tau * L)

    for k in range(-N, N + 1):
        ck = coefficients[k + N]

        if ck == 0:
            continue

        a = 2 * mp.pi * k / L
        denom = a - tau

        # For this first test gamma1 is safely away from
        # the removable singularities.
        g = (exp_tL - 1) / (1j * denom)

        gp = (
            -L * exp_tL / denom
            - 1j * (exp_tL - 1) / (denom ** 2)
        )

        H += ck * g
        Hp += ck * gp

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * (1j * L / 2 * H + Hp)
        / mp.sqrt(L)
    )


print("F'(gamma1) numerical =",
      mp.nstr(F_even_deriv_numeric(gamma1), 40))

print("F'(gamma1) analytic  =",
      mp.nstr(F_even_deriv_analytic(gamma1), 40))

print("difference            =",
      mp.nstr(
          F_even_deriv_numeric(gamma1)
          - F_even_deriv_analytic(gamma1),
          20
      ))

print("\nEigenvector symmetry and normalisation:")

for k in range(0, N + 1):
    ck = coefficients[N + k]

    if k == 0:
        vk = ck
    else:
        vk = mp.sqrt(2) * ck

    print(
        k,
        "c_k =", mp.nstr(ck, 20),
        "v_k =", mp.nstr(vk, 20)
    )

vk_values = [coefficients[N]]

for k in range(1, N + 1):
    vk_values.append(mp.sqrt(2) * coefficients[N + k])

beta = L / (4 * mp.pi)

P = (
    vk_values[0] / beta**2
    + mp.sqrt(2) * sum(
        vk_values[k] / (k**2 + beta**2)
        for k in range(1, N + 1)
    )
)

print("beta =", mp.nstr(beta, 30))
print("P_N(v) =", mp.nstr(P, 40))

# ------------------------------------------------------------
# Canonical paper coefficients from repository coefficients
# ------------------------------------------------------------

v = [coefficients[N]]

for k in range(1, N + 1):
    v.append(mp.sqrt(2) * coefficients[N + k])

print("\nCanonical v norm:")
print(mp.nstr(
    mp.sqrt(v[0]**2 + sum(v[k]**2 for k in range(1, N + 1))),
    40
))


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


# ------------------------------------------------------------
# Reconstruct F(v) from the canonical basis responses
# ------------------------------------------------------------

F_reconstructed = sum(
    v[k] * F_from_canonical_basis(k, gamma1)
    for k in range(N + 1)
)

Fp_reconstructed = sum(
    v[k] * Fprime_from_canonical_basis(k, gamma1)
    for k in range(N + 1)
)

print("\nF(gamma1):")
print("direct       =", mp.nstr(F_even(gamma1), 40))
print("reconstructed=", mp.nstr(F_reconstructed, 40))
print("difference   =", mp.nstr(
    F_even(gamma1) - F_reconstructed, 30
))

print("\nF'(gamma1):")
print("direct       =", mp.nstr(F_even_deriv_analytic(gamma1), 40))
print("reconstructed=", mp.nstr(Fp_reconstructed, 40))
print("difference   =", mp.nstr(
    F_even_deriv_analytic(gamma1) - Fp_reconstructed, 30
))

# ============================================================
# Build the four linear functionals
# ============================================================

gamma1 = mp.im(mp.zetazero(1))
gamma2 = mp.im(mp.zetazero(2))

beta = L / (4 * mp.pi)

a = []
b = []
p = []
d = []

for k in range(N + 1):

    # F response
    ak = F_from_canonical_basis(k, gamma1)
    bk = F_from_canonical_basis(k, gamma2)

    # F' response
    dk = Fprime_from_canonical_basis(k, gamma1)

    # Pole functional
    if k == 0:
        pk = 1 / beta**2
    else:
        pk = mp.sqrt(2) / (k**2 + beta**2)

    a.append(ak)
    b.append(bk)
    p.append(pk)
    d.append(dk)


print("gamma1 =", mp.nstr(gamma1, 40))
print("gamma2 =", mp.nstr(gamma2, 40))

print("\nConstraint row norms:")
print("||a|| =", mp.nstr(mp.sqrt(sum(x*x for x in a)), 30))
print("||b|| =", mp.nstr(mp.sqrt(sum(x*x for x in b)), 30))
print("||p|| =", mp.nstr(mp.sqrt(sum(x*x for x in p)), 30))
print("||d|| =", mp.nstr(mp.sqrt(sum(x*x for x in d)), 30))

print("\nGround-state checks:")

print("F(gamma1) direct =",
      mp.nstr(F_even(gamma1), 30))

print("F(gamma1) row    =",
      mp.nstr(sum(a[k] * v[k] for k in range(N + 1)), 30))

print("\nF(gamma2) row    =",
      mp.nstr(sum(b[k] * v[k] for k in range(N + 1)), 30))

print("\nP(v) direct      =",
      mp.nstr(P, 30))

print("P(v) row         =",
      mp.nstr(sum(p[k] * v[k] for k in range(N + 1)), 30))

print("\nF'(gamma1) direct =",
      mp.nstr(F_even_deriv_analytic(gamma1), 30))

print("F'(gamma1) row    =",
      mp.nstr(sum(d[k] * v[k] for k in range(N + 1)), 30))

# ============================================================
# Constraint geometry
# ============================================================

C2 = mp.matrix(3, N + 1)

for k in range(N + 1):
    C2[0, k] = a[k]
    C2[1, k] = b[k]
    C2[2, k] = p[k]

G2 = C2 * C2.T

print("C C^T =")
for i in range(3):
    print([
        mp.nstr(G2[i, j], 30)
        for j in range(3)
    ])

print("\ndeterminant(C C^T) =")
print(mp.nstr(mp.det(G2), 40))

print("\nconstraint-row cosine similarities:")

for i, j, name in [
    (0, 1, "a,b"),
    (0, 2, "a,p"),
    (1, 2, "b,p"),
]:
    dot = sum(C2[i,k] * C2[j,k] for k in range(N + 1))
    ni = mp.sqrt(sum(C2[i,k]**2 for k in range(N + 1)))
    nj = mp.sqrt(sum(C2[j,k]**2 for k in range(N + 1)))

    print(
        name, "=",
        mp.nstr(dot / (ni * nj), 30)
    )

# ============================================================
# Projection of d onto the constraint nullspace
# ============================================================

d_vec = mp.matrix(d)

C2d = C2 * d_vec

# Solve (C2 C2^T) y = C2 d
y = mp.lu_solve(G2, C2d)

# Orthogonal projection of d onto ker(C2)
d_perp2 = d_vec - C2.T * y

Dmax = mp.sqrt(mp.fdot(d_perp2, d_perp2))

# Maximising unit vector
v_star = d_perp2 / Dmax

print("D_max =")
print(mp.nstr(Dmax, 60))

print("\n||d_perp|| =")
print(mp.nstr(mp.sqrt(mp.fdot(d_perp2, d_perp2)), 60))

print("\n||v_star|| =")
print(mp.nstr(mp.sqrt(mp.fdot(v_star, v_star)), 60))

print("\nConstraint residuals:")
print("F(gamma1) =", mp.nstr((C2 * v_star)[0], 50))
print("F(gamma2) =", mp.nstr((C2 * v_star)[1], 50))
print("P(v)      =", mp.nstr((C2 * v_star)[2], 50))

print("\nDerivative:")
print("d . v_star =",
      mp.nstr(mp.fdot(d_vec, v_star), 60))

def F_from_canonical_vector(v, tau):
    total = mp.mpf('0')

    for k in range(N + 1):
        total += v[k] * F_from_canonical_basis(k, tau)

    return total

print("\nF around gamma1:")
for delta in ["-0.1", "-0.01", "0", "0.01", "0.1"]:
    x = gamma1 + mp.mpf(delta)
    print(
        delta,
        mp.nstr(F_from_canonical_vector(v_star, x), 40)
    )

print("\nF around gamma2:")
for delta in ["-0.1", "-0.01", "0", "0.01", "0.1"]:
    x = gamma2 + mp.mpf(delta)
    print(
        delta,
        mp.nstr(F_from_canonical_vector(v_star, x), 40)
    )

C1 = mp.matrix(2, N + 1)

for k in range(N + 1):
    C1[0, k] = a[k]
    C1[1, k] = b[k]

G1 = C1 * C1.T

print("C C^T =")
for i in range(2):
    print([
        mp.nstr(G1[i, j], 30)
        for j in range(2)
    ])

print("\ndeterminant(C C^T) =")
print(mp.nstr(mp.det(G1), 40))

print("\nconstraint-row cosine similarities:")

# ============================================================
# Projection of d onto the constraint nullspace
# ============================================================

d_vec = mp.matrix(d)

C1d = C1 * d_vec

# Solve (C C^T) y = C d
y = mp.lu_solve(G1, C1d)

# Orthogonal projection of d onto ker(C)
d_perp1 = d_vec - C1.T * y

Dmax = mp.sqrt(mp.fdot(d_perp1, d_perp1))

# Maximising unit vector
v_star = d_perp1 / Dmax

print("D_max =")
print(mp.nstr(Dmax, 60))

# ============================================================
# D_0: pole-neutral only
# ============================================================

C0 = mp.matrix(1, N + 1)

for k in range(N + 1):
    C0[0, k] = p[k]

G0 = C0 * C0.T

d_vec = mp.matrix(d)

Cd0 = C0 * d_vec

y0 = mp.lu_solve(G0, Cd0)

d_perp0 = d_vec - C0.T * y0

D0 = mp.sqrt(mp.fdot(d_perp0, d_perp0))

v0_star = d_perp0 / D0

print("D_0 =")
print(mp.nstr(D0, 60))

print("\nconstraint residual:")
print(mp.nstr(mp.fdot(mp.matrix(p), v0_star), 50))

print("\nnorm:")
print(mp.nstr(mp.sqrt(mp.fdot(v0_star, v0_star)), 50))

def constrained_D(m):
    """
    Maximum |F'(gamma1)| with:
        P(v) = 0
        F(gamma_j) = 0 for j=1,...,m
        ||v|| = 1
    """

    # Constraint rows
    rows = [mp.matrix(p)]

    for j in range(1, m + 1):
        gamma = mp.im(mp.zetazero(j))

        row = mp.matrix(
            [F_from_canonical_basis(k, gamma)
             for k in range(N + 1)]
        )

        rows.append(row)

    nr = len(rows)

    C = mp.matrix(nr, N + 1)

    for i in range(nr):
        for k in range(N + 1):
            C[i, k] = rows[i][k]

    G = C * C.T

    # Solve G y = C d
    d_vec = mp.matrix(d)
    y = mp.lu_solve(G, C * d_vec)

    # Project d onto ker(C)
    d_perp = d_vec - C.T * y

    D = mp.sqrt(mp.fdot(d_perp, d_perp))

    return D, C, G, d_perp

D3, C3, G3, dperp3 = constrained_D(3)

print("D_3 =")
print(mp.nstr(D3, 60))

print("\ndet(C3 C3^T) =")
print(mp.nstr(mp.det(G3), 50))

for m in range(4):
    Dm, Cm, Gm, dpm = constrained_D(m)

    print(
        f"D_{m} = "
        f"{mp.nstr(Dm, 60)}"
    )
    print(
        f"det(G_{m}) = "
        f"{mp.nstr(mp.det(Gm), 40)}"
    )
    print()

print("Constraint row norms:")

norms = []

for i in range(C3.rows):
    ni = mp.sqrt(mp.fdot(C3[i, :], C3[i, :]))
    norms.append(ni)
    print(i, mp.nstr(ni, 50))

print("\nNormalised Gram matrix:")

R3 = mp.matrix(C3.rows, C3.rows)

for i in range(C3.rows):
    for j in range(C3.rows):
        R3[i, j] = (
            mp.fdot(C3[i, :], C3[j, :])
            / (norms[i] * norms[j])
        )
        print(mp.nstr(R3[i, j], 30), end="  ")
    print()

print("\ndet(normalised Gram) =")
print(mp.nstr(mp.det(R3), 50))

# Third-zero row
c3 = mp.matrix([
    C3[3, k]
    for k in range(N + 1)
])

# Project c3 into ker(C2)
G2 = C2 * C2.T
y3 = mp.lu_solve(G2, C2 * c3)

c3_perp = c3 - C2.T * y3

c3_perp_norm = mp.sqrt(mp.fdot(c3_perp, c3_perp))
d2_norm = mp.sqrt(mp.fdot(d_perp2, d_perp2))

dot = mp.fdot(c3_perp, d_perp2)

cos_theta = dot / (c3_perp_norm * d2_norm)

print("||c3_perp|| =")
print(mp.nstr(c3_perp_norm, 60))

print("\n||d_perp2|| =")
print(mp.nstr(d2_norm, 60))

print("\n< c3_perp, dperp2 > =")
print(mp.nstr(dot, 60))

print("\ncos(theta) =")
print(mp.nstr(cos_theta, 60))

print("\n|cos(theta)|^2 =")
print(mp.nstr(cos_theta**2, 60))

