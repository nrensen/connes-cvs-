# ============================================================
# CELL 9 — FINITE-DIMENSIONAL WEIL DICTIONARY AUDIT
#
# Goal:
#
# Establish explicitly what a coefficient vector in the
# trigonometric Galerkin basis represents.
#
# We test the chain
#
#     coefficients
#          |
#          v
#     f_v(t)
#          |
#          v
#     Fourier transform F_v(tau)
#          |
#          v
#     translated test function
#          |
#          v
#     Weil quadratic form
#
# against the actual Galerkin matrix Q.
#
# This cell is deliberately independent of the ground-state
# machinery.  The first test uses arbitrary coefficient vectors,
# which makes it much harder for a coincidental ground-state
# property to hide a dictionary error.
# ============================================================

import mpmath as mp

from cell import (
    compute_L,
    canonical_to_full,
    canonical_pairs,
)

from connes_cvs import (
    build_galerkin_matrix,
    compute_ground_state,
)

from connes_cvs.operator import (
    extract_zeros,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 70

c = 13
N = 6
T = 40

L = compute_L(c)

omega = 2 * mp.pi / L

print("=" * 70)
print("CELL 9 — FINITE-DIMENSIONAL WEIL DICTIONARY AUDIT")
print("=" * 70)
print()
print("Parameters:")
print(f"c = {c}")
print(f"N = {N}")
print(f"T = {T}")
print(f"dps = {mp.mp.dps}")
print(f"L = {mp.nstr(L, 60)}")
print(f"2*pi/L = {mp.nstr(omega, 60)}")
print()


# ============================================================
# 1. CANONICAL BASIS FUNCTIONS
#
# The canonical real-even basis is
#
#   phi_0(t) = 1/sqrt(L)
#
#   phi_k(t) =
#       sqrt(2/L) cos(2*pi*k*t/L),  k >= 1.
#
# This follows from the repository's canonical/full convention:
#
#   u_0 = v_0
#   u_{+k} = u_{-k} = v_k/sqrt(2).
# ============================================================

def phi(k, t):
    k = int(k)

    if k == 0:
        return 1 / mp.sqrt(L)

    return mp.sqrt(2 / L) * mp.cos(2 * mp.pi * k * t / L)


# ============================================================
# 2. TEST VECTORS
#
# Use several vectors, including basis vectors and a generic
# linear combination.
# ============================================================

def make_basis_vector(k):
    v = mp.matrix(N + 1, 1)
    v[k] = 1
    return v


test_vectors = []

for k in range(min(N + 1, 4)):
    test_vectors.append(
        (f"basis k={k}", make_basis_vector(k))
    )

v_generic = mp.matrix(N + 1, 1)

raw = [
    "0.73",
    "-0.41",
    "0.19",
    "0.11",
    "-0.07",
    "0.05",
    "-0.03",
]

for k in range(N + 1):
    v_generic[k] = mp.mpf(raw[k])

test_vectors.append(("generic vector", v_generic))


# ============================================================
# 3. DIRECT REAL TEST FUNCTION
#
#     f_v(t) = sum_k v_k phi_k(t)
#
# This is the actual element of the finite-dimensional
# test-function space.
# ============================================================

def f_direct(v, t):
    total = mp.mpf("0")

    for k in range(N + 1):
        total += v[k] * phi(k, t)

    return total


# ============================================================
# 4. FULL FOURIER COEFFICIENTS
#
# In the complex exponential basis
#
#     e_k(t) = exp(2*pi*i*k*t/L)/sqrt(L)
#
# the real-even canonical vector corresponds to
#
#     c_0 = v_0
#
#     c_{+k} = c_{-k} = v_k/sqrt(2).
#
# This is exactly the convention used by cell.py.
# ============================================================

def full_coefficients(v):
    u = canonical_to_full(v)

    return {
        m: mp.mpf(u[m + N])
        for m in range(-N, N + 1)
    }


# ============================================================
# 5. FOURIER TRANSFORM
#
# Define
#
#     F_v(tau) = integral_0^L f_v(t)
#                         exp(i*tau*t) dt.
#
# For the individual exponential basis function this is
#
#     integral_0^L exp(i*(tau + 2*pi*k/L)t)/sqrt(L) dt.
#
# We evaluate the transform in TWO ways:
#
#   A. direct numerical t-quadrature
#   B. closed-form exponential expression
#
# This establishes the Fourier-side dictionary independently.
# ============================================================

def F_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda t: f_direct(v, t) * mp.exp(1j * tau * t),
        [0, L],
    )


def exponential_integral(a):
    a = mp.mpf(a)

    if a == 0:
        return L

    return (
        mp.exp(1j * a * L) - 1
    ) / (1j * a)


def F_closed(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = 0j

    for k in range(-N, N + 1):
        a = tau + 2 * mp.pi * k / L

        total += (
            u[k]
            * exponential_integral(a)
            / mp.sqrt(L)
        )

    return total


# ============================================================
# 6. FOURIER TRANSFORM AUDIT
# ============================================================

print("-" * 70)
print("1. FOURIER TRANSFORM DICTIONARY")
print("-" * 70)
print()

tau_values = [
    mp.mpf("0"),
    mp.mpf("0.3"),
    mp.mpf("1"),
    omega,
    mp.mpf("3"),
    mp.mpf("7"),
]

max_F_error = mp.mpf("0")

for name, v in test_vectors:

    print(name)

    for tau in tau_values:

        fd = F_direct(v, tau)
        fc = F_closed(v, tau)

        err = abs(fd - fc)

        max_F_error = max(max_F_error, err)

        print(
            f"  tau={mp.nstr(tau, 12)}"
        )
        print(
            f"    direct = {mp.nstr(fd, 35)}"
        )
        print(
            f"    closed = {mp.nstr(fc, 35)}"
        )
        print(
            f"    error  = {mp.nstr(err, 15)}"
        )

    print()


print(
    "Maximum |F_direct - F_closed| =",
    mp.nstr(max_F_error, 30),
)
print()


# ============================================================
# 7. BASIS NORMALISATION
#
# Verify directly that the canonical basis is orthonormal:
#
#     integral_0^L phi_j(t) phi_k(t) dt = delta_jk.
# ============================================================

print("-" * 70)
print("2. BASIS ORTHONORMALITY")
print("-" * 70)
print()

max_basis_error = mp.mpf("0")

for j in range(N + 1):
    for k in range(N + 1):

        value = mp.quad(
            lambda t: phi(j, t) * phi(k, t),
            [0, L],
        )

        target = mp.mpf(1) if j == k else mp.mpf(0)

        err = abs(value - target)

        max_basis_error = max(max_basis_error, err)

print(
    "Maximum orthonormality error =",
    mp.nstr(max_basis_error, 30),
)
print()


# ============================================================
# 8. GALERKIN MATRIX
# ============================================================

print("-" * 70)
print("3. GALERKIN MATRIX")
print("-" * 70)
print()

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=mp.mp.dps,
)

DIM = 2 * N + 1

print(f"Matrix dimension = {DIM} x {DIM}")
print()


# ============================================================
# 9. QUADRATIC FORM FROM Q
#
# For a canonical vector v, lift to full symmetric coordinates
# u and calculate
#
#     A_Q(v) = u^T Q u.
#
# This is the direct repository quadratic form.
# ============================================================

def Q_quadratic(v):
    u = canonical_to_full(v)
    return mp.fdot(u, Q * u)


# ============================================================
# 10. BASIS-BY-BASIS QUADRATIC FORM
#
# Compare the canonical quadratic form against the corresponding
# matrix projected into the canonical even sector.
#
# Constructing Q_even explicitly also gives a useful independent
# check on all sqrt(2) factors.
# ============================================================

V_even = mp.matrix(DIM, N + 1)

V_even[N, 0] = 1

inv_sqrt2 = 1 / mp.sqrt(2)

for k in range(1, N + 1):
    V_even[N + k, k] = inv_sqrt2
    V_even[N - k, k] = inv_sqrt2

Q_even = V_even.T * Q * V_even


def Q_even_quadratic(v):
    return mp.fdot(v, Q_even * v)


print("-" * 70)
print("4. CANONICAL / FULL QUADRATIC-FORM DICTIONARY")
print("-" * 70)
print()

max_Q_dictionary_error = mp.mpf("0")

for name, v in test_vectors:

    q_full = Q_quadratic(v)
    q_even = Q_even_quadratic(v)

    err = abs(q_full - q_even)

    max_Q_dictionary_error = max(
        max_Q_dictionary_error,
        err,
    )

    print(name)
    print(
        "  full-space Q form =",
        mp.nstr(q_full, 45),
    )
    print(
        "  canonical Q form =",
        mp.nstr(q_even, 45),
    )
    print(
        "  error             =",
        mp.nstr(err, 25),
    )
    print()


print(
    "Maximum Q dictionary error =",
    mp.nstr(max_Q_dictionary_error, 30),
)
print()


# ============================================================
# 11. GROUND STATE
#
# Now use the actual lowest even eigenvector.
#
# This is not yet an RH test. We simply reconstruct the function
# represented by the numerical ground state.
# ============================================================

print("-" * 70)
print("5. GROUND-STATE RECONSTRUCTION")
print("-" * 70)
print()

lambda_min, u_star = compute_ground_state(Q)

v_star = mp.matrix(N + 1, 1)

v_star[0] = u_star[N, 0]

for k in range(1, N + 1):
    v_star[k] = (
        mp.sqrt(2)
        * u_star[N + k, 0]
    )

v_norm = mp.sqrt(
    mp.fdot(v_star, v_star)
)

for k in range(N + 1):
    v_star[k] /= v_norm


print(
    "lambda_min =",
    mp.nstr(lambda_min, 50),
)

print(
    "||v_star|| =",
    mp.nstr(mp.sqrt(mp.fdot(v_star, v_star)), 30),
)

print()
print("Canonical coefficients:")

for k in range(N + 1):
    print(
        f"  v[{k}] =",
        mp.nstr(v_star[k], 40),
    )

print()


# ============================================================
# 12. GROUND-STATE FOURIER RESPONSE
#
# Evaluate F_v(tau) at a selection of frequencies.
#
# We report both the complex Fourier transform and its real
# part.  Because v is real-even, the relevant response has
# the expected conjugation/parity structure.
# ============================================================

print("-" * 70)
print("6. GROUND-STATE FOURIER RESPONSE")
print("-" * 70)
print()

for tau in [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("3"),
    mp.mpf("5"),
    mp.mpf("10"),
]:

    F = F_closed(v_star, tau)
    Fm = F_closed(v_star, -tau)

    print(
        "tau =",
        mp.nstr(tau, 15),
    )
    print(
        "  F(tau)   =",
        mp.nstr(F, 40),
    )
    print(
        "  F(-tau)  =",
        mp.nstr(Fm, 40),
    )
    print(
        "  F(-tau)-conj(F(tau)) =",
        mp.nstr(
            Fm - mp.conj(F),
            20,
        ),
    )
    print()


# ============================================================
# 13. REPOSITORY ZERO EXTRACTION
#
# The repository documents the spectral test function as
#
#   F_even(tau)
#      = Re[
#          exp(i*tau*L/2)
#          sum_k c_k g_k(tau)
#        ] / sqrt(L)
#
# where
#
#   g_k(tau)
#      = (exp(-i*tau*L)-1)
#        / (i*(2*pi*k/L - tau)).
#
# Reconstruct this explicitly and compare it with the equivalent
# Fourier representation above.
# ============================================================

def g_repo(k, tau):
    tau = mp.mpf(tau)

    denom = 2 * mp.pi * k / L - tau

    if denom == 0:
        return L

    return (
        mp.exp(-1j * tau * L) - 1
    ) / (
        1j * denom
    )


def F_repo_even(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = 0j

    for k in range(-N, N + 1):
        total += u[k] * g_repo(k, tau)

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * total
    ) / mp.sqrt(L)


print("-" * 70)
print("7. REPOSITORY SPECTRAL-FUNCTION DICTIONARY")
print("-" * 70)
print()

max_F_repo_error = mp.mpf("0")

for tau in [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("3"),
    mp.mpf("5"),
    mp.mpf("10"),
]:

    F = F_closed(v_star, tau)

    repo = F_repo_even(v_star, tau)

    # The repository object is a centred real-even transform.
    # Compare it with the corresponding centred real response.
    centred = mp.re(
        mp.exp(1j * tau * L / 2) * F
    )

    err = abs(centred / mp.sqrt(L) - repo)

    max_F_repo_error = max(
        max_F_repo_error,
        err,
    )

    print(
        f"tau={mp.nstr(tau, 15)}"
    )
    print(
        "  centred direct =",
        mp.nstr(centred / mp.sqrt(L), 40),
    )
    print(
        "  repository     =",
        mp.nstr(repo, 40),
    )
    print(
        "  error          =",
        mp.nstr(err, 20),
    )
    print()


print(
    "Maximum spectral dictionary error =",
    mp.nstr(max_F_repo_error, 30),
)
print()


# ============================================================
# SUMMARY
# ============================================================

print("=" * 70)
print("CELL 9 SUMMARY")
print("=" * 70)
print()

print(
    "Maximum Fourier-transform error =",
    mp.nstr(max_F_error, 30),
)

print(
    "Maximum basis orthonormality error =",
    mp.nstr(max_basis_error, 30),
)

print(
    "Maximum full/canonical Q-form error =",
    mp.nstr(max_Q_dictionary_error, 30),
)

print(
    "Maximum spectral-function dictionary error =",
    mp.nstr(max_F_repo_error, 30),
)

print()
print(
    "The purpose of Cell 9 is dictionary validation, not an RH"
)
print(
    "claim.  A successful run establishes exactly what finite"
)
print(
    "coefficient vectors represent on the Fourier/Weil side."
)
print()
print("=" * 70)
print("END CELL 9")
print("=" * 70)
