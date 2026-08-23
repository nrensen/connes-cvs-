# ============================================================
# CELL 11 — EXACT FOURIER / EXTRACTION DICTIONARY AUDIT
#
# Purpose:
#
# Determine exactly what mathematical Fourier transform is
# implemented by connes_cvs.operator.extract_zeros().
#
# In particular, distinguish:
#
#   F_+(tau) = integral_0^L f(t) exp(+i tau t) dt
#   F_-(tau) = integral_0^L f(t) exp(-i tau t) dt
#
# and their centred versions.
#
# This cell deliberately does NOT perform zero finding.
# ============================================================

import mpmath as mp

from cell import (
    compute_L,
    canonical_to_full,
)

from connes_cvs import (
    build_galerkin_matrix,
    compute_ground_state,
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
print("CELL 11 — EXACT FOURIER / EXTRACTION DICTIONARY AUDIT")
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
# 1. CANONICAL BASIS
# ============================================================

def phi(k, t):
    k = int(k)

    if k == 0:
        return 1 / mp.sqrt(L)

    return mp.sqrt(2 / L) * mp.cos(2 * mp.pi * k * t / L)


def f_direct(v, t):
    total = mp.mpf("0")

    for k in range(N + 1):
        total += v[k] * phi(k, t)

    return total


def full_coefficients(v):
    u = canonical_to_full(v, N)

    return {
        k: mp.mpf(u[k + N])
        for k in range(-N, N + 1)
    }


# ============================================================
# 2. DIRECT FOURIER TRANSFORMS
# ============================================================

def F_plus_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda t: f_direct(v, t) * mp.exp(1j * tau * t),
        [0, L],
    )


def F_minus_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda t: f_direct(v, t) * mp.exp(-1j * tau * t),
        [0, L],
    )


# ============================================================
# 3. CLOSED FORM F_+
#
# e_k(t) = exp(+i omega k t) / sqrt(L)
#
# F_+(tau)
#   = sum_k c_k / sqrt(L)
#       * integral exp(i(tau + omega k)t) dt
# ============================================================

def exponential_integral_plus(a):
    a = mp.mpf(a)

    if a == 0:
        return L

    return (
        mp.exp(1j * a * L) - 1
    ) / (1j * a)


def F_plus_closed(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = mp.mpc(0)

    for k in range(-N, N + 1):
        a = tau + omega * k

        total += (
            u[k]
            * exponential_integral_plus(a)
            / mp.sqrt(L)
        )

    return total


# ============================================================
# 4. CLOSED FORM F_-
#
# F_-(tau)
#   = sum_k c_k / sqrt(L)
#       * integral exp(i(omega k - tau)t) dt
#
# This is exactly the structure used by extract_zeros().
# ============================================================

def exponential_integral_minus(k, tau):
    tau = mp.mpf(tau)

    denom = omega * k - tau

    if denom == 0:
        return L

    return (
        mp.exp(1j * denom * L) - 1
    ) / (1j * denom)


def F_minus_closed(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = mp.mpc(0)

    for k in range(-N, N + 1):
        total += (
            u[k]
            * exponential_integral_minus(k, tau)
            / mp.sqrt(L)
        )

    return total


# ============================================================
# 5. REPOSITORY EXTRACTION EXPRESSION
#
# This reproduces the expression in operator.py:
#
#   exp(+i tau L/2) *
#       sum c_k g_k(tau) / sqrt(L)
#
# where
#
#   g_k(tau)
#       = (exp(-i tau L)-1)
#         / (i(omega k - tau)).
#
# ============================================================

def repository_expression(v, tau):
    tau = mp.mpf(tau)

    u = full_coefficients(v)

    total = mp.mpc(0)

    exp_tL = mp.exp(-1j * tau * L)

    for k in range(-N, N + 1):
        denom = omega * k - tau

        if denom == 0:
            term = mp.mpc(L, 0)

        else:
            term = (
                exp_tL - 1
            ) / (1j * denom)

        total += u[k] * term

    total /= mp.sqrt(L)

    return mp.exp(1j * tau * L / 2) * total


# ============================================================
# 6. CENTRED TRANSFORMS
#
# If
#
#   g(x) = f(x + L/2)
#
# then
#
#   integral_{-L/2}^{L/2}
#       g(x) exp(-i tau x) dx
#
# = exp(+i tau L/2) F_-(tau)
#
# and
#
#   integral_{-L/2}^{L/2}
#       g(x) exp(+i tau x) dx
#
# = exp(-i tau L/2) F_+(tau).
# ============================================================

def centred_minus_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda x:
            f_direct(v, x + L / 2)
            * mp.exp(-1j * tau * x),
        [-L / 2, L / 2],
    )


def centred_plus_direct(v, tau):
    tau = mp.mpf(tau)

    return mp.quad(
        lambda x:
            f_direct(v, x + L / 2)
            * mp.exp(+1j * tau * x),
        [-L / 2, L / 2],
    )


# ============================================================
# 7. TEST VECTORS
# ============================================================

test_vectors = []

for k in range(min(N + 1, 4)):
    v = mp.matrix(N + 1, 1)
    v[k] = 1
    test_vectors.append((f"basis k={k}", v))


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


tau_values = [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    omega,
    mp.mpf("3"),
    mp.mpf("5"),
]


# ============================================================
# 8. DIRECT VS CLOSED TRANSFORMS
# ============================================================

print("-" * 70)
print("1. DIRECT / CLOSED FOURIER TRANSFORMS")
print("-" * 70)
print()

max_plus_error = mp.mpf("0")
max_minus_error = mp.mpf("0")

for name, v in test_vectors:

    print(name)

    for tau in tau_values:

        dp = F_plus_direct(v, tau)
        cp = F_plus_closed(v, tau)

        dm = F_minus_direct(v, tau)
        cm = F_minus_closed(v, tau)

        ep = abs(dp - cp)
        em = abs(dm - cm)

        max_plus_error = max(max_plus_error, ep)
        max_minus_error = max(max_minus_error, em)

        print(
            f" tau={mp.nstr(tau, 12)}"
        )

        print(
            "   F+ direct  =",
            mp.nstr(dp, 30)
        )

        print(
            "   F+ closed  =",
            mp.nstr(cp, 30)
        )

        print(
            "   |error|    =",
            mp.nstr(ep, 12)
        )

        print(
            "   F- direct  =",
            mp.nstr(dm, 30)
        )

        print(
            "   F- closed  =",
            mp.nstr(cm, 30)
        )

        print(
            "   |error|    =",
            mp.nstr(em, 12)
        )

    print()

print(
    "Maximum F+ direct/closed error =",
    mp.nstr(max_plus_error, 30),
)

print(
    "Maximum F- direct/closed error =",
    mp.nstr(max_minus_error, 30),
)

print()


# ============================================================
# 9. REPOSITORY EXPRESSION VS F_-
# ============================================================

print("-" * 70)
print("2. REPOSITORY EXPRESSION VS CENTRED F_-")
print("-" * 70)
print()

max_repo_error = mp.mpf("0")

v = v_generic

for tau in tau_values:

    repo = repository_expression(v, tau)

    expected = (
        mp.exp(1j * tau * L / 2)
        * F_minus_direct(v, tau)
    )

    err = abs(repo - expected)

    max_repo_error = max(max_repo_error, err)

    print(
        "tau =",
        mp.nstr(tau, 20)
    )

    print(
        "  repository expression =",
        mp.nstr(repo, 35)
    )

    print(
        "  exp(+i*tau*L/2) F_-   =",
        mp.nstr(expected, 35)
    )

    print(
        "  error =",
        mp.nstr(err, 15)
    )

    print()

print(
    "Maximum repository-vs-centred-F_- error =",
    mp.nstr(max_repo_error, 30),
)

print()


# ============================================================
# 10. CENTRED TRANSFORM AUDIT
# ============================================================

print("-" * 70)
print("3. CENTRED TRANSFORM IDENTITIES")
print("-" * 70)
print()

max_centred_minus_error = mp.mpf("0")
max_centred_plus_error = mp.mpf("0")

v = v_generic

for tau in tau_values:

    cm_direct = centred_minus_direct(v, tau)

    cm_expected = (
        mp.exp(1j * tau * L / 2)
        * F_minus_direct(v, tau)
    )

    cp_direct = centred_plus_direct(v, tau)

    cp_expected = (
        mp.exp(-1j * tau * L / 2)
        * F_plus_direct(v, tau)
    )

    em = abs(cm_direct - cm_expected)
    ep = abs(cp_direct - cp_expected)

    max_centred_minus_error = max(
        max_centred_minus_error,
        em,
    )

    max_centred_plus_error = max(
        max_centred_plus_error,
        ep,
    )

    print(
        "tau =",
        mp.nstr(tau, 20)
    )

    print(
        "  centred minus error =",
        mp.nstr(em, 15)
    )

    print(
        "  centred plus error  =",
        mp.nstr(ep, 15)
    )

print()

print(
    "Maximum centred-minus identity error =",
    mp.nstr(max_centred_minus_error, 30),
)

print(
    "Maximum centred-plus identity error =",
    mp.nstr(max_centred_plus_error, 30),
)

print()


# ============================================================
# 11. THE CRITICAL COMPARISON
#
# Compare the repository expression against BOTH possible
# centred transforms.
# ============================================================

print("-" * 70)
print("4. CRITICAL REPOSITORY DICTIONARY TEST")
print("-" * 70)
print()

max_repo_minus = mp.mpf("0")
max_repo_plus = mp.mpf("0")

v = v_generic

for tau in tau_values:

    repo = repository_expression(v, tau)

    centred_minus = centred_minus_direct(v, tau)
    centred_plus = centred_plus_direct(v, tau)

    err_minus = abs(
        repo - centred_minus
    )

    err_plus = abs(
        repo - centred_plus
    )

    max_repo_minus = max(
        max_repo_minus,
        err_minus,
    )

    max_repo_plus = max(
        max_repo_plus,
        err_plus,
    )

    print(
        "tau =",
        mp.nstr(tau, 20)
    )

    print(
        "  |repo - centred F_-| =",
        mp.nstr(err_minus, 20)
    )

    print(
        "  |repo - centred F_+| =",
        mp.nstr(err_plus, 20)
    )

    print()

print(
    "Maximum |repo - centred F_-| =",
    mp.nstr(max_repo_minus, 30),
)

print(
    "Maximum |repo - centred F_+| =",
    mp.nstr(max_repo_plus, 30),
)

print()


# ============================================================
# 12. EVENNESS CHECK
#
# Because the canonical functions are real and even around the
# midpoint after translation, the correctly centred Fourier
# transform should be real and even.
# ============================================================

print("-" * 70)
print("5. CENTRED EVENNESS / REALITY CHECK")
print("-" * 70)
print()

max_real_minus = mp.mpf("0")
max_even_minus = mp.mpf("0")

for tau in tau_values:

    cm = centred_minus_direct(v_generic, tau)
    cm_neg = centred_minus_direct(v_generic, -tau)

    real_error = abs(mp.im(cm))
    even_error = abs(cm - cm_neg)

    max_real_minus = max(
        max_real_minus,
        real_error,
    )

    max_even_minus = max(
        max_even_minus,
        even_error,
    )

    print(
        "tau =",
        mp.nstr(tau, 20)
    )

    print(
        "  |Im centred F_-| =",
        mp.nstr(real_error, 15)
    )

    print(
        "  |F_-(tau)-F_-(-tau)| =",
        mp.nstr(even_error, 15)
    )

print()

print(
    "Maximum imaginary part =",
    mp.nstr(max_real_minus, 30),
)

print(
    "Maximum evenness error =",
    mp.nstr(max_even_minus, 30),
)

print()


# ============================================================
# 13. SYMBOLIC CONCLUSION PRINT
# ============================================================

print("=" * 70)
print("CELL 11 SUMMARY")
print("=" * 70)
print()
print(
    "The repository extract_zeros expression is mathematically"
)
print(
    "the centred NEGATIVE-frequency Fourier transform if"
)
print(
    "the corresponding errors above are at numerical precision."
)
print()
print(
    "This cell therefore distinguishes a genuine repository"
)
print(
    "dictionary error from a Cell-9 implementation error."
)
print()
print("=" * 70)
print("END CELL 11")
print("=" * 70)
