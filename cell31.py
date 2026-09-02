"""
Cell 31 — Large-T quadrature forensic

Purpose
-------
Cell 30 showed that the dyadic tail contributions

    D(T) = integral_T^(2T) I(r) dr

remain positive and decay approximately like T^(-1), but with
T*D(T) continuing to drift upward.

Before interpreting that drift asymptotically, we need to establish
that the large-T quadrature is actually resolving the oscillatory
integrand.

This cell therefore performs a controlled quadrature audit at three
large-T intervals. Precision and interval subdivision are varied
independently.

The mathematical integrand is unchanged throughout:

    I(r) = h_+(r) * K_fourier(v_star, r, L)

No asymptotic approximation is introduced here.

The forensic ground state is the canonical cached state specified by
FORENSIC_GROUND_STATE in cell.py. The cached state was generated at
the specified generation precision (currently 150 dps), and is
decoded at the working precision used by each individual test.
"""

from __future__ import annotations

import time
import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    K_fourier,
    compute_L,
    get_ground_state,
)


# ---------------------------------------------------------------------------
# Forensic parameters
# ---------------------------------------------------------------------------

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]
T_ground = FORENSIC_GROUND_STATE["T"]
GROUND_STATE_DPS = FORENSIC_GROUND_STATE["dps"]

# Large-T intervals selected from Cell 30.
T_VALUES = [
    1310720,
    5242880,
    20971520,
]

# Precision sweep.
PRECISIONS = [80, 100, 120]

# Equal subdivisions.
#
# subdivisions=1 reproduces the ordinary mp.quad call used in Cell 30.
# Increasing subdivision provides an independent test of whether the
# oscillatory structure is being resolved adequately.
SUBDIVISIONS = [1, 4, 16, 64]


# ---------------------------------------------------------------------------
# Archimedean factor
# ---------------------------------------------------------------------------

def h_plus_cell31(r: mp.mpf) -> mp.mpf:
    """
    Archimedean h_+(r):

        Re psi(1/4 + i r/2) - log(pi)
    """
    return (
        mp.re(
            mp.digamma(
                mp.mpf("0.25") + 0.5j * r
            )
        )
        - mp.log(mp.pi)
    )


# ---------------------------------------------------------------------------
# Fixed forensic ground state
# ---------------------------------------------------------------------------

def get_forensic_state(dps):
    """
    Retrieve the canonical forensic ground state at the requested
    working precision.

    The cache entry is generated/certified at GROUND_STATE_DPS.
    Because get_ground_state() decodes at the caller's current
    mp.mp.dps, this function must be called inside an appropriate
    mp.workdps() context.
    """
    lambda_min, v_canonical, metadata = get_ground_state(
        c=c,
        N=N,
        T=T_ground,
        dps=GROUND_STATE_DPS,
    )

    L = compute_L(c)

    return lambda_min, v_canonical, L, metadata


# ---------------------------------------------------------------------------
# Integrand factory
# ---------------------------------------------------------------------------

def make_integrand(v_star, L):
    """
    Construct the scalar Archimedean tail integrand for a fixed
    working-precision ground state.
    """
    def integrand(r):
        return (
            h_plus_cell31(r)
            * K_fourier(v_star, r, L)
        )

    return integrand


# ---------------------------------------------------------------------------
# Integration helpers
# ---------------------------------------------------------------------------

def integrate_equal_subdivisions(
    integrand,
    a,
    b,
    subdivisions,
):
    """
    Integrate over [a,b], explicitly splitting it into equal pieces.

    subdivisions=1 is the ordinary mp.quad call.
    """
    if subdivisions == 1:
        return mp.quad(
            integrand,
            [a, b],
        )

    width = (b - a) / subdivisions

    total = mp.mpf("0")

    for k in range(subdivisions):
        left = a + k * width
        right = a + (k + 1) * width

        total += mp.quad(
            integrand,
            [left, right],
        )

    return total


def timed_call(fn):
    start = time.perf_counter()

    value = fn()

    elapsed = (
        time.perf_counter()
        - start
    )

    return value, elapsed


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

print("=" * 78)
print("CELL 31 — LARGE-T QUADRATURE FORENSIC")
print("=" * 78)
print()
print(f"c                         = {c}")
print(f"N                         = {N}")
print(f"T_ground                 = {T_ground}")
print(f"ground-state generation  = {GROUND_STATE_DPS} dps")
print(f"precision sweep           = {PRECISIONS}")
print(f"subdivisions              = {SUBDIVISIONS}")
print()
print(
    "The ground state is retrieved separately inside each working-"
    "precision context."
)
print(
    "All retrievals refer to the same immutable forensic cache entry."
)
print()


# ---------------------------------------------------------------------------
# 1. Precision sweep
# ---------------------------------------------------------------------------

print("=" * 78)
print("1. PRECISION SWEEP — UNSUBDIVIDED mp.quad")
print("=" * 78)
print()

precision_results = {}

for dps in PRECISIONS:

    print(f"WORKING PRECISION = {dps} dps")
    print("-" * 78)

    with mp.workdps(dps):

        (
            lambda_min,
            v_star,
            L,
            ground_meta,
        ) = get_forensic_state(dps)

        print(
            f"ground-state working dps = "
            f"{ground_meta['working_dps']}"
        )
        print(
            f"ground-state generation dps = "
            f"{ground_meta['generation_dps']}"
        )
        print(
            f"lambda_min = "
            f"{mp.nstr(lambda_min, 25)}"
        )
        print(
            f"||v_star|| = "
            f"{mp.nstr(mp.norm(v_star), 25)}"
        )
        print(
            f"L = log(c) = "
            f"{mp.nstr(L, 30)}"
        )
        print()

        integrand = make_integrand(
            v_star,
            L,
        )

        precision_results[dps] = {}

        for T_int in T_VALUES:

            T = mp.mpf(T_int)
            a = T
            b = 2 * T

            value, elapsed = timed_call(
                lambda: mp.quad(
                    integrand,
                    [a, b],
                )
            )

            precision_results[dps][T_int] = value

            print(
                f"T={T_int:,}  "
                f"D(T)={mp.nstr(value, 35)}  "
                f"time={elapsed:8.2f}s"
            )

    print()


# ---------------------------------------------------------------------------
# Compare precision results
# ---------------------------------------------------------------------------

print("=" * 78)
print("2. PRECISION DIFFERENCES")
print("=" * 78)
print()

reference_dps = max(PRECISIONS)

for T_int in T_VALUES:

    reference = (
        precision_results[reference_dps][T_int]
    )

    print(f"T = {T_int:,}")
    print("-" * 78)

    for dps in PRECISIONS:

        value = precision_results[dps][T_int]

        delta = value - reference

        if reference != 0:
            relative = abs(
                delta / reference
            )
        else:
            relative = mp.inf

        print(
            f"{dps:3d} dps: "
            f"abs diff = {mp.nstr(abs(delta), 12)}   "
            f"rel diff = {mp.nstr(relative, 12)}"
        )

    print()


# ---------------------------------------------------------------------------
# 3. Subdivision sweep
# ---------------------------------------------------------------------------

print("=" * 78)
print("3. SUBDIVISION SWEEP")
print("=" * 78)
print()

subdivision_results = {}

for dps in [80, 120]:

    print(f"WORKING PRECISION = {dps} dps")
    print("-" * 78)

    with mp.workdps(dps):

        (
            lambda_min,
            v_star,
            L,
            ground_meta,
        ) = get_forensic_state(dps)

        integrand = make_integrand(
            v_star,
            L,
        )

        subdivision_results[dps] = {}

        for T_int in T_VALUES:

            T = mp.mpf(T_int)
            a = T
            b = 2 * T

            print(
                f"Interval [{T_int:,}, {2*T_int:,}]"
            )

            results_at_dps = {}

            for subdivisions in SUBDIVISIONS:

                value, elapsed = timed_call(
                    lambda s=subdivisions:
                        integrate_equal_subdivisions(
                            integrand,
                            a,
                            b,
                            s,
                        )
                )

                results_at_dps[subdivisions] = value

                print(
                    f"  subdivisions={subdivisions:3d}  "
                    f"D(T)={mp.nstr(value, 35)}  "
                    f"time={elapsed:8.2f}s"
                )

            subdivision_results[dps][T_int] = (
                results_at_dps
            )

            print()

    print()


# ---------------------------------------------------------------------------
# 4. Subdivision convergence
# ---------------------------------------------------------------------------

print("=" * 78)
print("4. SUBDIVISION CONVERGENCE")
print("=" * 78)
print()

for dps in [80, 120]:

    print(f"WORKING PRECISION = {dps} dps")
    print("-" * 78)

    for T_int in T_VALUES:

        results = (
            subdivision_results[dps][T_int]
        )

        reference = results[
            max(SUBDIVISIONS)
        ]

        print(f"T = {T_int:,}")

        for subdivisions in SUBDIVISIONS:

            value = results[subdivisions]

            delta = (
                value
                - reference
            )

            if reference != 0:
                relative = abs(
                    delta / reference
                )
            else:
                relative = mp.inf

            print(
                f"  {subdivisions:3d} -> "
                f"{max(SUBDIVISIONS):3d}: "
                f"rel diff = "
                f"{mp.nstr(relative, 12)}"
            )

        print()

    print()


# ---------------------------------------------------------------------------
# 5. Cross-comparison of the two strongest calculations
# ---------------------------------------------------------------------------

print("=" * 78)
print("5. 80/120-dps REFINED CROSS-CHECK")
print("=" * 78)
print()

for T_int in T_VALUES:

    value_80 = (
        subdivision_results[80][T_int][64]
    )

    value_120 = (
        subdivision_results[120][T_int][64]
    )

    delta = (
        value_80
        - value_120
    )

    relative = abs(
        delta / value_120
    )

    print(f"T = {T_int:,}")
    print(
        f"  80 dps / 64 subdivisions  = "
        f"{mp.nstr(value_80, 35)}"
    )
    print(
        f" 120 dps / 64 subdivisions  = "
        f"{mp.nstr(value_120, 35)}"
    )
    print(
        f"  relative difference       = "
        f"{mp.nstr(relative, 15)}"
    )
    print()


# ---------------------------------------------------------------------------
# Interpretation
# ---------------------------------------------------------------------------

print("=" * 78)
print("6. INTERPRETATION GUIDE")
print("=" * 78)
print(
    """
The purpose of this cell is purely numerical validation.

The strongest evidence that the Cell 30 tail values are trustworthy
would be:

  1. 80, 100 and 120 dps give essentially the same D(T);
  2. increasing subdivision gives essentially the same D(T);
  3. the 80-dps and 120-dps refined calculations agree.

If increasing dps changes D(T) materially while subdivision does not,
that points toward precision loss.

If subdivision changes D(T) materially while dps does not,
that points toward unresolved oscillatory structure.

If both precision and subdivision changes become negligible, then
the Cell 30 large-T drift is much less likely to be a straightforward
quadrature artefact.

No asymptotic conclusion is drawn by this cell.
"""
)
