#!/usr/bin/env python3
"""
================================================================================
CELL 73 — AUDIT OF POSITIVE REGULARIZED STIELTJES FUNCTION H(mu),
         EXACT MODAL RATIO Q_j = H(mu_1)/H(mu_j), AND INTERLACING POLE SPECTRUM
================================================================================

PURPOSE:
--------
Perform a dedicated numerical audit of the positive regularized Stieltjes function
H(mu) = (mu - lam_0)^2 * G_d'(mu) and verify the exact modal ratio representation:
    Q_j = H(mu_1) / H(mu_j)   <===>   T_j = K_j * (H(mu_1) / H(mu_j))
established in Paper NR2 Proposition 8.24 (Milestone M20).

Background & Analytical Framework:
----------------------------------
In Cell 71 and Cell 72, the three-factor decomposition:
    T_j = K_j * G_j * S_j = K_j * Q_j
revealed that overlap growth G_2 ~ 10^5 and spectral suppression S_2 ~ 10^{-11}
compete across 16 orders of magnitude. 

Proposition 8.24 eliminates both intermediate scales simultaneously by introducing
the positive regularized Stieltjes function:
    H(mu) = D_0^2 + sum_{k=1}^N d_k^2 * ((mu - lam_0) / (E_k - mu))^2,
which is manifestly positive and finite for all mu not in sigma(Q_even \\ {lam_0}).

This script audits:
1. Exact Identity Verification:
       H(mu_j) == (mu_j - lam_0)^2 * G_d'(mu_j)
       Q_j     == H(mu_1) / H(mu_j)
       T_j     == K_j * (H(mu_1) / H(mu_j))
   verifying that the residual is zero to within mpmath numerical precision.

2. Term-by-Term Pole Decomposition of H(mu_j):
   Decompose H(mu_j) across all even modes k in {0, ..., N}:
       H_0(mu_j) = D_0^2,
       H_k(mu_j) = d_k^2 * ((mu_j - lam_0) / (E_k - mu_j))^2   (k >= 1).
   Determine:
   - The percentage contribution of each pole k.
   - The dominant pole k_dom(j) = argmax_k H_k(mu_j) and its share.
   - The nearest even pole k_near(j) = argmin_{k >= 1} |E_k - mu_j| and the
     spectral distance |E_{k_near} - mu_j|.
   - The sharpness of the single-pole lower bound:
         H(mu_j) >= H_{k_near}(mu_j) = d_{k_near}^2 * ((mu_j - lam_0) / (E_{k_near} - mu_j))^2.

3. Semiclassical Inverse Filtering Growth:
   Evaluate the inverse filtering ratios across discrete dimensions N in {8, 12, 16, 20, 24}:
       R_{21}(N) = H(mu_2) / H(mu_1) = Q_2^{-1},
       R_{31}(N) = H(mu_3) / H(mu_1) = Q_3^{-1},
   and fit consecutive logarithmic rates:
       sigma_{H21} = Delta log(R_{21}) / Delta N,
       sigma_{H31} = Delta log(R_{31}) / Delta N.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 73 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix
from cell import get_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]
BARRIER_V_STAR = mp.mpf('1.0')


# -----------------------------------------------------------------------------
# Parity Basis & Operator Construction
# -----------------------------------------------------------------------------

def full_parity_basis(N: int) -> tuple[mp.matrix, mp.matrix]:
    """Construct orthonormal basis matrices E (even) and O (odd) for R^{2N+1}."""
    dim = 2 * N + 1
    E = mp.matrix(dim, N + 1)
    O = mp.matrix(dim, N)

    E[N, 0] = mp.mpf(1)

    inv_sqrt2 = 1 / mp.sqrt(2)
    for m in range(1, N + 1):
        E[N + m, m] = inv_sqrt2
        E[N - m, m] = inv_sqrt2
        O[N + m, m - 1] = inv_sqrt2
        O[N - m, m - 1] = -inv_sqrt2

    return E, O


def solve_parity_eigensystems(
    Q_full: mp.matrix, N: int
) -> tuple[
    mp.mpf, mp.matrix, mp.matrix,
    list[mp.mpf], mp.matrix,
    list[mp.mpf], mp.matrix
]:
    """Compute even and odd spectra and eigenvectors."""
    E, O = full_parity_basis(N)

    Q_even = E.T * Q_full * E
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)

    Q_odd = O.T * Q_full * O
    Q_odd = mp.mpf('0.5') * (Q_odd + Q_odd.T)

    evals_e, V_e = mp.eigsy(Q_even)
    evals_o, V_o = mp.eigsy(Q_odd)

    # Sort even
    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    # Sort odd
    idx_o = sorted(range(N), key=lambda i: evals_o[i])
    sorted_evals_o = [evals_o[i] for i in idx_o]
    sorted_V_o = mp.matrix(N, N)
    for col_idx, orig_col in enumerate(idx_o):
        for row_idx in range(N):
            sorted_V_o[row_idx, col_idx] = V_o[row_idx, orig_col]

    lam_0 = sorted_evals_e[0]
    return lam_0, E, O, sorted_evals_e, sorted_V_e, sorted_evals_o, sorted_V_o


# -----------------------------------------------------------------------------
# Regularized Stieltjes Decomposition Function
# -----------------------------------------------------------------------------

def compute_H_decomposition(
    mu: mp.mpf,
    lam_0: mp.mpf,
    D_0_sq: mp.mpf,
    d_k_list: list[mp.mpf],
    evals_e: list[mp.mpf],
    N: int,
) -> dict:
    """
    Compute term-by-term pole distribution of H(mu):
        H(mu) = D_0^2 + sum_{k=1}^N d_k^2 * ((mu - lam_0) / (E_k - mu))^2
    """
    gap = mu - lam_0
    gap_sq = gap ** 2

    # Term k = 0
    H_k_terms: list[mp.mpf] = [D_0_sq]

    # Terms k = 1, ..., N
    for k in range(1, N + 1):
        E_k = evals_e[k]
        denom = (E_k - mu) ** 2
        term_k = (d_k_list[k] ** 2) * (gap_sq / denom)
        H_k_terms.append(term_k)

    H_sum = sum(H_k_terms)

    # Cross-check via Stieltjes derivative definition
    # G_d'(mu) = sum_{k=0}^N d_k^2 / (E_k - mu)^2
    G_prime = sum((d_k_list[k] ** 2) / ((evals_e[k] - mu) ** 2) for k in range(N + 1))
    H_stieltjes = gap_sq * G_prime
    id_residual = abs(H_sum - H_stieltjes)

    # Percentage shares
    shares = [(H_k_terms[k] / H_sum) * 100 for k in range(N + 1)]

    # Dominant pole (over all k in 0..N)
    k_dom = max(range(N + 1), key=lambda k: H_k_terms[k])
    H_dom = H_k_terms[k_dom]
    share_dom = shares[k_dom]

    # Nearest even pole (over excited poles k in 1..N)
    k_near = min(range(1, N + 1), key=lambda k: abs(evals_e[k] - mu))
    dist_near = abs(evals_e[k_near] - mu)
    H_near = H_k_terms[k_near]
    share_near = shares[k_near]
    bound_ratio = H_sum / H_near if H_near > 0 else mp.inf

    return {
        "H_sum": H_sum,
        "H_stieltjes": H_stieltjes,
        "id_residual": id_residual,
        "H_k_terms": H_k_terms,
        "shares": shares,
        "k_dom": k_dom,
        "H_dom": H_dom,
        "share_dom": share_dom,
        "k_near": k_near,
        "dist_near": dist_near,
        "H_near": H_near,
        "share_near": share_near,
        "bound_ratio": bound_ratio,
    }


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell73() -> None:
    print("=" * 105)
    print("CELL 73 — AUDIT OF POSITIVE REGULARIZED STIELTJES FUNCTION H(mu),")
    print("         EXACT MODAL RATIO Q_j = H(mu_1)/H(mu_j), AND INTERLACING POLE SPECTRUM")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
    print()

    synthesis_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Retrieve cached Galerkin matrix
        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

        # Solve parity eigensystems
        lam_0, E, O, evals_e, V_e, evals_o, V_o = solve_parity_eigensystems(Q_full, N)

        # Ground state c in R^{N+1}
        c_even = V_e[:, 0]

        # Boundary vector d in R^{N+1}: d_even = (1, sqrt(2), ..., sqrt(2))^T
        d_even = mp.matrix(N + 1, 1)
        d_even[0, 0] = mp.mpf(1)
        for m in range(1, N + 1):
            d_even[m, 0] = mp.sqrt(2)

        D_0 = sum(d_even[m, 0] * c_even[m, 0] for m in range(N + 1))
        D_0_sq = D_0 ** 2

        # Boundary overlaps in even sector: d_k = <u_k^{even}, d>
        d_k_list: list[mp.mpf] = []
        for k in range(N + 1):
            d_val = sum(d_even[m, 0] * V_e[m, k] for m in range(N + 1))
            d_k_list.append(d_val)

        # Coordinate transition dipoles b_{0j} = <c, K u_j>
        b_0j_list: list[mp.mpf] = []
        for j in range(N):
            b_j = sum(c_even[m, 0] * m * V_o[m - 1, j] for m in range(1, N + 1))
            b_0j_list.append(b_j)

        # Coordinate kinetic energies: ||K u_j||^2 = sum_{m=1}^N m^2 (u_{j, m-1})^2
        ke_list: list[mp.mpf] = []
        for j in range(N):
            ke_val = sum((m ** 2) * (V_o[m - 1, j] ** 2) for m in range(1, N + 1))
            ke_list.append(ke_val)

        ke_1 = ke_list[1]
        b_01_sq = b_0j_list[1] ** 2

        # Compute H decomposition for j = 1, 2, 3
        h_decomp: dict[int, dict] = {}
        for j in [1, 2, 3]:
            if j < N:
                mu_j = evals_o[j]
                h_decomp[j] = compute_H_decomposition(mu_j, lam_0, D_0_sq, d_k_list, evals_e, N)

        H_1 = h_decomp[1]["H_sum"]

        # Audit mode 2 and mode 3
        mode_data: dict[int, dict] = {}
        for j in [2, 3]:
            if j < N:
                H_j = h_decomp[j]["H_sum"]
                Q_j_H = H_1 / H_j
                R_j1 = H_j / H_1

                ke_j = ke_list[j]
                K_j = ke_j / ke_1
                T_j_direct = (b_0j_list[j] ** 2) / b_01_sq
                T_j_H = K_j * Q_j_H
                tail_res = abs(T_j_direct - T_j_H)

                mode_data[j] = {
                    "H_j": H_j,
                    "Q_j_H": Q_j_H,
                    "R_j1": R_j1,
                    "K_j": K_j,
                    "T_j_direct": T_j_direct,
                    "T_j_H": T_j_H,
                    "tail_res": tail_res,
                }

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "lam_0": lam_0,
            "D_0_sq": D_0_sq,
            "evals_e": evals_e,
            "evals_o": evals_o,
            "d_k_list": d_k_list,
            "h_decomp": h_decomp,
            "mode_data": mode_data,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  H(mu_1)                      = {mp.nstr(H_1, 9)}")
        print(f"  H(mu_1) Stieltjes residual   = {mp.nstr(h_decomp[1]['id_residual'], 4)}")
        if 2 in mode_data:
            m2 = mode_data[2]
            hd2 = h_decomp[2]
            print(f"  H(mu_2)                      = {mp.nstr(m2['H_j'], 9)}")
            print(f"  H(mu_2) Stieltjes residual   = {mp.nstr(hd2['id_residual'], 4)}")
            print(f"  Ratio Q_2 = H(mu_1)/H(mu_2)  = {mp.nstr(m2['Q_j_H'], 9)}")
            print(f"  Inv ratio H(mu_2)/H(mu_1)    = {mp.nstr(m2['R_j1'], 9)}")
            print(f"  Tail residual |T_2 - K_2*Q_2|= {mp.nstr(m2['tail_res'], 4)}")
            print(f"  Dominant pole k_dom          = {hd2['k_dom']} (Share: {float(hd2['share_dom']):.2f}%)")
            print(f"  Nearest pole k_near          = {hd2['k_near']} (dist = {mp.nstr(hd2['dist_near'], 5)}, bound ratio = {float(hd2['bound_ratio']):.4f})")
        if 3 in mode_data:
            m3 = mode_data[3]
            hd3 = h_decomp[3]
            print(f"  H(mu_3)                      = {mp.nstr(m3['H_j'], 9)}")
            print(f"  Ratio Q_3 = H(mu_1)/H(mu_3)  = {mp.nstr(m3['Q_j_H'], 9)}")
            print(f"  Inv ratio H(mu_3)/H(mu_1)    = {mp.nstr(m3['R_j1'], 9)}")
            print(f"  Dominant pole k_dom          = {hd3['k_dom']} (Share: {float(hd3['share_dom']):.2f}%)")
            print(f"  Nearest pole k_near          = {hd3['k_near']} (dist = {mp.nstr(hd3['dist_near'], 5)}, bound ratio = {float(hd3['bound_ratio']):.4f})")
        print()

    # -------------------------------------------------------------------------
    # Synthesis Tables
    # -------------------------------------------------------------------------

    print("=" * 105)
    print("SYNTHESIS TABLE 1: POSITIVE REGULARIZED STIELTJES FUNCTION H(mu) AND EXACT MODAL RATIOS")
    print("=" * 105)
    header1 = (
        f"{'N':>3} | {'H(mu_1)':>15} | {'H(mu_2)':>15} | {'Q_2 = H1/H2':>15} | "
        f"{'H(mu_2)/H(mu_1)':>15} | {'Q_3 = H1/H3':>15} | {'Max Id Res':>10}"
    )
    print(header1)
    print("-" * len(header1))
    for rec in synthesis_records:
        N = rec["N"]
        H_1 = rec["h_decomp"][1]["H_sum"]
        m2 = rec["mode_data"][2]
        m3 = rec["mode_data"][3]
        max_res = max(
            rec["h_decomp"][1]["id_residual"],
            rec["h_decomp"][2]["id_residual"],
            rec["h_decomp"][3]["id_residual"],
            m2["tail_res"],
            m3["tail_res"],
        )
        print(
            f"{N:3d} | {mp.nstr(H_1, 9):>15} | {mp.nstr(m2['H_j'], 9):>15} | "
            f"{mp.nstr(m2['Q_j_H'], 9):>15} | {mp.nstr(m2['R_j1'], 9):>15} | "
            f"{mp.nstr(m3['Q_j_H'], 9):>15} | {mp.nstr(max_res, 3):>10}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: TERM-BY-TERM POLE DECOMPOSITION OF H(mu_2)")
    print("=" * 105)
    header2 = (
        f"{'N':>3} | {'H_0 = D_0^2':>13} | {'H_1 (pole 1)':>14} | {'H_2 (pole 2)':>14} | "
        f"{'H_3 (pole 3)':>14} | {'H_4 (pole 4)':>14} | {'H_{>=5}':>14} | {'H_total(mu_2)':>15}"
    )
    print(header2)
    print("-" * len(header2))
    for rec in synthesis_records:
        N = rec["N"]
        hd2 = rec["h_decomp"][2]
        H_terms = hd2["H_k_terms"]
        H_0 = H_terms[0]
        H_1 = H_terms[1] if len(H_terms) > 1 else mp.mpf(0)
        H_2 = H_terms[2] if len(H_terms) > 2 else mp.mpf(0)
        H_3 = H_terms[3] if len(H_terms) > 3 else mp.mpf(0)
        H_4 = H_terms[4] if len(H_terms) > 4 else mp.mpf(0)
        H_ge5 = sum(H_terms[k] for k in range(5, len(H_terms))) if len(H_terms) > 5 else mp.mpf(0)
        H_tot = hd2["H_sum"]
        print(
            f"{N:3d} | {mp.nstr(H_0, 7):>13} | {mp.nstr(H_1, 7):>14} | {mp.nstr(H_2, 7):>14} | "
            f"{mp.nstr(H_3, 7):>14} | {mp.nstr(H_4, 7):>14} | {mp.nstr(H_ge5, 7):>14} | {mp.nstr(H_tot, 9):>15}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: POLE SHARES (%), DOMINANCE, AND NEAREST-POLE BOUND FOR MODE j = 2")
    print("=" * 105)
    header3 = (
        f"{'N':>3} | {'k_dom':>5} | {'Share_dom(%)':>12} | {'k_near':>6} | "
        f"{'|E_near - mu_2|':>17} | {'H_{k_near}':>15} | {'Bound Ratio H/H_near':>20}"
    )
    print(header3)
    print("-" * len(header3))
    for rec in synthesis_records:
        N = rec["N"]
        hd2 = rec["h_decomp"][2]
        print(
            f"{N:3d} | {hd2['k_dom']:5d} | {float(hd2['share_dom']):12.4f} | {hd2['k_near']:6d} | "
            f"{mp.nstr(hd2['dist_near'], 9):>17} | {mp.nstr(hd2['H_near'], 9):>15} | {float(hd2['bound_ratio']):20.4f}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: COMPARATIVE POLE ARCHITECTURE FOR MODE j = 3")
    print("=" * 105)
    header4 = (
        f"{'N':>3} | {'k_dom':>5} | {'Share_dom(%)':>12} | {'k_near':>6} | "
        f"{'|E_near - mu_3|':>17} | {'H_{k_near}':>15} | {'Bound Ratio H/H_near':>20}"
    )
    print(header4)
    print("-" * len(header4))
    for rec in synthesis_records:
        N = rec["N"]
        hd3 = rec["h_decomp"][3]
        print(
            f"{N:3d} | {hd3['k_dom']:5d} | {float(hd3['share_dom']):12.4f} | {hd3['k_near']:6d} | "
            f"{mp.nstr(hd3['dist_near'], 9):>17} | {mp.nstr(hd3['H_near'], 9):>15} | {float(hd3['bound_ratio']):20.4f}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 5: CONSECUTIVE LOGARITHMIC GROWTH RATES FOR INVERSE RATIOS R_{j1} = H(mu_j)/H(mu_1)")
    print("=" * 105)
    header5 = (
        f"{'Interval':>14} | {'Delta N':>7} | {'R_21(start)':>15} | {'R_21(end)':>15} | "
        f"{'sigma_H21':>12} | {'sigma_H31':>12}"
    )
    print(header5)
    print("-" * len(header5))
    for idx in range(len(synthesis_records) - 1):
        r_a = synthesis_records[idx]
        r_b = synthesis_records[idx + 1]
        n_a, n_b = r_a["N"], r_b["N"]
        dn = n_b - n_a

        r21_a = r_a["mode_data"][2]["R_j1"]
        r21_b = r_b["mode_data"][2]["R_j1"]
        sigma_h21 = (mp.log(r21_b) - mp.log(r21_a)) / dn

        r31_a = r_a["mode_data"][3]["R_j1"]
        r31_b = r_b["mode_data"][3]["R_j1"]
        sigma_h31 = (mp.log(r31_b) - mp.log(r31_a)) / dn

        interval_str = f"N = {n_a:2d} -> {n_b:2d}"
        print(
            f"{interval_str:>14} | {dn:7d} | {mp.nstr(r21_a, 9):>15} | {mp.nstr(r21_b, 9):>15} | "
            f"{float(sigma_h21):12.5f} | {float(sigma_h31):12.5f}"
        )
    print()

    print("=" * 105)
    print("CELL 73 EXECUTION COMPLETE")
    print("=" * 105)


if __name__ == "__main__":
    run_cell73()
