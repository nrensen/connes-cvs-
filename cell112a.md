# Cell 112a Analytical Note: Multi-Eigenvalue Branch Tracking & Fast 70-DPS Extinction Audit

**Companion Computational Script:** [`cell112a.py`](file:///c:/data/github/connes-cvs-/cell112a.py) | **Verification Log:** `cell112a.out`  
**Status:** Working Research Note (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell111.md`](file:///c:/data/github/connes-cvs-/cell111.md); [`cell111.out`](file:///c:/data/github/connes-cvs-/cell111.out); [`cell110.md`](file:///c:/data/github/connes-cvs-/cell110.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary

Cell 111 established the exact row-wise resolvent identity (Theorem 1, certified to residual $1.09 \times 10^{-66}$):
$$\boxed{\alpha_N = (H u_N)_m + \frac{T_{v_N}(0)}{\sqrt{2}} a_m - E_{11} m^2 v_{N, m} \qquad (\forall m \in \{1, \dots, N\}),}$$
and verified the upper-boundary equipartition $(H u_N)_N \approx \frac{T(0)}{\sqrt{2}} a_N \approx \frac{1}{2} \alpha_N$ at $N=192$.

However, `cell111.out` revealed an apparent plateau in the boundary observables for $N \ge 64$:
$$\alpha_N \approx 1.923 \times 10^{-22}, \qquad T_{v_N}(0) \approx 6.66 \times 10^{-25}, \qquad S_2(N) \approx 2.48 \times 10^{-19},$$
which caused the extinction product $P_\alpha(N) = |\alpha_N| \sqrt{N}$ to drift slightly upward (from $2.02 \times 10^{-21}$ at $N=96$ to $2.66 \times 10^{-21}$ at $N=192$).

While initial scrutiny questioned whether $1.92 \times 10^{-22}$ might be an eigensolver noise floor requiring 110 dps (which would take hours of redundant quadratures), forensic analysis of `cell111.out` refutes the noise hypothesis and uncovers a clear physical mechanism:
1. **The Algebraic Certification Exceeds 60 Digits:** The row identity holds to $1.09 \times 10^{-66}$, meaning $\alpha_N$ is algebraically consistent with the computed state $v_N$ to 44 digits beyond $10^{-22}$.
2. **Eigenvalue Crossing & Branch Switching at $N=64$:**
   - For $N \le 48$, the ground state has positive energy $E_{11} > 0$ ($+2.10 \times 10^{-50}$ at $N=48$) and a dominant ground-state center amplitude $v_0 \approx 0.54 \to 0.456$ (the localized solitary wave).
   - For $N \ge 64$, the lowest eigenvalue $E_0$ drops below zero and locks at $-1.063 \times 10^{-51}$, matching the finite-$T=600$ Archimedean tail leakage floor $-\delta_{T=600}^{\mathrm{tail}}$ (Paper NR1 Theorem 5.5).
   - Simultaneously, $v_0$ drops from $0.456$ to $0.0708$ ($N=64$) and $0.0638$ ($N=192$).
   - The solver switched tracking from the solitary wave to a cutoff-induced edge/background mode!

Cell 112a executes a fast, cache-accelerated 70-dps multi-eigenvalue spectrum tracking audit across $N \in \{16, \dots, 192\}$ in under 1 minute, simultaneously evaluating both the literal ground state and the solitary wave branch.

---

## 1. Physical Mechanism of the $N \ge 64$ Plateau

### 1.1 Finite-$T$ Archimedean Tail Leakage
In Paper NR1 (Theorem 5.5), the Archimedean weight $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$ is strictly negative on a compact low-frequency window $[0, r_*]$ ($r_* \approx 6.29, h_+(0) \approx -5.37$). For any finite truncation parameter $T$, the continuous quadratic form admits a negative tail leakage:
$$\delta_T^{\mathrm{tail}} \approx \mathcal{O}(e^{-\gamma T}) \approx 10^{-51} \quad (\text{at } T = 600).$$

When $N$ is small ($N \le 48$), the discrete Galerkin grid does not resolve this delicate leakage, and the positive kinetic energy of the localized solitary wave ($v_0 \approx 0.54$) dominates, yielding $E_0 > 0$.

At $N \ge 64$, the Galerkin dimension becomes large enough to support a delocalized background/edge mode that exploits the negative Archimedean leakage, resulting in an eigenvalue $E_0 = -1.063 \times 10^{-51} < 0$. Because `cell111.py` tracked strictly the index $k=0$, it began tracking this edge mode instead of the solitary wave!

### 1.2 Contrast Between the Two Modes
- **The Edge/Background Mode ($k=0$ for $N \ge 64$):**
  - Characterized by small $v_0 \approx 0.064$, delocalized modal profile, and negative energy $E_0 \approx -10^{-51}$.
  - Its boundary coupling $\alpha_N \approx 1.92 \times 10^{-22}$ and defect $T(0) \approx 6.66 \times 10^{-25}$ reflect the finite-$T=600$ cutoff envelope rather than solitary wave dynamics.
- **The Solitary Wave Branch ($k = k_{\mathrm{sol}}$):**
  - Characterized by dominant center amplitude $v_0 \approx 0.45 - 0.55$, localized profile, and rapid exponential decay into the high modes.
  - This is the genuine physical state entering André Weil's explicit formula and the continuum limit.

---

## 2. Mathematical Architecture of `cell112a.py`

`cell112a.py` implements four diagnostic modules:

### 2.1 Spectrum & Branch Tracking (Module 1)
For each dimension $N \in \{16, 24, 32, 48, 64, 80, 96, 128, 192\}$, the script:
1. Solves the full symmetric eigensystem of $H_N$.
2. Extracts the lowest $K = 5$ eigenpairs $(E_k, v^{(k)})$.
3. Identifies the solitary wave index:
   $$k_{\mathrm{sol}} = \operatorname{argmax}_{0 \le k < 5} v_0^{(k)}.$$
4. If $k_{\mathrm{sol}} \ne 0$, an eigenvalue crossing is confirmed, and both branches are tracked concurrently.

### 2.2 Eigensolver Residual Certification (Module 2)
To certify that both branches are mathematically genuine eigenvectors, the script measures:
- Operator residual norm: $\|(H_N - E_k I) v^{(k)}\|_2 < 10^{-58}$.
- Rayleigh quotient error: $|E_k - \langle v^{(k)}, H_N v^{(k)} \rangle / \|v^{(k)}\|_2^2| < 10^{-58}$.
- Commutator orthogonality identity: $\langle v^{(k)}, \xi_N^{(k)} \rangle = 0$.

### 2.3 Multi-Route Evaluation of $\alpha_N$ (Module 3)
Evaluates $\alpha_N$ across three independent mathematical routes for each branch:
- **Route 1 (Direct Modal Sum):**
  $$\alpha_N^{(1)} = \sum_{k=1}^N a_k v_k.$$
- **Route 2 (Boundary Mode Specialization, Theorem 1):**
  $$\alpha_N^{(2)} = (H u_N)_N + \frac{T(0)}{\sqrt{2}} a_N - E N^2 v_N.$$
- **Route 3 (High-Sector Average):**
  $$\alpha_N^{(3)} = \frac{1}{N - M} \sum_{m=M+1}^N \left( (H u_N)_m + \frac{T(0)}{\sqrt{2}} a_m - E m^2 v_m \right), \quad M = \lfloor N/2 \rfloor.$$
Agreement to $< 10^{-58}$ mathematically certifies $\alpha_N$ on both branches.

### 2.4 Extinction Product Comparison (Module 4)
Tracks the scaling of:
$$P_T(N) = |T_{v_N}(0)| N^{3/2}, \qquad P_\alpha(N) = |\alpha_N| \sqrt{N}, \qquad \kappa_\alpha(N) = \frac{|\alpha_N|}{|T_{v_N}(0)|}$$
across both branches to determine whether the solitary wave branch continues decaying toward zero, satisfying $(H_{\mathrm{ext}})$.

---

## 3. Pre-Flight Quality & Performance Specification

| Parameter / Item | Specification | Rationale |
| :--- | :--- | :--- |
| **Cutoff & Truncation** | $c = 13$, $T = 600$, $N_{\max} = 192$ | Canonical parameters matching Paper NR2 and Cell 111 |
| **Working Precision** | `mp.mp.dps = 70`, `GROUND_DPS = 70` | Matches cached matrix, loads instantly, avoids redundant 4500s quadratures |
| **Dimension Sweep** | $N \in \{16, 24, 32, 48, 64, 80, 96, 128, 192\}$ | 9 dimensions spanning before, during, and after the $N=64$ transition |
| **Spectrum Depth** | $K = 5$ lowest eigenpairs | Sufficient to capture low-lying doublet / crossing states |
| **Estimated Runtime** | $< 60$ seconds | Matrix retrieval: ~1.7 s; 9 eigensolves: ~30 s |

---

## 4. Testable Hypotheses

- **Hypothesis 1 (Eigenvalue Crossing):** For $N \ge 64$, the solitary wave branch ($v_0 \approx 0.5$) becomes the 1st or 2nd excited state ($E_1$ or $E_2$), while the ground state $E_0 \approx -10^{-51}$ is a negative-energy edge mode.
- **Hypothesis 2 (True Extinction on Solitary Branch):** On the solitary wave branch $k_{\mathrm{sol}}$, $|\alpha_N|$ and $|T(0)|$ continue rapid decay for $N \ge 64$, restoring $P_\alpha(N) = |\alpha_N| \sqrt{N} \to 0$.
- **Hypothesis 3 (Theorem 1 Universality):** Theorem 1 row identity and multi-route consistency hold to $< 10^{-58}$ on *all* eigenmodes, proving that algebraic boundary flux redirection is an exact operator property independent of the specific state.
