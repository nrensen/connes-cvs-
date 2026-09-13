# Cell 123 Analytical Note: Operator Decomposition Audit ($Q_{\mathrm{arch}}, Q_{\mathrm{prime}}, Q_{\mathrm{pole}}$) & High-Mode Coercivity

**Companion Document:** Analytical Research Note (Gate 1 / Milestone M-G1.4)  
**Status:** Pre-Flight Research Note / Theoretical Reduction & Diagnostic Hypotheses  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.4: Continuum Spectral Threshold Proof](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell122.md`](file:///c:/data/github/connes-cvs-/cell122.md); [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md); [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out); [`cell120.out`](file:///c:/data/github/connes-cvs-/cell120.out); [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py)  
**Date:** September 2026  

---

## 1. Executive Summary & Audit Mandate

In Cell 122, an analytical mechanism was proposed to secure the continuum spectral lower bound $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0$ for core sizes $L \ge 12$. The goal is to eliminate spectral crowding in the remote denominator $R_{\mathrm{spec}}(N, L)$ and isolate Gate 1 to the asymptotic decay of the parity-doublet tunneling splitting $\Delta_j(N)$.

However, rigorous analysis of the proposal revealed that the original premise of Route C ($C_3 \succeq 0.386 I$ on $\mathcal{H}_{\mathrm{high}} = \operatorname{span}\{e_3, \dots, e_N\}$) was blocked by two distinct mathematical obstructions:
1. **The Multiplier-to-Matrix Representation Distinction:** In [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py), $Q_{\mathrm{arch}}$ is a dense divided-difference matrix, not a diagonal matrix $\operatorname{diag}(h_+(a_m))$. Furthermore, $Q_{\mathrm{prime}}$ carries an explicit minus sign and has a negative expectation value on the ground state ($Q_{\mathrm{prime}}[u_0] \approx -0.0775$). Its sign on high modes is unproven.
2. **The Cauchy Interlacing Bound-State Obstruction:** Cauchy's interlacing theorem for deleting $M = 3$ modes rigorously implies:
   $$\lambda_{\min}(C_3) = \lambda_0(C_3) \le E_3(Q_{\mathrm{even}}).$$
   Because the potential well carries $\bar{N}_{\mathrm{bound}} \approx 11$ bound states clustered near zero, $E_3(64) \approx 2.98 \times 10^{-38}$. Thus $\lambda_{\min}(C_3) \le 10^{-38} \ll 0.386$. The submatrix $C_3$ cannot be coercive with constant $0.386$ because $\sim 8$ bound states of the well remain supported inside $C_3$.

### The Audit Mandate for Cell 123
Cell 123 executes an exact, high-precision operator decomposition audit directly from the code definitions to:
- Dissect $Q_{\mathrm{even}} = Q_{\mathrm{even}, \mathrm{arch}} + Q_{\mathrm{even}, \mathrm{prime}} + Q_{\mathrm{even}, \mathrm{pole}}$ into its constituent symmetric components.
- Numerically test the $M = 3$ block: determine whether $\lambda_{\min}(C_3) \le E_3$ holds, measure the component spectra, and evaluate the interior eigenvalue $\lambda_{10}(C_3)$ against $E_{13}$.
- Test the $M = 12$ continuum block ($C_{12}$ on $\operatorname{span}\{e_{12}, \dots, e_N\}$): determine whether removing the entire 12-dimensional bound-state core restores macroscopic coercivity $\lambda_{\min}(C_{12}) \ge c_{12} > 0$ ($c_{12} \approx 0.50$).
- Sweep the sub-block cut $M \in \{3, 4, 6, 8, 10, 12, 14, 16\}$ to directly observe the phase transition of $\lambda_{\min}(C_M)$ from bound-state collapse ($M \le 10$) to macroscopic continuum coercivity ($M \ge 12$).

---

## 2. Tripartite Galerkin Decomposition in Canonical Coordinates

### 2.1 Exact Operator Definitions
In [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py), the symbol $\psi(x)$ is the sum of three odd functions:
$$\psi(x) = \psi_{\mathrm{prime}}(x) + \psi_{\mathrm{pole}}(x) + \psi_{\mathrm{arch}}(x),$$
defined for cutoff $c = 13$ and $L = \log 13$ by:
1. **Prime Piece:**
   $$\psi_{\mathrm{prime}}(x) = -\frac{1}{\pi} \sum_{n \text{ prime power} \le c} \frac{\Lambda(n)}{\sqrt{n}} \sin\left(2\pi x \left(1 - \frac{\log n}{L}\right)\right).$$
2. **Pole Piece:**
   $$\psi_{\mathrm{pole}}(x) = \frac{1}{\pi} \int_0^L \sin\left(2\pi x \left(1 - \frac{y}{L}\right)\right) \cdot 2\cosh\left(\frac{y}{2}\right) dy.$$
3. **Archimedean Piece:**
   $$\psi_{\mathrm{arch}}(x) = \frac{1}{2\pi^2} \int_{-T}^T h_+(\tau) \operatorname{Re}(\hat{S}(\tau, x)) d\tau.$$

Because the divided-difference operation is strictly linear:
$$\frac{\psi(m) - \psi(n)}{m - n} = \frac{\psi_{\mathrm{prime}}(m) - \psi_{\mathrm{prime}}(n)}{m - n} + \frac{\psi_{\mathrm{pole}}(m) - \psi_{\mathrm{pole}}(n)}{m - n} + \frac{\psi_{\mathrm{arch}}(m) - \psi_{\mathrm{arch}}(n)}{m - n},$$
and on the diagonal $\psi'(n) = \psi_{\mathrm{prime}}'(n) + \psi_{\mathrm{pole}}'(n) + \psi_{\mathrm{arch}}'(n)$, the full $(2N+1) \times (2N+1)$ Galerkin matrix decomposes identically as:
$$Q_{\mathrm{total}} = Q_{\mathrm{prime}} + Q_{\mathrm{pole}} + Q_{\mathrm{arch}}.$$

### 2.2 Parity Projection onto the Even Sector
All three components $\psi_{\mathrm{prime}}, \psi_{\mathrm{pole}}, \psi_{\mathrm{arch}}$ are strictly odd functions: $\psi(-x) = -\psi(x)$.
Consequently, each matrix $Q_{\mathrm{comp}} \in \{Q_{\mathrm{prime}}, Q_{\mathrm{pole}}, Q_{\mathrm{arch}}\}$ satisfies:
$$Q_{\mathrm{comp}}[-m, -n] = Q_{\mathrm{comp}}[m, n] \qquad \text{(centrosymmetry)}.$$
Under the canonical even-basis isometry $V_{\mathrm{even}}: \mathbb{R}^{N+1} \to \mathbb{R}^{2N+1}$:
$$V_{\mathrm{even}}[N, 0] = 1, \qquad V_{\mathrm{even}}[N \pm m, m] = \frac{1}{\sqrt{2}} \quad (m \ge 1),$$
the projected even-sector matrix splits linearly without cross-terms:
$$Q_{\mathrm{even}} = V_{\mathrm{even}}^T Q_{\mathrm{total}} V_{\mathrm{even}} = Q_{\mathrm{even}, \mathrm{prime}} + Q_{\mathrm{even}, \mathrm{pole}} + Q_{\mathrm{even}, \mathrm{arch}}.$$

---

## 3. The Dual Submatrix Hypotheses: $M = 3$ vs $M = 12$

For any cutoff index $M \in \{1, \dots, N\}$, let $\mathcal{H}_{\mathrm{high}}^{(M)} = \operatorname{span}\{e_M, \dots, e_N\}$ denote the subspace spanned by modes $m \ge M$. The corresponding principal submatrix of $Q_{\mathrm{even}}$ is:
$$C_M = Q_{\mathrm{even}}[M..N, M..N] \in \mathbb{R}^{(N-M+1) \times (N-M+1)}.$$

Cauchy's interlacing theorem for removing $M$ rows and columns from an $(N+1) \times (N+1)$ symmetric matrix states:
$$E_k^{(N)} \le \lambda_k(C_M) \le E_{k+M}^{(N)} \qquad (\forall k \in \{0, \dots, N-M\}).$$

### 3.1 The $M = 3$ Submatrix (Multiplier Sign-Change Boundary)
At $c = 13$, the discrete Fourier frequencies $a_m = \frac{2\pi m}{\log 13}$ cross the Archimedean sign-change point $r_* \approx 6.28984$ between $m = 2$ and $m = 3$:
$$h_+(a_2) \approx -0.914 < 0, \qquad h_+(a_3) \approx +0.386 > 0.$$

- **Hypothesis H-C3.1 (Bound-State Collapse of $\lambda_{\min}(C_3)$):**
  Because $M = 3$ removes only 3 modes from an $\sim 11$-dimensional bound-state well, the remaining $\sim 8$ bound states reside inside $C_3$. By Cauchy interlacing with $k = 0$:
  $$\lambda_{\min}(C_3) \le E_3^{(N)} \approx 2.98 \times 10^{-38} \ll 0.386.$$
  $C_3$ is **not coercive** with constant $0.386$.
- **Hypothesis H-C3.2 (Interior Interlacing Lower Bound):**
  By Cauchy interlacing with $k = 10$:
  $$E_{13}^{(N)} \ge \lambda_{10}(C_3) \ge E_{10}^{(N)}.$$
  The 11th eigenvalue $\lambda_{10}(C_3)$ lies in the scattering continuum ($\ge 0.40$), while the ground state $\lambda_0(C_3)$ remains collapsed.

### 3.2 The $M = 12$ Submatrix (Bound-State Cluster Boundary)
Cell 121 demonstrated that the bound-state cluster saturates at $L \approx 11$. For $L \ge 12$, the Ritz gap $g_{2, L}(N)$ stabilizes to macroscopic values ($g_{2, 12} \ge 0.582$).

- **Hypothesis H-C12 (Continuum Submatrix Coercivity):**
  Projecting out all $M = 12$ modes ($m \in \{0, \dots, 11\}$) eliminates the entire bound-state cluster. The principal submatrix $C_{12} = Q_{\mathrm{even}}[12..N, 12..N]$ satisfies:
  $$C_{12} \succeq c_{12} I \qquad \text{with } c_{12} \approx 0.50 > 0 \quad \text{uniformly across } N.$$
- **Cauchy Consequence for Gate 1:**
  By Cauchy interlacing with $M = 12$ and $k = 1$:
  $$\boxed{E_{13}^{(N)} \ge \lambda_1(C_{12}) \ge \lambda_{\min}(C_{12}) \ge c_{12} > 0 \qquad (\forall N \ge 13).}$$
  Because $E_{L+1}^{(N)} \ge E_{13}^{(N)}$ for all $L \ge 12$, establishing $C_{12} \succeq c_{12} I > 0$ **unconditionally proves the Denominator Theorem** for all $L \ge 12$:
  $$E_{L+1}^{(N)} - E_{j+1}^{(N)} \ge c_{12} - E_{j+1}^{(N)} \ge \delta > 0.$$

---

## 4. Component Sign Structures & Audit Questions

For both $M = 3$ and $M = 12$, Cell 123 evaluates the tripartite components:
$$C_M = C_{M, \mathrm{arch}} + C_{M, \mathrm{prime}} + C_{M, \mathrm{pole}}.$$

| Component | Construction in Code | Analytical Sign Expectation on High Modes | Open Audit Question |
| :--- | :--- | :--- | :--- |
| **Archimedean ($C_{\mathrm{arch}}$)** | Smeared Fejér integral against $h_+(\tau)$ | $h_+(a_m) > 0$ for $m \ge 3$, but sinc tails overlap $[0, r_*]$ | Does off-diagonal leakage allow $\lambda_{\min}(C_{\mathrm{arch}}) < 0$ on $M = 3$? Does it become strictly positive on $M = 12$? |
| **Prime ($C_{\mathrm{prime}}$)** | Finite cosine/sine sum over $p^k \le 13$ with explicit minus sign | Negative on ground state ($Q_{\mathrm{prime}}[u_0] \approx -0.078$) | Is $C_{\mathrm{prime}}$ negative, positive, or indefinite on $M = 3$ and $M = 12$? What is its spectral norm relative to $C_{\mathrm{arch}}$? |
| **Pole ($C_{\mathrm{pole}}$)** | Integral of $2\cosh(y/2) \sin(\dots)$ | Positive definite on full space ($Q_{\mathrm{pole}}[u_0] \approx +1.572$) | Does $C_{\mathrm{pole}}$ remain strictly positive definite on high-mode submatrices? |

---

## 5. Quantitative Falsification & Verification Criteria

1. **Falsification of Route C $M = 3$ Coercivity:**
   If $\lambda_{\min}(C_3) \le 10^{-20}$ across tested dimensions $N \in [16, 64]$, the naive Route C coercivity premise ($C_3 \succeq 0.386 I$) is **conclusively falsified**.
2. **Verification of Cauchy Interlacing Bounds:**
   - For $M = 3$: Confirm $\lambda_{\min}(C_3) \le E_3^{(N)}$ and $E_{13}^{(N)} \ge \lambda_{10}(C_3)$ to machine precision.
   - For $M = 12$: Confirm $E_{13}^{(N)} \ge \lambda_1(C_{12}) \ge \lambda_{\min}(C_{12})$ to machine precision.
3. **Certification of Core-Submatrix Coercivity ($M = 12$):**
   If $\lambda_{\min}(C_{12}) \ge c_{12} > 0$ with $c_{12} \approx 0.50$ across all tested dimensions $N \in [16, 64]$, the **Core-Submatrix Coercivity Conjecture** is empirically secured, providing the exact target for analytical certification.
4. **Sharpness of the Bound-State-to-Continuum Transition:**
   Across the sweep $M \in \{3, 4, 6, 8, 10, 12, 14, 16\}$, $\lambda_{\min}(C_M)$ must exhibit a sharp transition from exponential collapse ($M \le 10$) to macroscopic positivity ($M \ge 12$), mirroring the Ritz gap transition from Cell 121.
