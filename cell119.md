# Cell 119 Analytical Note: Exact Algebraic Modal Decomposition ($A_{N, k}$ vs $B_{N, k}$ Cancellation Balance)

**Companion Computational Script:** [`cell119.py`](file:///c:/data/github/connes-cvs-/cell119.py) | **Verification Log:** [`cell119.out`](file:///c:/data/github/connes-cvs-/cell119.out) (runtime: 1940.57 s at 70 dps)  
**Status:** Executed & Audited / Definitive Capstone of Phase IX (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell118.md`](file:///c:/data/github/connes-cvs-/cell118.md); [`cell118.out`](file:///c:/data/github/connes-cvs-/cell118.out); [`cell117.md`](file:///c:/data/github/connes-cvs-/cell117.md); [`cell116.md`](file:///c:/data/github/connes-cvs-/cell116.md); [`cell112.md`](file:///c:/data/github/connes-cvs-/cell112.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Headline Verdict

Cell 119 was authored as an exact algebraic diagnostic probe to test whether the 20 decades of modal cancellation observed in Cells 117–118 arise from an inter-component cancellation between the boundary symbol term $A_{N, k} = a_N \frac{k^2 v_{N, k}}{N^2 - k^2}$ and the modal arithmetic symbol term $B_{N, k} = a_k \frac{k^2 v_{N, k}}{N^2 - k^2}$.

### Headline Verdict

> **Classification: Decisive Negative Result and Definitive Capstone to Phase IX.**
>
> 1. **Hypothesis H1 (Inter-Sum Cancellation) Is Falsified:**  
>    The $10^{20}$-scale cancellation is **not** an inter-component cancellation between two macroscopic sums $A$ and $B$. At every tested dimension $N \in [64, 192]$, both total off-diagonal sums are already microscopic:
>    $$\Sigma_A^{\mathrm{off}} = \mathcal{O}(10^{-23}), \qquad \Sigma_B^{\mathrm{off}} = \mathcal{O}(10^{-22}).$$
>    At $N = 192$: $\Sigma_A = +4.28 \times 10^{-23}$ and $\Sigma_B = +1.57 \times 10^{-22}$, giving $\Sigma_A - \Sigma_B = -1.14 \times 10^{-22}$ and net flux $(H u_N)_N = -8.35 \times 10^{-23}$. The cancellation factor between the two sums is only $\mathcal{C}_{AB} \approx 1.38$, not $10^{20}$.
>
> 2. **Hypothesis H2 (Naive Moments Pre-Cancellation) Is Falsified:**  
>    Although the solitary moments are small ($M_2 \sim -2.2 \times 10^{-19}$, $M_{\psi, 3} \sim 5.4 \times 10^{-17}$), the exact sums are *two to three orders of magnitude smaller still*. The leading Taylor approximations fail to control the sums, exhibiting large relative errors ($8\times$ to $83\times$), proving that a naive termwise $k/N$ expansion does not uniformly control the oscillatory sum.
>
> 3. **The Real Discovery — Internal Oscillatory Cancellation within the $A$-Trajectory:**  
>    In the cumulative trajectory, $A_{N, k}$ itself undergoes virtually the entire $10^{20}$-scale cancellation:
>    - At $x = 0.01$: $S_A = +3.61 \times 10^{-3}$, while $S_B = -1.52 \times 10^{-6}$.
>    - At $x = 0.02$: $S_A = -1.14 \times 10^{-2}$, while $S_B = +5.67 \times 10^{-5}$.
>    - At $x = 0.10$: $S_A = +4.93 \times 10^{-10}$, while $S_B = -6.10 \times 10^{-12}$.
>    - At $x = 0.20$: $S_A = -2.82 \times 10^{-19}$, while $S_B = +5.58 \times 10^{-20}$.
>    $B$ is comparatively negligible throughout the entire early excursion. The low-mode peak at $k = 3$ is created by $A$, and is subsequent extinguished *internally within the $A$-sum itself across the modal spectrum*.
>
> 4. **Exact Theorem 1 Residual Certified:**  
>    The row-wise identity $(H u_N)_N = \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + E_1 N^2 v_{N, N}$ closes to $2.33 \times 10^{-69}$ residual at $N = 192$ (and $\le 8.9 \times 10^{-70}$ across all smaller $N$). The resolvent term $E_1 N^2 v_{N, N} \sim 10^{-73}$ is strictly negligible.
>
> 5. **Strategic Conclusion for Phase IX:**  
>    The $K^2$-commutator / boundary-defect route ($\alpha_N = o(N^{-1/2})$) has reached the point of diminishing returns. The observed empirical decay $\alpha_N \sim N^{-0.29}$ at fixed $T = 600$ is sub-critical. Phase IX is officially closed, and the project pivots to the direct Gate 1 target: $\Delta_j(N) R_{\mathrm{spec}}(N, L) \to 0$.

---

## 1. Mathematical Structure & Exact Definitions

For $k \in \{1, \dots, N-1\}$, the off-diagonal divided-difference kernel decomposes identically as:
$$F_{N, k} \equiv H_{N, k} k^2 v_{N, k} = \frac{2(N\psi(N) - k\psi(k))}{N^2 - k^2} k^2 v_{N, k} \equiv A_{N, k} - B_{N, k},$$
where:
$$A_{N, k} \equiv a_N \frac{k^2 v_{N, k}}{N^2 - k^2} = \frac{2 N \psi(N) k^2 v_{N, k}}{N^2 - k^2}, \qquad B_{N, k} \equiv a_k \frac{k^2 v_{N, k}}{N^2 - k^2} = \frac{2 k^3 \psi(k) v_{N, k}}{N^2 - k^2},$$
with $a_m \equiv 2 m \psi(m) = \sqrt{2} m^2 H_{0, m}$.

At $k = N$, the diagonal entry is:
$$F_{N, N} \equiv H_{N, N} N^2 v_{N, N} = \left(\psi'(N) + \frac{\psi(N)}{N}\right) N^2 v_{N, N}.$$

The total boundary flux is:
$$(H u_N)_N = \Sigma_A^{\mathrm{off}}(N) - \Sigma_B^{\mathrm{off}}(N) + F_{N, N},$$
where:
$$\Sigma_A^{\mathrm{off}}(N) \equiv \sum_{k=1}^{N-1} A_{N, k}, \qquad \Sigma_B^{\mathrm{off}}(N) \equiv \sum_{k=1}^{N-1} B_{N, k}.$$

---

## 2. Audited Numerical Results (`cell119.out`)

### 2.1 Module 1: Low-Mode Anatomy & Peak Mechanism ($k \le 10$)

At low modes $k \ll N$:
$$\frac{B_{N, k}}{A_{N, k}} = \frac{a_k}{a_N} = \frac{k \psi(k)}{N \psi(N)} \sim \mathcal{O}\left(\frac{k}{N}\right) \ll 1.$$

| Dimension $N$ | $A_{N, 3}$ | $B_{N, 3}$ | $F_{N, 3} = A - B$ | $B_{N, 3} / A_{N, 3}$ | $A_{N, 3}$ Share $\%$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $+4.4909 \times 10^{-3}$ | $+5.1306 \times 10^{-4}$ | $+3.9778 \times 10^{-3}$ | $+0.1142$ | $112.9\%$ |
| **96** | $-6.4804 \times 10^{-3}$ | $+2.2348 \times 10^{-4}$ | $-6.7039 \times 10^{-3}$ | $-0.0345$ | $96.7\%$ |
| **128** | $-9.2987 \times 10^{-3}$ | $+1.2526 \times 10^{-4}$ | $-9.4239 \times 10^{-3}$ | $-0.0135$ | $98.7\%$ |
| **160** | $+3.2127 \times 10^{-3}$ | $+7.9985 \times 10^{-5}$ | $+3.1327 \times 10^{-3}$ | $+0.0249$ | $102.6\%$ |
| **192** | $-1.3409 \times 10^{-2}$ | $+5.5512 \times 10^{-5}$ | $-1.3465 \times 10^{-2}$ | **$-0.00414$** | **$99.59\%$** |

- **Validation:** At $N = 192$, $A_{N, 3}$ accounts for $99.59\%$ of $F_{N, 3}$, with $B$ suppressed to $< 0.5\%$. This rigorously confirms that the origin of the initial $k = 3$ excursion is driven almost entirely by the boundary symbol $a_N$ acting on $k^2 v_k$.

---

### 2.2 Module 2: Total Off-Diagonal Sums Separation Across All 8 Dimensions

| $N$ | $\Sigma_A^{\mathrm{off}}(N)$ | $\Sigma_B^{\mathrm{off}}(N)$ | $\Delta \Sigma = \Sigma_A - \Sigma_B$ | $F_{N, N}$ (diag) | Net Flux $(H u_N)_N$ | $\mathcal{C}_{AB}$ Factor |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $-8.3033 \times 10^{-24}$ | $+1.9765 \times 10^{-22}$ | $-2.0596 \times 10^{-22}$ | $-4.0826 \times 10^{-23}$ | $-2.4678 \times 10^{-22}$ | $0.9597$ |
| **80** | $+4.7282 \times 10^{-24}$ | $+2.0517 \times 10^{-22}$ | $-2.0044 \times 10^{-22}$ | $+6.6140 \times 10^{-24}$ | $-1.9383 \times 10^{-22}$ | $1.0236$ |
| **96** | $+7.1195 \times 10^{-24}$ | $+1.7131 \times 10^{-22}$ | $-1.6419 \times 10^{-22}$ | $-6.2197 \times 10^{-24}$ | $-1.7041 \times 10^{-22}$ | $1.0434$ |
| **112** | $+7.2292 \times 10^{-24}$ | $+2.1866 \times 10^{-22}$ | $-2.1143 \times 10^{-22}$ | $+4.0159 \times 10^{-23}$ | $-1.7127 \times 10^{-22}$ | $1.0342$ |
| **128** | $+3.8451 \times 10^{-23}$ | $+1.7295 \times 10^{-22}$ | $-1.3449 \times 10^{-22}$ | $-1.1275 \times 10^{-23}$ | $-1.4577 \times 10^{-22}$ | $1.2859$ |
| **144** | $+3.6094 \times 10^{-23}$ | $+1.4005 \times 10^{-22}$ | $-1.0396 \times 10^{-22}$ | $-4.4512 \times 10^{-23}$ | $-1.4847 \times 10^{-22}$ | $1.3472$ |
| **160** | $-1.5218 \times 10^{-23}$ | $+1.3433 \times 10^{-22}$ | $-1.4955 \times 10^{-22}$ | $-3.5977 \times 10^{-23}$ | $-1.8553 \times 10^{-22}$ | $0.8982$ |
| **192** | $+4.2807 \times 10^{-23}$ | $+1.5679 \times 10^{-22}$ | $-1.1398 \times 10^{-22}$ | $+3.0492 \times 10^{-23}$ | $-8.3491 \times 10^{-23}$ | **$1.3756$** |

- **Decisive Falsification of H1:** $\Sigma_A^{\mathrm{off}}$ and $\Sigma_B^{\mathrm{off}}$ are **not** large cancelling numbers ($\sim 10^{-2}$). Both sums are individually already $\sim 10^{-23} - 10^{-22}$! The cancellation ratio $\mathcal{C}_{AB} = \max(|\Sigma_A|, |\Sigma_B|) / |\Delta \Sigma|$ remains bounded between $0.89$ and $1.38$.

---

### 2.3 Module 3: Cumulative Trajectory Tracking at $N = 192$

| $x = k/N$ | Mode $m$ | $S_A(x) = \sum_{k \le m} A_k$ | $S_B(x) = \sum_{k \le m} B_k$ | $S_F(x) = S_A - S_B$ | Ratio $|S_F| / \max(|S_A|, |S_B|)$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **0.01** | $1$ | $+3.6082 \times 10^{-3}$ | $-1.5241 \times 10^{-6}$ | $+3.6097 \times 10^{-3}$ | $1.0004$ |
| **0.02** | $3$ | $-1.1369 \times 10^{-2}$ | $+5.6714 \times 10^{-5}$ | $-1.1426 \times 10^{-2}$ | $1.0050$ |
| **0.03** | $5$ | $-8.0626 \times 10^{-4}$ | $+6.0692 \times 10^{-5}$ | $-8.6695 \times 10^{-4}$ | $1.0753$ |
| **0.05** | $9$ | $-9.5589 \times 10^{-5}$ | $+3.7942 \times 10^{-6}$ | $-9.9384 \times 10^{-5}$ | $1.0397$ |
| **0.10** | $19$ | $+4.9263 \times 10^{-10}$ | $-6.0977 \times 10^{-12}$ | $+4.9873 \times 10^{-10}$ | $1.0124$ |
| **0.15** | $28$ | $-2.1259 \times 10^{-14}$ | $+8.1036 \times 10^{-16}$ | $-2.2070 \times 10^{-14}$ | $1.0381$ |
| **0.20** | $38$ | $-2.8233 \times 10^{-19}$ | $+5.5818 \times 10^{-20}$ | $-3.3815 \times 10^{-19}$ | $1.1977$ |
| **0.30** | $57$ | $+4.9643 \times 10^{-23}$ | $+1.7809 \times 10^{-22}$ | $-1.2844 \times 10^{-22}$ | $0.7212$ |
| **0.50** | $96$ | $+4.8850 \times 10^{-23}$ | $+1.7839 \times 10^{-22}$ | $-1.2954 \times 10^{-22}$ | $0.7262$ |
| **0.80** | $153$ | $+5.2146 \times 10^{-23}$ | $+1.7867 \times 10^{-22}$ | $-1.2653 \times 10^{-22}$ | $0.7081$ |
| **0.95** | $182$ | $+5.7112 \times 10^{-23}$ | $+1.7942 \times 10^{-22}$ | $-1.2231 \times 10^{-22}$ | $0.6817$ |

- **Anatomy of the Cancellation:**  
  Notice that $S_A(x)$ drops from $-1.14 \times 10^{-2}$ at $x = 0.02$ down to $-2.82 \times 10^{-19}$ at $x = 0.20$ and $+4.96 \times 10^{-23}$ at $x = 0.30$. **The entire 20 decades of cancellation occur internally within $A$ alone.** $B$ never exceeds $6 \times 10^{-5}$ and settles to $+1.78 \times 10^{-22}$ by $x = 0.30$.

---

### 2.4 Module 4: Moments Balance and Failure of Naive Taylor Expansion

| $N$ | $M_2(N) = \sum k^2 v_k$ | $M_{\psi, 3}(N) = 2 \sum k^3 \psi v_k$ | $\Sigma_A^{\mathrm{exact}}$ | $\Sigma_A^{\mathrm{lead}}$ | Rel Err $A$ | $\Sigma_B^{\mathrm{exact}}$ | $\Sigma_B^{\mathrm{lead}}$ | Rel Err $B$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $-2.8402 \times 10^{-19}$ | $+6.8625 \times 10^{-17}$ | $-8.3033 \times 10^{-24}$ | $-5.1659 \times 10^{-22}$ | $61.2\times$ | $+1.9765 \times 10^{-22}$ | $+1.6754 \times 10^{-20}$ | $83.8\times$ |
| **96** | $-2.2259 \times 10^{-19}$ | $+5.4319 \times 10^{-17}$ | $+7.1195 \times 10^{-24}$ | $+5.9612 \times 10^{-22}$ | $82.7\times$ | $+1.7131 \times 10^{-22}$ | $+5.8940 \times 10^{-21}$ | $33.4\times$ |
| **128** | $-2.1632 \times 10^{-19}$ | $+5.3025 \times 10^{-17}$ | $+3.8451 \times 10^{-23}$ | $+8.3426 \times 10^{-22}$ | $20.7\times$ | $+1.7295 \times 10^{-22}$ | $+3.2364 \times 10^{-21}$ | $17.7\times$ |
| **160** | $-2.1932 \times 10^{-19}$ | $+5.4003 \times 10^{-17}$ | $-1.5218 \times 10^{-23}$ | $-2.9288 \times 10^{-22}$ | $18.3\times$ | $+1.3433 \times 10^{-22}$ | $+2.1095 \times 10^{-21}$ | $14.7\times$ |
| **192** | $-2.1727 \times 10^{-19}$ | $+5.3603 \times 10^{-17}$ | $+4.2807 \times 10^{-23}$ | $+1.2117 \times 10^{-21}$ | $27.3\times$ | $+1.5679 \times 10^{-22}$ | $+1.4541 \times 10^{-21}$ | $8.27\times$ |

- **Warning on Asymptotics:** Naive termwise Taylor expansion fails to approximate the sums because the sum is dominated by delicate phase cancellation across higher modes, rather than the low-mode power series.

---

### 2.5 Module 5: Exact Theorem 1 Synthesis

| $N$ | $\alpha_N$ | $\frac{a_N}{\sqrt{2}} T_v(0)$ | $E_1 N^2 v_{N, N}$ | Net Flux $(H u_N)_N$ | Theorem 1 Identity Residual |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $-2.4219 \times 10^{-22}$ | $+4.5943 \times 10^{-24}$ | $-4.60 \times 10^{-73}$ | $-2.4678 \times 10^{-22}$ | $9.55 \times 10^{-71}$ |
| **80** | $-1.9705 \times 10^{-22}$ | $-3.2247 \times 10^{-24}$ | $+7.78 \times 10^{-74}$ | $-1.9383 \times 10^{-22}$ | $1.51 \times 10^{-70}$ |
| **96** | $-1.8163 \times 10^{-22}$ | $-1.1216 \times 10^{-23}$ | $-1.52 \times 10^{-73}$ | $-1.7041 \times 10^{-22}$ | $7.75 \times 10^{-70}$ |
| **112** | $-1.7951 \times 10^{-22}$ | $-8.2443 \times 10^{-24}$ | $+6.08 \times 10^{-73}$ | $-1.7127 \times 10^{-22}$ | $4.16 \times 10^{-70}$ |
| **128** | $-1.7287 \times 10^{-22}$ | $-2.7097 \times 10^{-23}$ | $-1.29 \times 10^{-73}$ | $-1.4577 \times 10^{-22}$ | $8.77 \times 10^{-70}$ |
| **144** | $-1.7476 \times 10^{-22}$ | $-2.6297 \times 10^{-23}$ | $-5.15 \times 10^{-73}$ | $-1.4847 \times 10^{-22}$ | $8.91 \times 10^{-70}$ |
| **160** | $-1.7116 \times 10^{-22}$ | $+1.4369 \times 10^{-23}$ | $-2.90 \times 10^{-73}$ | $-1.8553 \times 10^{-22}$ | $9.04 \times 10^{-71}$ |
| **192** | $-1.6798 \times 10^{-22}$ | $-8.4490 \times 10^{-23}$ | $+3.07 \times 10^{-73}$ | $-8.3491 \times 10^{-23}$ | **$2.33 \times 10^{-69}$** |

---

## 3. Epistemic Synthesis & Formal Closure of Phase IX

1. **The Observed Finite-$(T, N)$ Cancellation Structure:**  
   We have **not** discovered the continuum analytical mechanism behind the $10^{20}$-scale cancellation. What we have established is a robust finite-$T$, finite-$N$ numerical cancellation structure:
   - Airy scaling collapse fails;
   - Normalized master-curve collapse fails;
   - Simple arithmetic prime-power phase locking fails;
   - $A/B$ inter-component cancellation fails;
   - The large low-mode excursion is extinguished internally within the $A$-sum before $x \sim 0.2 - 0.3$, leaving a $10^{-22}$ boundary residual.
2. **Termination of the $K^2$-Commutator / Boundary Defect Programme:**  
   The boundary defect rate $\alpha_N = o(N^{-1/2})$ was an auxiliary sufficient condition introduced to prove $\sup_N \mathcal{K}_2(N) < \infty$. The empirical data shows $\alpha_N \sim N^{-0.29}$, which does not clear the $N^{-1/2}$ hurdle at fixed $T = 600$. Continuing with further boundary-defect micro-analysis represents diminishing returns.
3. **Formal Closure:** Phase IX (Cells 111–119) is hereby completed and closed.

---

## 4. The Strategic Road Ahead: Attacking Gate 1 Directly

The project now pivots to the core target proposition of Gate 1:
$$\lim_{L \to \infty} \limsup_{N \to \infty} \Delta_j(N) R_{\mathrm{spec}}(N, L) = 0 \implies \Pi_{j, \mathrm{tail}} \to 1.$$

The strategic fork to be audited in Cell 120 is:
- **Route 1B Proof Audit:** Does fixed-$T$ Loewner smoothness $\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T)$ rigorously guarantee $R_{\mathrm{spec}} = \mathcal{O}(1)$, reducing the entire problem to showing $\Delta_j(N) \to 0$? What is the status of the lower Ritz gap separation $\inf_N (E_{L+1}^{(N)} - E_{j+1}^{(N)}) \ge g_* > 0$?
- **Route 1A Discrete Barrier Mechanics:** If Route 1B leaves an open gap, can we derive a discrete Agmon/WKB barrier estimate $\Delta_j(N) \le C_j e^{-2 S_j N}$ ($S_j > 0$) directly from the Galerkin operator to crush any soft spectral growth?
