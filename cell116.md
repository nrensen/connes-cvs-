# Cell 116 Analytical Note: Boundary-Flux Kernel Asymptotics & The Airy Boundary Layer Reconnaissance

**Companion Computational Script:** [`cell116.py`](file:///c:/data/github/connes-cvs-/cell116.py) | **Verification Log:** [`cell116.out`](file:///c:/data/github/connes-cvs-/cell116.out) (runtime: 543.20 s at 70 dps)  
**Status:** Executed & Audited / Calibrated Analytical Reconnaissance (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell115.md`](file:///c:/data/github/connes-cvs-/cell115.md); [`cell115.out`](file:///c:/data/github/connes-cvs-/cell115.out); [`cell112.md`](file:///c:/data/github/connes-cvs-/cell112.md); [`cell111.md`](file:///c:/data/github/connes-cvs-/cell111.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Calibrated Verdict

Cell 116 executed an analytical and numerical reconnaissance of the exact boundary-flux identity:
$$(H u_N)_N = \alpha_N - \frac{T_v(0)}{\sqrt{2}} a_N + E_{11} N^2 v_{N, N}.$$
The purpose was to investigate whether the sub-critical power-law decay observed in Cell 115 ($\hat{\beta}_\alpha = -0.2885 \pm 0.0631, \hat{\beta}_T = -0.3244 \pm 0.0642$) could be traced to a specific structural mechanism in the boundary row of the Connes–CvS Galerkin operator.

### Headline Verdict

> **Classification: Promising analytical reconnaissance, not yet a derived analytical theorem.**
>
> 1. **Confirmed Discovery — Empirical Boundary Proportionality:**
>    On the localized solitary branch across $N \in [64, 192]$, the boundary coupling $\alpha_N$ and contact defect $T_v(0)$ maintain an empirical proportionality:
>    $$\kappa_\alpha(N) \equiv \frac{|\alpha_N|}{|T_v(0)|} = 284.18 \pm 3.51 \quad (\text{only } 3.99\% \text{ relative variation over a 3-fold span in } N).$$
>    The exact boundary identity was verified to a residual of $2.33 \times 10^{-69}$ at $N = 192$.
>
> 2. **Confirmed Empirical Candidate — $N^{-1/3}$ Scaling:**
>    The invariant products $N^{1/3}|\alpha_N|$ and $N^{1/3}|T_v(0)|$ exhibit modest spreads of $15.3\%$ and $17.2\%$, showing that $N^{-1/3}$ is a **highly plausible empirical model** consistent with Cell 115's regressions.
>
> 3. **Correction of Analytical Flaws & Retractions:**
>    - *The $a_N$ Explanation Is Refuted:* The hypothesis that $\kappa_\alpha \approx 2 (a_N/\sqrt{2})$ universally across all $N$ is **empirically false**. The sequence $a_N / \sqrt{2} = \sqrt{2} N \psi(N)$ oscillates wildly from $+24.17$ to $-145.38$ across $N \in [64, 192]$ due to high-frequency trigonometric terms in $\psi(N)$. The equipartition $(H u_N)_N \approx \frac{1}{2}\alpha_N \approx \frac{a_N}{\sqrt{2}} T_v(0)$ observed at $N = 192$ is an empirical feature at $N=192$, not an exact identity for all $N$.
>    - *Airy Model Is Heuristic, Not Derived:* The Airy boundary layer model ($\delta t \sim N^{-2/3}$) was imported from model Schrödinger equations with a hard wall, not derived from the Connes–CvS Galerkin kernel. The oscillatory integral factor does not tend to a constant.
>    - *Massive Cancellation in Boundary Flux:* Decomposing $(H u_N)_N = \sum_{k=1}^N H_{Nk} k^2 v_k$ reveals that the low-mode core ($k \le 24$) contributes $-5.03 \times 10^{-13}$, which is **ten orders of magnitude larger** than the net flux ($-8.35 \times 10^{-23}$). The boundary flux is a delicate cancellation problem, not a simple boundary-layer dominated sum.
>    - *Retraction of Premature Conclusions:* We retract claims that fixed-$T$ boundary extinction is mathematically impossible and that the joint limit $T(N) \ge 2\pi N/L$ eliminates the defect. These remain open working hypotheses.

---

## 1. Audited Numerical Results (`cell116.out`)

### 1.1 Module 1: The Boundary Proportionality Ratio

| $N$ | $|\alpha_N|$ | $|T_v(0)|$ | $\kappa_\alpha = |\alpha_N| / |T_v(0)|$ | $a_N / \sqrt{2}$ | $\kappa / (a_N / \sqrt{2})$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $2.4219 \times 10^{-22}$ | $8.7212 \times 10^{-25}$ | $277.70$ | $+5.27$ | $+52.71$ |
| **80** | $1.9705 \times 10^{-22}$ | $7.0171 \times 10^{-25}$ | $280.82$ | $-4.60$ | $-61.11$ |
| **96** | $1.8163 \times 10^{-22}$ | $6.4269 \times 10^{-25}$ | $282.61$ | $-17.45$ | $-16.19$ |
| **112** | $1.7951 \times 10^{-22}$ | $6.3201 \times 10^{-25}$ | $284.03$ | $-13.04$ | $-21.77$ |
| **128** | $1.7287 \times 10^{-22}$ | $6.0647 \times 10^{-25}$ | $285.04$ | $-44.68$ | $-6.38$ |
| **144** | $1.7476 \times 10^{-22}$ | $6.1060 \times 10^{-25}$ | $286.22$ | $-43.07$ | $-6.65$ |
| **160** | $1.7116 \times 10^{-22}$ | $5.9439 \times 10^{-25}$ | $287.96$ | $+24.17$ | $+11.91$ |
| **192** | $1.6798 \times 10^{-22}$ | $5.8117 \times 10^{-25}$ | $289.04$ | $-145.38$ | $-1.99$ |

- **Mean $\kappa_\alpha$:** $284.18 \pm 3.51$ (relative variation: $3.99\%$).
- **Crucial Diagnostic Correction:**
  Notice that $a_N / \sqrt{2}$ changes sign and magnitude across dimensions ($-145.38$ to $+24.17$). Consequently, $\kappa / (a_N / \sqrt{2})$ ranges from $-61$ to $+53$. The earlier hypothesis that $\kappa_\alpha \approx 2 (a_N/\sqrt{2})$ was an accidental feature of $N = 192$ and is **not supported** across $N \in [64, 160]$.

### 1.2 Module 2: Invariant Product Evaluation ($N^{-1/3}$ Model)

| $N$ | $N^{1/3} |\alpha_N|$ | $N^{1/3} |T_v(0)|$ | $N^{-1/6} P_\alpha(N)$ | $N^{-7/6} P_T(N)$ |
| :---: | :---: | :---: | :---: | :---: |
| **64** | $9.688 \times 10^{-22}$ | $3.488 \times 10^{-24}$ | $9.688 \times 10^{-22}$ | $3.488 \times 10^{-24}$ |
| **80** | $8.491 \times 10^{-22}$ | $3.024 \times 10^{-24}$ | $8.491 \times 10^{-22}$ | $3.024 \times 10^{-24}$ |
| **96** | $8.317 \times 10^{-22}$ | $2.943 \times 10^{-24}$ | $8.317 \times 10^{-22}$ | $2.943 \times 10^{-24}$ |
| **112** | $8.653 \times 10^{-22}$ | $3.046 \times 10^{-24}$ | $8.653 \times 10^{-22}$ | $3.046 \times 10^{-24}$ |
| **128** | $8.712 \times 10^{-22}$ | $3.056 \times 10^{-24}$ | $8.712 \times 10^{-22}$ | $3.056 \times 10^{-24}$ |
| **144** | $9.160 \times 10^{-22}$ | $3.200 \times 10^{-24}$ | $9.160 \times 10^{-22}$ | $3.200 \times 10^{-24}$ |
| **160** | $9.292 \times 10^{-22}$ | $3.227 \times 10^{-24}$ | $9.292 \times 10^{-22}$ | $3.227 \times 10^{-24}$ |
| **192** | $9.691 \times 10^{-22}$ | $3.353 \times 10^{-24}$ | $9.691 \times 10^{-22}$ | $3.353 \times 10^{-24}$ |

- $N^{1/3}|\alpha_N|$ mean: $9.00 \times 10^{-22}$ (spread: $15.27\%$).
- $N^{1/3}|T_v(0)|$ mean: $3.17 \times 10^{-24}$ (spread: $17.23\%$).
This indicates that $N^{-1/3}$ is an excellent empirical envelope over $64 \le N \le 192$, but does not yet establish a derivation.

### 1.3 Module 3: Boundary-Row Modal Flux Cancellation at $N = 192$

| Sector / Partition | Index Range | Partial Sum $\sum_{k \le M} H_{Nk} k^2 v_k$ | $\%$ of Total Flux |
| :--- | :---: | :---: | :---: |
| **Core Sector $k \le 24$** | $k \in [1, 24]$ | **$-5.0305 \times 10^{-13}$** | **$+6.025 \times 10^{11} \%$** |
| **Sector $k \le 48$** | $k \in [1, 48]$ | $-5.1524 \times 10^{-23}$ | $+61.71 \%$ |
| **Sector $k \le 96$** | $k \in [1, 96]$ | $-1.2954 \times 10^{-22}$ | $+155.15 \%$ |
| **Sector $k \le 144$** | $k \in [1, 144]$ | $-1.2960 \times 10^{-22}$ | $+155.22 \%$ |
| **Sector $k \le 180$** | $k \in [1, 180]$ | $-1.2278 \times 10^{-22}$ | $+147.05 \%$ |
| **Sector $k \le 191$** | $k \in [1, 191]$ | $-1.1398 \times 10^{-22}$ | $+136.52 \%$ |
| **Endpoint Diagonal $k = 192$** | $k = 192$ | $+3.0492 \times 10^{-23}$ | $-36.52 \%$ |
| **Total Flux $(H u_N)_N$** | $k \in [1, 192]$ | **$-8.3491 \times 10^{-23}$** | **$100.00 \%$** |

**Crucial Finding on Massive Cancellation:**  
The partial sum of the low-mode core ($k \le 24$) is $-5.03 \times 10^{-13}$, while the total boundary flux is $-8.35 \times 10^{-23}$. There is a **ten-order-of-magnitude destructive cancellation** across intermediate modes ($25 \le k \le 48$) that quenches the core flux down to the $10^{-23}$ scale. Bounding the boundary flux as a simple "boundary layer dominance" without treating this cancellation is incomplete.

### 1.4 Module 4: Equipartition at $N = 192$
At $N = 192$:
- Computed boundary flux $(H u_N)_N = -8.349068 \times 10^{-23}$
- Boundary coupling $\alpha_N = -1.679805 \times 10^{-22}$
- Contact term $\frac{a_N}{\sqrt{2}} T_v(0) = -8.448987 \times 10^{-23}$
- Eigenvalue term $E_{11} N^2 v_{N, N} = 3.07 \times 10^{-73}$
- Exact Identity Residual: $2.33 \times 10^{-69}$.
$$\frac{(H u_N)_N}{\alpha_N} = 0.4970, \qquad \frac{(H u_N)_N}{\frac{a_N}{\sqrt{2}} T_v(0)} = 0.9882.$$
The equipartition $(H u_N)_N \approx \frac{1}{2} \alpha_N \approx \frac{a_N}{\sqrt{2}} T_v(0)$ is exact to within $1.2\%$ at $N = 192$.

---

## 2. Theoretical Analysis & Identified Gaps

### 2.1 Why $a_N$ Oscillates
Recall that $a_N = 2 N \psi(N)$. The symbol $\psi(N)$ is given by (cf. `operator.py`):
$$\psi(N) = \psi_{\mathrm{prime}}(N) + \psi_{\mathrm{pole}}(N) + \psi_{\mathrm{arch}}(N).$$
The prime and pole terms contain discrete Fourier sums over prime powers $p^k \le c$:
$$\psi_{\mathrm{prime}}(N) \sim \sum_{p^k \le c} \frac{\log p}{p^{k/2}} \sin\left( \frac{2\pi N \log(p^k)}{L} \right).$$
These are rapidly oscillating trigonometric sums in $N$. Furthermore, the Archimedean piece $\psi_{\mathrm{arch}}(N)$ involves an integral with a singularity at $\alpha_N = \frac{2\pi N}{L}$.
Consequently, $\psi(N)$ is an oscillating sequence with changing sign and magnitude. Multiplying by $2N$ produces the wild swings in $a_N / \sqrt{2}$ seen in Module 1.

### 2.2 Mathematical Gaps in the Heuristic Airy Model
1. **Model Schrödinger Equation vs Galerkin Kernel:** Writing $-\frac{1}{N^2} \phi'' + F t \phi = E \phi$ imports a local differential equation that has not been derived from the non-local Connes–CvS integral operator.
2. **Identification of $\hbar$:** Treating $1/N$ as an effective semiclassical Planck constant assumes a continuous WKB limit that has not been justified for discrete finite-$N$ Galerkin truncations.
3. **Oscillatory Truncation Factor:** The integral $\int_0^{\delta t} \frac{\sin(2\pi N t)}{t} \phi(t) dt$ produces a factor $\frac{1 - \cos(2\pi N^{1/3})}{2\pi}$, which does not converge to a constant as $N \to \infty$.
4. **Boundary Condition Assumption:** Dirichlet boundary condition $T_\infty(0) = 0$ does not imply $T_\infty'(0) = 0$. Assuming the vanishing of the derivative without analytical derivation from the continuum integral equation is unfounded.
5. **Bulk Expansion Range Warning:** In Module 3, partition $k \le 96$ at $N = 192$ reaches $k = N/2$, which is not $k \ll N$. The expansion $\frac{k^2}{N^2 - k^2} = \frac{k^2}{N^2} + \mathcal{O}(N^{-4})$ holds only for $k/N \to 0$, so intermediate modes cannot be treated as asymptotic bulk modes.

---

## 3. The Three Possibilities for Gate 1

Rather than declaring Gate 1 closed or impossible at fixed $T$, Cell 116 isolates three clear, testable branches:

- **Possibility A (Genuine Fixed-$T$ Airy Barrier):** If the discrete Galerkin operator truly possesses an $N^{-1/3}$ boundary-layer scaling law, then $P_\alpha(N) \sim N^{1/6}$ diverges, proving that fixed-$T$ boundary extinction cannot occur. This would mandate the joint limit $(N, T) \to \infty$ with $T \ge 2\pi N / L$.
- **Possibility B (Pre-Asymptotic Transient):** The observed $N^{-0.29}$ regime over $N \in [64, 192]$ may be a pre-asymptotic transient that crosses over to super-critical decay ($N^{-\beta}$ with $\beta > 0.5$) at larger dimensions.
- **Possibility C (Cancellation-Driven Extinction):** The ten-order-of-magnitude cancellation seen in Module 3 between core and tail modes may produce a net boundary flux that decays faster than the individual sector envelopes asymptotically.

---

## 4. Forward Mandate for Cell 117

To determine whether the Airy boundary layer is an actual physical feature of the operator or merely an empirical coincidence:

1. **Modal Boundary-Row Profiling:**
   Compute individual modal flux terms:
   $$F_{N, k} \equiv H_{Nk} k^2 v_{N, k}$$
   and re-index by distance from the boundary:
   $$j \equiv N - k \qquad (j \in \{0, 1, 2, \dots, 24\}).$$
2. **Multi-Dimension Scaling Collapse Test:**
   Evaluate $F_{N, N-j}$ across $N \in \{64, 96, 128, 160, 192\}$.
   Test whether the boundary terms collapse onto a universal profile:
   $$F_{N, N-j} \stackrel{?}{\approx} N^{-\mu} f(j / N^{1/3}).$$
3. **Analytical Dissection of $H_{0, N}$:**
   Analytically derive the formula for $H_{0, N}$ and $\psi(N)$ to explain the exact oscillatory structure of $a_N$.
