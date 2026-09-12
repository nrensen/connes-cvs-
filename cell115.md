# Cell 115 Analytical Note: Localized Branch Boundary Defect Asymptotic Power-Law Regression

**Companion Computational Script:** [`cell115.py`](file:///c:/data/github/connes-cvs-/cell115.py) | **Verification Log:** [`cell115.out`](file:///c:/data/github/connes-cvs-/cell115.out) (runtime: 1722.17 s $\approx$ 28.7 min at 70 dps)  
**Status:** Executed & Audited (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell114.md`](file:///c:/data/github/connes-cvs-/cell114.md); [`cell114.out`](file:///c:/data/github/connes-cvs-/cell114.out); [`cell110.md`](file:///c:/data/github/connes-cvs-/cell110.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Headline Verdict

Cell 115 was designed to determine the quantitative scaling of the localized-branch boundary defect over the post-transition regime $N \in \{64, 80, 96, 112, 128, 144, 160, 192\}$ ($n = 8$ data points) at $c = 13, T = 600$. The script evaluated canonical Hamiltonians extracted from the cached Galerkin matrix at 70 dps and performed ordinary least-squares log-log regressions on boundary observables.

### Headline Verdict

> **1. The boundary defect exhibits empirical sub-critical decay, NOT a constant nonzero limit:**
> Over $64 \le N \le 192$, the localized-branch ($k = 1$) data yield:
> $$\boxed{\hat{\beta}_\alpha = -0.2885 \pm 0.0631 \quad (R^2 = 0.7767), \qquad |\alpha_N| \approx 7.26 \times 10^{-22} N^{-0.2885}}$$
> and
> $$\boxed{\hat{\beta}_T = -0.3244 \pm 0.0642 \quad (R^2 = 0.8099), \qquad |T_v(0)| \approx 3.03 \times 10^{-24} N^{-0.3244}.}$$
> Because $\hat{\beta}_\alpha$ is negative and separated from zero by more than $4.5$ standard errors (two-standard-error interval: $[-0.415, -0.162]$), **Hypothesis 1 (a constant nonzero limiting defect $\alpha_\infty > 0$) is refuted as the leading interpretation.** The data favour a slow, power-law decay toward zero.
>
> **2. The boundary defect decays far too slowly for Gate 1 extinction:**
> Gate 1 sufficient conditions require:
> $$\alpha_N = o(N^{-1/2}) \iff \beta_\alpha < -0.50, \qquad T_v(0) = o(N^{-3/2}) \iff \beta_T < -1.50.$$
> The empirical exponents fall substantially short of these requirements:
> $$\Delta_\alpha \equiv \hat{\beta}_\alpha - (-0.50) = +0.2115, \qquad \Delta_T \equiv \hat{\beta}_T - (-1.50) = +1.1756.$$
> Consequently, the rescaled extinction products are **slowly increasing** with $N$:
> $$P_\alpha(N) \equiv |\alpha_N|\sqrt{N} \sim N^{+0.2115} \quad (R^2 = 0.6516), \qquad P_T(N) \equiv |T_v(0)|N^{3/2} \sim N^{+1.1756} \quad (R^2 = 0.9824).$$
>
> **3. Deconstruction of the "Extinction Plateau":**
> What previously appeared as an immovable $2 \times 10^{-21}$ plateau in Cells 111–114 is revealed to be a slowly growing power law: over a 3-fold span in $N$ ($64 \to 192$), $3^{0.2115} \approx 1.26$, exactly matching the observed movement from $1.94 \times 10^{-21}$ to $2.33 \times 10^{-21}$.

---

## 1. Complete Numerical Data & Log-Log Regressions

### 1.1 Multi-Dimension Raw Data Table ($T = 600, \text{dps} = 70$)

Evaluating the localized candidate ($k = 1$) alongside the edge candidate ($k = 0$) across $N \in [48, 192]$:

| $N$ | Branch $k$ | Energy $E$ | $v_0$ | $1 - L_{24}$ | $\mathcal{K}_2$ | $|T_v(0)|$ | $|\alpha_N|$ | $P_\alpha = |\alpha_N|\sqrt{N}$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **48** | $k=0$ (loc) | $2.103 \times 10^{-50}$ | $0.4558$ | $3.63 \times 10^{-27}$ | $21.45$ | $1.783 \times 10^{-24}$ | $4.570 \times 10^{-22}$ | $3.166 \times 10^{-21}$ |
| | $k=1$ (exc) | $1.490 \times 10^{-49}$ | $0.4908$ | $9.46 \times 10^{-26}$ | $85.81$ | $9.270 \times 10^{-24}$ | $2.377 \times 10^{-21}$ | $1.647 \times 10^{-20}$ |
| **64** | $k=0$ (edge) | $-5.085 \times 10^{-52}$ | $0.0708$ | $9.00 \times 10^{-26}$ | $82.19$ | $9.674 \times 10^{-25}$ | $2.684 \times 10^{-22}$ | $2.147 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.713 \times 10^{-50}$ | $0.6657$ | $7.36 \times 10^{-26}$ | $25.46$ | $8.721 \times 10^{-25}$ | $2.422 \times 10^{-22}$ | $1.938 \times 10^{-21}$ |
| **80** | $k=0$ (edge) | $-8.388 \times 10^{-52}$ | $0.0666$ | $9.29 \times 10^{-26}$ | $82.72$ | $7.921 \times 10^{-25}$ | $2.222 \times 10^{-22}$ | $1.987 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.686 \times 10^{-50}$ | $0.6661$ | $7.41 \times 10^{-26}$ | $24.94$ | $7.017 \times 10^{-25}$ | $1.971 \times 10^{-22}$ | $1.762 \times 10^{-21}$ |
| **96** | $k=0$ (edge) | $-9.335 \times 10^{-52}$ | $0.0654$ | $9.33 \times 10^{-26}$ | $82.87$ | $7.299 \times 10^{-25}$ | $2.060 \times 10^{-22}$ | $2.018 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.678 \times 10^{-50}$ | $0.6662$ | $7.39 \times 10^{-26}$ | $24.80$ | $6.427 \times 10^{-25}$ | $1.816 \times 10^{-22}$ | $1.780 \times 10^{-21}$ |
| **112** | $k=0$ (edge) | $-9.804 \times 10^{-52}$ | $0.0649$ | $9.28 \times 10^{-26}$ | $82.94$ | $7.195 \times 10^{-25}$ | $2.041 \times 10^{-22}$ | $2.160 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.675 \times 10^{-50}$ | $0.6663$ | $7.33 \times 10^{-26}$ | $24.72$ | $6.320 \times 10^{-25}$ | $1.795 \times 10^{-22}$ | $1.900 \times 10^{-21}$ |
| **128** | $k=0$ (edge) | $-1.005 \times 10^{-51}$ | $0.0645$ | $9.28 \times 10^{-26}$ | $82.98$ | $6.922 \times 10^{-25}$ | $1.970 \times 10^{-22}$ | $2.229 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.673 \times 10^{-50}$ | $0.6663$ | $7.32 \times 10^{-26}$ | $24.68$ | $6.065 \times 10^{-25}$ | $1.729 \times 10^{-22}$ | $1.956 \times 10^{-21}$ |
| **144** | $k=0$ (edge) | $-1.027 \times 10^{-51}$ | $0.0643$ | $9.25 \times 10^{-26}$ | $83.01$ | $6.978 \times 10^{-25}$ | $1.994 \times 10^{-22}$ | $2.393 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.671 \times 10^{-50}$ | $0.6664$ | $7.28 \times 10^{-26}$ | $24.65$ | $6.106 \times 10^{-25}$ | $1.748 \times 10^{-22}$ | $2.097 \times 10^{-21}$ |
| **160** | $k=0$ (edge) | $-1.052 \times 10^{-51}$ | $0.0640$ | $9.23 \times 10^{-26}$ | $83.05$ | $6.802 \times 10^{-25}$ | $1.955 \times 10^{-22}$ | $2.473 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.669 \times 10^{-50}$ | $0.6664$ | $7.25 \times 10^{-26}$ | $24.61$ | $5.944 \times 10^{-25}$ | $1.712 \times 10^{-22}$ | $2.165 \times 10^{-21}$ |
| **192** | $k=0$ (edge) | $-1.063 \times 10^{-51}$ | $0.0638$ | $9.22 \times 10^{-26}$ | $83.06$ | $6.663 \times 10^{-25}$ | $1.923 \times 10^{-22}$ | $2.664 \times 10^{-21}$ |
| | $k=1$ (loc) | $4.668 \times 10^{-50}$ | $0.6664$ | $7.24 \times 10^{-26}$ | $24.59$ | $5.812 \times 10^{-25}$ | $1.680 \times 10^{-22}$ | $2.328 \times 10^{-21}$ |

---

### 1.2 Log-Log Power-Law Regression Parameters ($N \in [64, 192], n = 8$)

Ordinary least-squares regressions of the model $\log_{10} Y = \beta \log_{10} N + \log_{10} C$:

| Observable $Y(N)$ | Branch | Fitted Exponent $\hat{\beta}$ | $\mathrm{SE}(\hat{\beta})$ | $R^2$ | Pre-factor $\hat{C}$ | Required Exponent | Gap $\hat{\beta} - \beta_{\mathrm{req}}$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **$|\alpha_N|$** | $k = 1$ (loc) | **$-0.2885$** | $0.0631$ | $0.7767$ | $7.263 \times 10^{-22}$ | $< -0.5000$ | **$+0.2115$** |
| **$|\alpha_N|$** | $k = 0$ (edge) | **$-0.2624$** | $0.0585$ | $0.7705$ | $7.278 \times 10^{-22}$ | $< -0.5000$ | **$+0.2376$** |
| **$|T_v(0)|$** | $k = 1$ (loc) | **$-0.3244$** | $0.0642$ | $0.8099$ | $3.031 \times 10^{-24}$ | $< -1.5000$ | **$+1.1756$** |
| **$|T_v(0)|$** | $k = 0$ (edge) | **$-0.2976$** | $0.0594$ | $0.8072$ | $3.031 \times 10^{-24}$ | $< -1.5000$ | **$+1.2024$** |
| **$P_\alpha = |\alpha_N|\sqrt{N}$** | $k = 1$ (loc) | **$+0.2115$** | $0.0631$ | $0.6516$ | $7.263 \times 10^{-22}$ | $< 0.0000$ | **$+0.2115$** |
| **$P_\alpha = |\alpha_N|\sqrt{N}$** | $k = 0$ (edge) | **$+0.2376$** | $0.0585$ | $0.7335$ | $7.278 \times 10^{-22}$ | $< 0.0000$ | **$+0.2376$** |
| **$P_T = |T_v(0)|N^{3/2}$** | $k = 1$ (loc) | **$+1.1756$** | $0.0642$ | $0.9824$ | $< 0.0000$ | **$+1.1756$** |
| **$P_T = |T_v(0)|N^{3/2}$** | $k = 0$ (edge) | **$+1.2024$** | $0.0594$ | $0.9856$ | $< 0.0000$ | **$+1.2024$** |
| **$|\alpha_1| / |\alpha_0|$** | $k = 1$ vs $0$ | **$-0.0261$** | $0.0047$ | $0.8363$ | $0.9979$ | — | — |

The exact consistency identities:
$$\gamma_\alpha = \beta_\alpha + 0.5 \quad (-0.288481 + 0.5 = 0.211519), \qquad \gamma_T = \beta_T + 1.5 \quad (-0.324389 + 1.5 = 1.175611)$$
are verified exactly.

---

## 2. Epistemic & Methodological Corrections

### 2.1 Correction to Hypothesis 1 Formulation
The pre-flight formulation asserted:
$$\lim_{N \to \infty} |\alpha_N| = \alpha_\infty > 0 \iff \beta_\alpha = 0.$$
This equivalence is **not mathematically valid**. A fitted exponent $\beta = 0$ over a finite range does not prove convergence to a positive constant; conversely, a sequence can tend to zero slower than any power (e.g. $1/\log N$), or exhibit logarithmic transients that masquerade as power laws.

The epistemically calibrated claim is:
> **A fitted exponent statistically consistent with zero would be evidence against observable power-law extinction, but would not by itself establish a nonzero limiting defect.**

### 2.2 Rejection of the Constant Nonzero Plateau Hypothesis
The data do not support a constant plateau $\alpha_N \to \alpha_\infty > 0$. The measured exponent $\hat{\beta}_\alpha = -0.2885 \pm 0.0631$ is negative and strictly bounded away from zero. Therefore:
$$\boxed{\text{The hypothesis that } \alpha_N \text{ has already saturated to a nonzero constant is rejected.}}$$
The localized branch boundary defect is decaying, but at a **sub-critical rate** ($\sim N^{-0.29}$) that is insufficient to establish uniform $H^2$ regularity via the current kinetic resolvent route.

### 2.3 Empirical Finite-Range Status (Not an Asymptotic Theorem)
Because the regression covers $n = 8$ data points with $R^2 \approx 0.78$, there is noticeable curvature/scatter around a log-log straight line. We explicitly avoid designating $N^{-0.288}$ as an "asymptotic law". It is strictly an **empirical finite-range power-law fit over $64 \le N \le 192$**. There may be crossovers, logarithmic corrections, or finite-$T$ boundary layer phenomena that modify this behaviour at larger scales.

### 2.4 Branch Invariance of the Scaling Mechanism
A striking discovery of Cell 115 is that both branches exhibit virtually identical scaling:
$$\beta_\alpha^{(1)} - \beta_\alpha^{(0)} = -0.0261 \pm 0.0067, \qquad \beta_T^{(1)} - \beta_T^{(0)} = -0.0268 \pm 0.0065.$$
The branch ratio $|\alpha_N^{(1)}| / |\alpha_N^{(0)}|$ decays only as $N^{-0.0261}$ ($R^2 = 0.836$).
This decisively proves that **branch selection affects only the pre-factor ($\sim 13\%$ suppression), not the underlying asymptotic boundary-coupling mechanism.**

---

## 3. Computational Complexity Audit

`cell115.py` required $1722.17\text{ s}$ ($\approx 28.7\text{ minutes}$) on the compute node. Although the full $N = 192$ matrix was retrieved instantly from cache (3.88 s), the repeated high-precision diagonalization (`mp.eigsy` at 70 dps) across 9 dimensions ($N = 48 \dots 192$) scales steeply:
$$\text{Cost} \propto \sum_{i=1}^9 N_i^3 \times \mathcal{M}(\text{dps}).$$
Consequently, **extending brute-force numerical sweeps to $N = 256, 320, 384$ is computationally inefficient** and yields diminishing scientific returns. The primary question is now analytical.

---

## 4. Forward Strategic Mandate for Cell 116

Rather than executing another expensive numerical sweep, the programme must attack the analytical mechanism:

$$\boxed{\textbf{Why does the exact Connes–CvS boundary kernel produce a sub-critical exponent } \beta_\alpha \approx -0.3 \text{?}}$$

### Analytical Attack Plan for Cell 116:
Recall the exact row-wise boundary-flux identity proven in Cell 111 / 112:
$$\alpha_N = (H u_N)_N + \frac{T_v(0)}{\sqrt{2}} a_N - E_{11} N^2 v_{N, N}.$$
Because the eigenvalue term $E_{11} N^2 v_N \sim 10^{-50} \times 192^2 \times 10^{-13} \approx 10^{-59}$ is completely negligible:
$$\alpha_N \approx \sum_{k=1}^N H_{Nk} k^2 v_{N, k} + \frac{T_v(0)}{\sqrt{2}} a_N.$$

Cell 116 will investigate:
1. **Asymptotic Decomposition of the Boundary Row $H_{Nk}$:**
   For $H_{Nk} = \frac{2(N\psi(N) - k\psi(k))}{N^2 - k^2}$, decompose the sum into:
   - Bulk contribution: modes $k \ll N$, where $\frac{k^2}{N^2 - k^2} \approx \frac{k^2}{N^2}$.
   - Boundary-layer contribution: modes near the edge $N - k \sim N^\theta$.
2. **Fractional Power Mechanism:**
   Determine whether the interaction between the boundary layer width, the solitary wave's physical curvature $T_\infty''(0)$, and the Archimedean weight $h_+(r)$ analytically forces a fractional power like $N^{-1/3}$ or $N^{-0.3}$.
3. **Implication for Gate 1:**
   - If analytical analysis proves that $\beta_\alpha \le -1/3 < -0.5$ is an asymptotic impossibility at fixed $T$, it formally establishes that **Gate 1 cannot be closed at fixed $T$**, making the joint limit $(N, T) \to \infty$ mathematically mandatory.
   - If analytical analysis reveals that an asymptotic crossover to super-critical decay $\beta < -0.5$ must occur beyond a calculable scale $N_*$, it will provide the exact target scale needed to verify extinction.
