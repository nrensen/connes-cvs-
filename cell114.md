# Cell 114 Analytical Note: Two-State Spectral Reordering Anatomy & Localization Invariants

**Companion Computational Script:** [`cell114.py`](file:///c:/data/github/connes-cvs-/cell114.py) | **Verification Log:** [`cell114.out`](file:///c:/data/github/connes-cvs-/cell114.out)  
**Status:** Executed & Audited (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell113.md`](file:///c:/data/github/connes-cvs-/cell113.md); [`cell113.out`](file:///c:/data/github/connes-cvs-/cell113.out); [`cell112a.md`](file:///c:/data/github/connes-cvs-/cell112a.md); [`cell112a.out`](file:///c:/data/github/connes-cvs-/cell112a.out); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Headline Verdict

Cell 114 was designed to resolve the low-energy spectral reordering diagnosed in Cells 112a and 113. The calculation executed on external compute (`cell114.out`, runtime 1048.25 s at 70 dps) tracks the two-dimensional low-energy subspace across the transition region $N \in [40, 60]$ ($\Delta N = 2$), audits physical localization invariants, compares boundary coupling observables side-by-side across $N \in \{64, 80, 96, 128, 192\}$, and measures finite-$T$ sensitivity at $N=48$ across $T \in \{400, 500, 600\}$.

### Headline Verdict

> **1. The branch ambiguity has effectively been resolved numerically:** For $N \ge 64$, the $k = 1$ branch is overwhelmingly the natural localized continuation, characterized by a stationary central amplitude $v_0(k=1) \to 0.6664$ and stable energy $E_1 \to 4.67 \times 10^{-50}$. Conversely, the $k = 0$ state is a delocalized edge/background mode with $v_0(k=0) \to 0.0638$ and energy sinking into the negative Archimedean leakage floor $E_0 \to -1.063 \times 10^{-51}$.
>
> **2. Crucial Negative Finding — The Branch Escape Route Has Failed:** On the localized $k = 1$ branch, the boundary extinction product $P_\alpha(k=1) \equiv |\alpha_N| \sqrt{N}$ does **NOT** decay. It remains locked in essentially the same $10^{-21}$ regime as the $k = 0$ state:
> $$P_\alpha^{(1)}(64) = 1.938 \times 10^{-21} \quad \longrightarrow \quad P_\alpha^{(1)}(192) = 2.328 \times 10^{-21}.$$
> Thus the hoped-for explanation *"the $2 \times 10^{-21}$ extinction plateau was merely an artifact of tracking the wrong spectral branch"* is **firmly refuted**.

---

## 1. Quantitative High-$N$ Branch Comparison & Extinction Audit

Cell 114 evaluated both low-energy states concurrently across $N \in \{64, 80, 96, 128, 192\}$ at $c = 13, T = 600$:

| $N$ | Branch | Energy $E$ | $v_0$ | $\mathcal{K}_2$ | $|T_v(0)|$ | $|\alpha_N|$ | $P_\alpha = |\alpha_N|\sqrt{N}$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $k=0$ (Edge) | $-5.085 \times 10^{-52}$ | $0.0708$ | $82.19$ | $9.674 \times 10^{-25}$ | $2.684 \times 10^{-22}$ | $2.147 \times 10^{-21}$ |
| | $k=1$ (Localized) | $4.713 \times 10^{-50}$ | $0.6657$ | $25.46$ | $8.721 \times 10^{-25}$ | $2.422 \times 10^{-22}$ | $1.938 \times 10^{-21}$ |
| **80** | $k=0$ (Edge) | $-8.388 \times 10^{-52}$ | $0.0666$ | $82.72$ | $7.921 \times 10^{-25}$ | $2.222 \times 10^{-22}$ | $1.987 \times 10^{-21}$ |
| | $k=1$ (Localized) | $4.686 \times 10^{-50}$ | $0.6661$ | $24.94$ | $7.017 \times 10^{-25}$ | $1.971 \times 10^{-22}$ | $1.762 \times 10^{-21}$ |
| **96** | $k=0$ (Edge) | $-9.335 \times 10^{-52}$ | $0.0654$ | $82.87$ | $7.299 \times 10^{-25}$ | $2.060 \times 10^{-22}$ | $2.018 \times 10^{-21}$ |
| | $k=1$ (Localized) | $4.678 \times 10^{-50}$ | $0.6662$ | $24.80$ | $6.427 \times 10^{-25}$ | $1.816 \times 10^{-22}$ | $1.780 \times 10^{-21}$ |
| **128** | $k=0$ (Edge) | $-1.005 \times 10^{-51}$ | $0.0645$ | $82.98$ | $6.922 \times 10^{-25}$ | $1.970 \times 10^{-22}$ | $2.229 \times 10^{-21}$ |
| | $k=1$ (Localized) | $4.673 \times 10^{-50}$ | $0.6663$ | $24.68$ | $6.065 \times 10^{-25}$ | $1.729 \times 10^{-22}$ | $1.956 \times 10^{-21}$ |
| **192** | $k=0$ (Edge) | $-1.063 \times 10^{-51}$ | $0.0638$ | $83.06$ | $6.663 \times 10^{-25}$ | $1.923 \times 10^{-22}$ | $2.664 \times 10^{-21}$ |
| | $k=1$ (Localized) | $4.668 \times 10^{-50}$ | $0.6664$ | $24.59$ | $5.812 \times 10^{-25}$ | $1.680 \times 10^{-22}$ | $2.328 \times 10^{-21}$ |

### 1.1 The Boundary Coupling Ratio
The ratio of boundary coupling between the localized state and the edge state:
$$\frac{P_\alpha^{(1)}}{P_\alpha^{(0)}} = \frac{|\alpha_N^{(1)}|}{|\alpha_N^{(0)}|} = 0.9022 \; (N=64) \;\longrightarrow\; 0.8870 \; (80) \;\longrightarrow\; 0.8818 \; (96) \;\longrightarrow\; 0.8775 \; (128) \;\longrightarrow\; 0.8737 \; (192).$$
The localized branch has a consistently smaller boundary coupling ($\sim 13\%$ suppression relative to $k=0$). However, this mild suppression is **not asymptotic extinction**: $P_\alpha^{(1)}(192) = 2.328 \times 10^{-21}$ is actually slightly *larger* than $P_\alpha^{(1)}(64) = 1.938 \times 10^{-21}$.

### 1.2 The Contact Defect $|T_v(0)|$
For $k=1$, the physical contact defect exhibits a very slow decline:
$$|T_v(0)|: \quad 8.72 \times 10^{-25} \;\to\; 7.02 \times 10^{-25} \;\to\; 6.43 \times 10^{-25} \;\to\; 6.06 \times 10^{-25} \;\to\; 5.81 \times 10^{-25}.$$
This slow decline does not approach the required $o(N^{-3/2})$ rate. In fact, scaling by $N^{3/2}$:
$$N^{3/2}|T_v(0)|: \quad 4.46 \times 10^{-22} \;\longrightarrow\; 1.54 \times 10^{-21},$$
which increases with $N$.

### 1.3 Epistemic Implication for Gate 1
The two sufficient conditions identified in Cell 110 for closing Gate 1:
$$\alpha_N = o(N^{-1/2}) \qquad \text{and} \qquad T_v(0) = o(N^{-3/2})$$
are **not supported numerically at the observed scale ($T = 600, N \le 192$)**.

Consequently, Gate 1 must now be split epistemically into:
$$\boxed{\text{Regularity Mechanism (Theorem 9.16): Conditionally Closed}}$$
$$\boxed{\text{Boundary-Defect Extinction: Currently Unsupported by Numerical Data}}$$

---

## 2. Analytical & Diagnostic Corrections to Cell 114

Four specific analytical, diagnostic, and conceptual corrections must be recorded to maintain audit-proof integrity:

### Correction 1: Re-Interpretation of the $\mathcal{K}_2$ Kinetic Sobolev Invariant
Prior to Cell 114, working notes (Cells 106, 110) hypothesized that $\mathcal{K}_2 \approx 83.1$ was the characteristic invariant of the continuum solitary wave.
The numerical data decisively refute this assignment:
$$\mathcal{K}_2(k=0) \longrightarrow 83.06, \qquad \mathcal{K}_2(k=1) \longrightarrow 24.59.$$
The $\mathcal{K}_2 \approx 83.1$ plateau observed in earlier large-$N$ runs belongs entirely to the **edge/background branch ($k=0$)**, not the localized solitary wave ($k=1$).

Furthermore, the fine sweep across $N \in [40, 60]$ reveals that the two states **exchange their $\mathcal{K}_2$ character**:

| $N$ | $\mathcal{K}_2(k=0)$ | $\mathcal{K}_2(k=1)$ | $v_0(k=0)$ | $v_0(k=1)$ |
| :---: | :---: | :---: | :---: | :---: |
| **40** | $9.23$ | $97.14$ | $0.5410$ | $0.3963$ |
| **44** | $9.98$ | $96.96$ | $0.5350$ | $0.4036$ |
| **48** | $21.45$ | $85.81$ | $0.4558$ | $0.4908$ |
| **50** | $44.97$ | $62.39$ | $0.3149$ | $0.5910$ |
| **52** | $69.40$ | $38.09$ | $0.1629$ | $0.6495$ |
| **54** | $76.91$ | $30.63$ | $0.1104$ | $0.6604$ |
| **56** | $79.92$ | $27.67$ | $0.0881$ | $0.6637$ |
| **60** | $81.46$ | $26.17$ | $0.0764$ | $0.6651$ |

As $v_0(k=0)$ collapses from $0.541$ to $0.076$, its kinetic moment $\mathcal{K}_2$ climbs from $9.23$ to $81.46$. Simultaneously, $v_0(k=1)$ climbs from $0.396$ to $0.665$, and its kinetic moment drops from $97.14$ to $26.17$.
Because the physical meaning of $\mathcal{K}_2$ differs substantially between branches, $\mathcal{K}_2$ is **withdrawn as a branch discriminator**. Branch identification is governed by the central amplitude $v_0 \approx 0.666$ and energy stability $E_1 \approx 4.67 \times 10^{-50}$.

### Correction 2: Re-Interpretation of the SVD Projection Diagnostic $\phi$
Cell 114 defined:
$$M = V_A^T V_B, \qquad \phi \equiv \arctan\left(\frac{|M_{01}|}{|M_{00}|}\right).$$
This $\phi$ is **not** an internal rotation angle between the two eigenstates at a fixed dimension $N$. Rather, it is a **cross-$N$ basis-overlap diagnostic** $\phi_{\mathrm{overlap}}$ measuring how the first basis vector at dimension $N_A$ projects onto the two eigenvectors at dimension $N_B$.

The numerical data:
$$\phi_{\mathrm{overlap}}(40 \to 42) = 0.22^\circ, \quad \phi_{\mathrm{overlap}}(48 \to 50) = 14.85^\circ, \quad \phi_{\mathrm{overlap}}(50 \to 52) = 13.99^\circ, \quad \phi_{\mathrm{overlap}}(58 \to 60) = 0.72^\circ.$$
The angle peaks modestly near $15^\circ$ during the state-mixing window and returns to near zero. It does **not** reach $45^\circ$ or $90^\circ$.
Crucially, the principal singular values satisfy:
$$\sigma_1 \approx 1.0000000, \qquad \sigma_2 \ge 0.9999986 \quad (\forall N \in [40, 60]).$$
This confirms that the **two-dimensional low-energy subspace $\operatorname{span}\{v_0, v_1\}$ is remarkably invariant under $N$-refinement**, while the eigenbasis undergoes a localized rotation near $N \approx 50$.

### Correction 3: Scientific Precision for Core Mass Tails
Reporting $L_{24} = 1.0$ in terminal tables at fixed decimal places obscured the modal tail leakage. In future cells and formal records, the complement:
$$1 - L_{24} \equiv \sum_{m=25}^N v_m^2$$
must be reported in scientific notation to distinguish true exponential localization from low-level polynomial tails.

### Correction 4: Epistemic Status of Branch Identification
We avoid designating $k = 1$ as a "certified continuum solitary branch", as no continuum convergence theorem is yet proved. It is rigorously designated as the **candidate localized/solitary branch for $N \ge 64$**, characterized by the empirical signature:
$$v_0 \approx 0.6664, \qquad E \approx 4.67 \times 10^{-50}.$$

---

## 3. Finite-$T$ Sensitivity at Fixed $N = 48$ (Module 4)

Evaluating both states across cutoffs $T \in \{400, 500, 600\}$ at $N = 48$:

| $T$ | $k$ | Energy $E$ | $v_0$ | $\mathcal{K}_2$ | $|T_v(0)|$ | $|\alpha_{48}|$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **400** | $0$ | $8.295 \times 10^{-57}$ | $0.5404$ | $9.31$ | $2.687 \times 10^{-27}$ | $6.168 \times 10^{-25}$ |
| | $1$ | $1.119 \times 10^{-49}$ | $0.3958$ | $97.94$ | $1.008 \times 10^{-23}$ | $2.168 \times 10^{-21}$ |
| **500** | $0$ | $8.530 \times 10^{-57}$ | $0.5404$ | $9.31$ | $2.831 \times 10^{-27}$ | $7.169 \times 10^{-25}$ |
| | $1$ | $1.149 \times 10^{-49}$ | $0.3958$ | $97.95$ | $9.593 \times 10^{-24}$ | $2.288 \times 10^{-21}$ |
| **600** | $0$ | $2.103 \times 10^{-50}$ | $0.4558$ | $21.45$ | $1.783 \times 10^{-24}$ | $4.570 \times 10^{-22}$ |
| | $1$ | $1.490 \times 10^{-49}$ | $0.4908$ | $85.81$ | $9.270 \times 10^{-24}$ | $2.377 \times 10^{-21}$ |

### Physical Interpretation:
At $T = 400$ and $T = 500$, the state $k = 0$ is unperturbed ($v_0 = 0.5404, \mathcal{K}_2 = 9.31, E \sim 10^{-56}$). At $T = 600$, the energy of $k=0$ jumps to $2.10 \times 10^{-50}$ and $v_0$ drops to $0.4558$, while $k=1$ shifts to $v_0 = 0.4908$.
This demonstrates that **increasing the Archimedean cutoff $T$ moves the system into the exact same state-mixing regime as increasing $N$**. The spectral reordering is intrinsically tied to the interaction between the high-$T$ Archimedean integral cutoff and the finite-$N$ Fourier Galerkin grid.

---

## 4. Forward Research Path & Cell 115 Mandate

With the branch-selection problem settled ($k = 1$ is the localized branch for $N \ge 64$), the research programme must confront the core question without assuming in advance that the boundary defect must vanish:

### The Asymptotic Scaling Question:
Does the localized branch possess a **nonzero limiting boundary defect** $\alpha_\infty > 0$, or is there a genuine, slow power-law decay?

### Cell 115 Plan:
1. Evaluate $N \in \{48, 64, 80, 96, 112, 128, 144, 160, 192\}$ on the cached $T = 600$ matrix.
2. For $k = 1$, compute:
   $$|\alpha_N|, \quad |T_v(0)|, \quad N^{1/2}|\alpha_N|, \quad N|\alpha_N|, \quad N^{3/2}|T_v(0)|, \quad N|T_v(0)|.$$
3. Perform logarithmic power-law regressions:
   $$\log_{10} |\alpha_N| = \beta_\alpha \log_{10} N + C_\alpha,$$
   $$\log_{10} |T_v(0)| = \beta_T \log_{10} N + C_T,$$
   $$\log_{10} P_\alpha(N) = \gamma_\alpha \log_{10} N + D_\alpha.$$
4. If $\beta_\alpha \approx 0$ (or $\gamma_\alpha \approx 0.5$), the boundary defect $\alpha_N$ does not vanish, confirming a nonzero limiting defect $\alpha_\infty > 0$ for the finite-$T=600$ Galerkin model. This would formally falsify naive finite-$T$ boundary extinction and guide the required formulation of the joint $(N, T) \to \infty$ limit in Gate 1.
