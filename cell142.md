# CELL 142 — Precision Robustness, Extended Well Spectrum, and Energy-Deficit Participation Rank $r_\eta$

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction)  
**Milestone:** Milestone M-G1.6 (Variational Lower Bound & Coupled Operator Geometry)  
**Target Operators:** Certified Continuum Truncations $K_{\mathrm{rest}}(N)$, $W_\perp(N)$, and $H_1(N) = K_{\mathrm{rest}} - W_\perp$ on $\mathcal{B}_{11}^\perp$  
**Pre-Flight Invariants ($N = 64$):** $\omega_0 = 2.9315259531$, $\nu_0 = 4.2604954421$, $\mu_0 = -0.4869792197$  
**Computational Output:** Certified via [`cell142.out`](file:///c:/data/github/connes-cvs-/cell142.out) (runtime: 328.89s)  

---

## 1. Executive Summary & Epistemic Synthesis

Cell 142 delivers a decisive multi-precision audit of the projected prime well operator $W_\perp \equiv P_{\mathcal{B}_{11}^\perp} \widetilde{W} P_{\mathcal{B}_{11}^\perp}$, surveys the first 20 eigenvalues of its continuum spectrum, certifies the modal deficit conservation identity, and evaluates the energy-deficit participation rank $r_\eta(N)$.

### 1.1 Summary of Established Results

1. **Precision Stability Certified (Scenario C Refuted at Finite $N$):**  
   Across dimensions $N \in \{48, 56, 64\}$, the top five splittings $\delta\nu_1 \dots \delta\nu_5$ agree between 50 dps and 70 dps to the full displayed precision (relative discrepancy $\rho_j = 0.00 \times 10^{-0}$). The observed $10^{-8}$-scale splitting $\delta\nu_1(64) = 3.7391 \times 10^{-8}$ is a **numerically stable eigenvalue splitting of the specified finite-$N$ projected operator under $50 \to 70$ dps precision escalation**, not a numerical precision artefact. (This does not rule out modeling/truncation phenomena as $N \to \infty$).
2. **Extended Top-20 Spectrum Profile:**  
   The spectrum at $N=64$ ($q = 54$) exhibits a highly compressed top region ($j = 0, 1, 2, 3$), followed by a transitional sector ($j = 4, 5, 6$), merging into a macroscopic bulk ($j \ge 7$). Successive splitting ratios $\delta\nu_j / \delta\nu_{j-1}$ drop below the diagnostic threshold $2.0$ at index $j_{\mathrm{cross}}^{(2)} = 7$ ($\delta\nu_8/\delta\nu_7 = 1.3118$).
3. **Exact Modal Deficit Conservation Certified:**  
   The physical coordinate well deficit $\Delta W_{\mathrm{phys}} = \nu_0 - v(1)^T W_\perp v(1)$ and the spectral modal sum $\Delta W_{\mathrm{spec}} = \sum_{j=1}^{q-1} \delta\nu_j P_W(j)$ agree to the precision floor:
   $$\Delta W_{\mathrm{phys}} = 0.5507123831, \quad \Delta W_{\mathrm{spec}} = 0.5507123831, \quad \mathcal{R}_{\mathrm{cons}} = 0.00 \times 10^{-0}.$$
4. **Decoupling of State Occupation from Energy Sacrifice:**  
   Modes $j \in \{0, 1, 2, 3\}$ contain **$49.2933\%$ of the ground-state probability mass**, yet contribute only $1.087 \times 10^{-5}$ ($0.00197\%$) to the well-energy deficit $\Delta W(1)$. Outside the top four modes ($j \ge 4$), over **$99.998\%$** of the well sacrifice is paid.
5. **Decisive Falsification of Low-Dimensional Reduction (Celebrated Negative Result):**  
   The energy-deficit participation rank $r_{0.05}(N)$ (modes required to account for $95\%$ of $\Delta W(1)$) scales as approximately **$90\%$ of the entire continuum dimension $q(N)$**:
   $$\frac{r_{0.05}(N)}{q(N)} \in \{0.909, \; 0.895, \; 0.913, \; 0.907\} \approx 0.90.$$
   The data **decisively falsify** the hoped-for low-dimensional active subspace reduction ($r_\eta \ll q$). The top cluster is **not** an energetic active subspace; rather, the system exhibits a dual structure:
   $$\boxed{\text{Top Cluster } (j \le 3) = \text{low-cost reservoir of probability mass}}$$
   $$\boxed{\text{Deep Spectral Tail } (j \ge 4) = \text{distributed reservoir paying the energy deficit}}.$$

---

## 2. Mathematical Framework & Calibrated Definitions

### 2.1 Multi-Precision Invariance & Scenario C Protocol

Let $W_\perp^{(d)}(N)$ denote the projected potential well operator evaluated at numerical precision $d \in \{50, 70\}$ decimal digits. Let $\nu_j^{(d)}(N)$ be its eigenvalues in descending order:
$$\nu_0^{(d)}(N) > \nu_1^{(d)}(N) > \dots > \nu_{q-1}^{(d)}(N), \qquad \delta\nu_j^{(d)}(N) \equiv \nu_0^{(d)}(N) - \nu_j^{(d)}(N).$$

> **Proposition 142.1 (Multi-Precision Stability).**  
> Across tested dimensions $N \in \{48, 56, 64\}$ and modes $j \in \{1, \dots, 5\}$, the relative discrepancy satisfies:
> $$\rho_j(N) \equiv \frac{\left|\delta\nu_j^{(50)}(N) - \delta\nu_j^{(70)}(N)\right|}{\delta\nu_j^{(70)}(N)} = 0.00 \times 10^{-0}.$$
> This confirms that the observed $10^{-8}$-scale splittings are numerically stable eigenvalue splittings of the specified finite-$N$ projected operator under $50 \to 70$ dps precision escalation.

---

### 2.2 The Extended Spectrum Profile & Diagnostic Crossing Index $j_{\mathrm{cross}}^{(2)}$

Let $\gamma_j \equiv \delta\nu_j / \delta\nu_{j-1}$ denote the successive splitting ratio for $j \ge 2$.

> **Definition 142.2 (Diagnostic Ratio-Crossing Index $j_{\mathrm{cross}}^{(2)}$).**  
> The diagnostic ratio-crossing index under threshold $2.0$ is defined as:
> $$j_{\mathrm{cross}}^{(2)} \equiv \min \left\{ j \ge 2 : \frac{\delta\nu_j}{\delta\nu_{j-1}} < 2.0 \right\}.$$
> At $N = 64$, the sequence of ratios begins:
> $$\gamma_2 = 178.79, \quad \gamma_3 = 19.22, \quad \gamma_4 = 59.38, \quad \gamma_5 = 9.19, \quad \gamma_6 = 4.55, \quad \gamma_7 = 1.69.$$
> Consequently, $j_{\mathrm{cross}}^{(2)} = 7$ (since $\delta\nu_7 / \delta\nu_6 = 1.6858 < 2$).  
> *Calibration note:* $j_{\mathrm{cross}}^{(2)}$ is an empirical diagnostic indicator of the transition toward bulk mode spacing; it is **not** an intrinsic mathematical cluster boundary and does **not** define the active subspace.

---

### 2.3 Modal Deficit Conservation Identity

The physical ground state $v(1)$ of $H(1) = K_{\mathrm{rest}} - W_\perp$ induces the discrete modal projection weights:
$$P_W(j) \equiv |\langle y_j, v(1) \rangle|^2, \qquad \sum_{j=0}^{q-1} P_W(j) = 1.$$
The well-energy deficit decomposes exactly as:
$$\Delta W(1) \equiv \nu_0 - \langle v(1), W_\perp v(1) \rangle = \sum_{j=1}^{q-1} \delta\nu_j P_W(j) \equiv \sum_{j=1}^{q-1} \Delta W_j.$$

> **Proposition 142.3 (Modal Deficit Identity Certified).**  
> The physical coordinate deficit $\Delta W_{\mathrm{phys}}$ and the spectral modal sum $\Delta W_{\mathrm{spec}}$ satisfy:
> $$\mathcal{R}_{\mathrm{cons}} \equiv \left| \Delta W_{\mathrm{phys}} - \sum_{j=1}^{q-1} \delta\nu_j P_W(j) \right| = 0.00 \times 10^{-0},$$
> certified to the precision floor across all tested dimensions.

---

### 2.4 Energy-Deficit Participation Rank $r_\eta(N)$

Rather than imposing an arbitrary spectral cutoff $\epsilon$, the distribution of well-energy deficit across modes is tracked by the energy-deficit participation rank.

> **Definition 142.4 (Energy-Deficit Participation Rank $r_\eta$).**  
> For any deficit tolerance fraction $\eta \in (0, 1)$, the energy-deficit participation rank $r_\eta(N)$ is the minimal number of top-ordered $W$-modes required to account for at least a $(1 - \eta)$ fraction of the total well sacrifice $\Delta W(1)$:
> $$r_\eta(N) \equiv \min \left\{ r \ge 1 : \frac{\sum_{j=0}^{r-1} \delta\nu_j(N) P_W(j; N)}{\Delta W(1; N)} \ge 1 - \eta \right\}.$$
> *Epistemic scope:* $r_\eta$ is a diagnostic of the spectral breadth of the well-energy sacrifice in the $W$-eigenbasis; it does **not** establish a low-dimensional variational reduction, nor does it measure the dimension of the state in a basis-independent sense.

---

## 3. Certified Computational Results (`cell142.out`)

### 3.1 Table 1: Dual-Precision Comparison (50 dps vs 70 dps)
Pre-flight audit at $N=64$ passed in 97.24s ($\omega_0 = 2.9315259531$, $\nu_0 = 4.2604954421$, $\mu_0 = -0.4869792197$; residuals $< 1.1 \times 10^{-7}$).

| $N$ | Mode $j$ | $\delta\nu_j$ (50 dps) | $\delta\nu_j$ (70 dps) | Relative Discrepancy $\rho_j$ | Agreement Digits | Status |
|---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 48 | 1 | $3.135858546630 \times 10^{-5}$ | $3.135858546630 \times 10^{-5}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 48 | 2 | $4.825916646516 \times 10^{-3}$ | $4.825916646516 \times 10^{-3}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 48 | 3 | $6.599076477002 \times 10^{-2}$ | $6.599076477002 \times 10^{-2}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 48 | 4 | $3.613086866983 \times 10^{-1}$ | $3.613086866983 \times 10^{-1}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 48 | 5 | $5.614516012888 \times 10^{-1}$ | $5.614516012888 \times 10^{-1}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 56 | 1 | $2.024887509137 \times 10^{-6}$ | $2.024887509137 \times 10^{-6}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 56 | 2 | $8.224136363627 \times 10^{-5}$ | $8.224136363627 \times 10^{-5}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 56 | 3 | $6.485978834467 \times 10^{-3}$ | $6.485978834467 \times 10^{-3}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 56 | 4 | $6.932862018960 \times 10^{-2}$ | $6.932862018960 \times 10^{-2}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 56 | 5 | $3.378588789026 \times 10^{-1}$ | $3.378588789026 \times 10^{-1}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 64 | 1 | $3.739144503596 \times 10^{-8}$ | $3.739144503596 \times 10^{-8}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 64 | 2 | $6.685161220052 \times 10^{-6}$ | $6.685161220052 \times 10^{-6}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 64 | 3 | $1.285020586882 \times 10^{-4}$ | $1.285020586882 \times 10^{-4}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 64 | 4 | $7.630214382405 \times 10^{-3}$ | $7.630214382405 \times 10^{-3}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |
| 64 | 5 | $7.012847752789 \times 10^{-2}$ | $7.012847752789 \times 10^{-2}$ | $0.0000 \times 10^{-0}$ | 50.0 | **CERTIFIED** |

---

### 3.2 Table 2: Extended Top-20 Spectrum Splittings at $N = 64$ ($q = 54$)

| Mode $j$ | Eigenvalue $\nu_j$ | Splitting $\delta\nu_j$ | Successive Ratio $\delta\nu_j / \delta\nu_{j-1}$ | Diagnostic Classification |
|---:|:---:|:---:|:---:|:---:|
| 0 | 4.2604954421 | $0.000000000000 \times 10^{-0}$ | — | Top Well Mode ($y_0$) |
| 1 | 4.2604954047 | $3.739144503596 \times 10^{-8}$ | $\infty$ | Compressed Cluster |
| 2 | 4.2604887569 | $6.685161220052 \times 10^{-6}$ | 178.7885 | Compressed Cluster |
| 3 | 4.2603669400 | $1.285020586882 \times 10^{-4}$ | 19.2220 | Compressed Cluster |
| 4 | 4.2528652277 | $7.630214382405 \times 10^{-3}$ | 59.3781 | Transition / Crossover |
| 5 | 4.1903669646 | $7.012847752789 \times 10^{-2}$ | 9.1909 | Transition / Crossover |
| 6 | 3.9412059005 | $3.192895415718 \times 10^{-1}$ | 4.5529 | Transition / Crossover |
| 7 | 3.7222312748 | $5.382641673274 \times 10^{-1}$ | 1.6858 | Macroscopic Bulk ($j_{\mathrm{cross}}^{(2)} = 7$) |
| 8 | 3.5543793988 | $7.061160433045 \times 10^{-1}$ | 1.3118 | Macroscopic Bulk |
| 9 | 3.5405949516 | $7.199004904627 \times 10^{-1}$ | 1.0195 | Macroscopic Bulk |
| 10 | 3.5375133800 | $7.229820621203 \times 10^{-1}$ | 1.0043 | Macroscopic Bulk |
| 11 | 3.5366923803 | $7.238030618054 \times 10^{-1}$ | 1.0011 | Macroscopic Bulk |
| 12 | 3.5279059533 | $7.325894887960 \times 10^{-1}$ | 1.0121 | Macroscopic Bulk |
| 13 | 3.4630674695 | $7.974279726061 \times 10^{-1}$ | 1.0885 | Macroscopic Bulk |
| 14 | 3.3792693083 | $8.812261338446 \times 10^{-1}$ | 1.1051 | Macroscopic Bulk |
| 15 | 3.2561612483 | $1.004334193854 \times 10^{0}$ | 1.1397 | Macroscopic Bulk |
| 16 | 3.1853244769 | $1.075170965217 \times 10^{0}$ | 1.0705 | Macroscopic Bulk |
| 17 | 3.1713276775 | $1.089167764604 \times 10^{0}$ | 1.0130 | Macroscopic Bulk |
| 18 | 3.1488784015 | $1.111617040634 \times 10^{0}$ | 1.0206 | Macroscopic Bulk |
| 19 | 3.0772224526 | $1.183272989553 \times 10^{0}$ | 1.0645 | Macroscopic Bulk |

---

### 3.3 Table 3: Joint Spectral/Weight Profile & Modal Deficit Flow ($N = 64$)

| Mode $j$ | Splitting $\delta\nu_j$ | Modal Mass $P_W(j)$ | Cum Mass $\Pi(j)$ | Deficit $\Delta W_j$ | Cum Deficit | Deficit Share $\%$ |
|---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | $0.0000 \times 10^{-0}$ | 21.0581% | 21.0581% | $0.0000 \times 10^{-0}$ | 0.00000000 | 0.00% |
| 1 | $3.7391 \times 10^{-8}$ | 10.4038% | 31.4619% | $3.8901 \times 10^{-9}$ | 0.00000000 | 0.00% |
| 2 | $6.6852 \times 10^{-6}$ | 9.8925% | 41.3543% | $6.6133 \times 10^{-7}$ | 0.00000067 | 0.00% |
| 3 | $1.2850 \times 10^{-4}$ | 7.9390% | **49.2933%** | $1.0202 \times 10^{-5}$ | **0.00001087** | **0.00197%** |
| 4 | $7.6302 \times 10^{-3}$ | 3.6702% | 52.9635% | $2.8004 \times 10^{-4}$ | 0.00029091 | 0.05% |
| 5 | $7.0128 \times 10^{-2}$ | 6.0328% | 58.9962% | $4.2307 \times 10^{-3}$ | 0.00452159 | 0.82% |
| 6 | $3.1929 \times 10^{-1}$ | 0.0439% | 59.0401% | $1.4010 \times 10^{-4}$ | 0.00466170 | 0.85% |
| 7 | $5.3826 \times 10^{-1}$ | 0.2001% | 59.2402% | $1.0773 \times 10^{-3}$ | 0.00573895 | 1.04% |
| 8 | $7.0612 \times 10^{-1}$ | 0.7985% | 60.0388% | $5.6385 \times 10^{-3}$ | 0.01137740 | 2.07% |
| 9 | $7.1990 \times 10^{-1}$ | 1.1338% | 61.1725% | $8.1619 \times 10^{-3}$ | 0.01953934 | 3.55% |
| 10 | $7.2298 \times 10^{-1}$ | 0.0127% | 61.1853% | $9.2119 \times 10^{-5}$ | 0.01963146 | 3.56% |
| 11 | $7.2380 \times 10^{-1}$ | 0.0000% | 61.1853% | $2.4650 \times 10^{-8}$ | 0.01963149 | 3.56% |
| 12 | $7.3259 \times 10^{-1}$ | 0.0005% | 61.1857% | $3.4086 \times 10^{-6}$ | 0.01963490 | 3.57% |
| 13 | $7.9743 \times 10^{-1}$ | 4.8846% | 66.0704% | $3.8951 \times 10^{-2}$ | 0.05858636 | 10.64% |
| 14 | $8.8123 \times 10^{-1}$ | 2.7950% | 68.8654% | $2.4630 \times 10^{-2}$ | 0.08321649 | 15.11% |
| 15 | $1.0043 \times 10^{0}$ | 1.5655% | 70.4308% | $1.5723 \times 10^{-2}$ | 0.09893930 | 17.97% |
| 16 | $1.0752 \times 10^{0}$ | 1.1317% | 71.5626% | $1.2168 \times 10^{-2}$ | 0.11110704 | 20.18% |
| 17 | $1.0892 \times 10^{0}$ | 0.3459% | 71.9085% | $3.7675 \times 10^{-3}$ | 0.11487454 | 20.86% |
| 18 | $1.1116 \times 10^{0}$ | 3.5753% | 75.4838% | $3.9744 \times 10^{-2}$ | 0.15461852 | 28.08% |
| 19 | $1.1833 \times 10^{0}$ | 2.0220% | 77.5058% | $2.3926 \times 10^{-2}$ | 0.17854453 | 32.42% |

*Total Physical Well Deficit:* $\Delta W_{\mathrm{phys}} = 0.5507123831$  
*Total Spectral Modal Deficit:* $\Delta W_{\mathrm{spec}} = 0.5507123831$  
*Conservation Residual:* $\mathcal{R}_{\mathrm{cons}} = 0.00 \times 10^{-0}$  
*Cumulative probability mass on $j \le 19$:* $77.51\%$, paying only $32.42\%$ of the deficit! The remaining $67.58\%$ of the deficit is paid across modes $j \in \{20, \dots, 53\}$.

---

### 3.4 Table 4: Energy-Deficit Participation Rank $r_\eta(N)$

| $N$ | $q$ | $\Delta W(1)$ | $r_{0.50}$ (50% Deficit) | $r_{0.25}$ (75% Deficit) | $r_{0.10}$ (90% Deficit) | $r_{0.05}$ (95% Deficit) | $r_{0.01}$ (99% Deficit) | $\mathbf{r_{0.05} / q}$ |
|---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 32 | 22 | 0.570807 | 11 | 17 | 19 | 20 | 22 | **0.909** |
| 48 | 38 | 0.553818 | 18 | 30 | 33 | 34 | 38 | **0.895** |
| 56 | 46 | 0.552398 | 22 | 35 | 40 | 42 | 45 | **0.913** |
| 64 | 54 | 0.550712 | 26 | 42 | 44 | 49 | 52 | **0.907** |

---

## 4. Analytical Interpretation & Falsification Analysis

### 4.1 Falsification of Low-Dimensional Reduction
The central strategic hypothesis of Cell 142 was that the energy-weighted deficit would be concentrated in a small, low-dimensional subset of modes ($r_\eta \ll q$).  
**The data decisively refute this hypothesis.**
- To capture $50\%$ of the well sacrifice, the system requires $r_{0.50} \approx 0.48 q$.
- To capture $95\%$ of the well sacrifice, the system requires $r_{0.05} \approx 0.90 q$.
- The ratio $r_{0.05}(N) / q(N)$ is extraordinarily stable across dimensions:
  $$0.909, \quad 0.895, \quad 0.913, \quad 0.907.$$

### 4.2 The Extensive Scaling Hypothesis
This stability leads directly to the formulation of a new empirical hypothesis:

> **Hypothesis 142.5 (Extensive Deficit Scaling).**  
> In the continuum Galerkin limit $N \to \infty$, the fraction of modes required to pay a fixed portion $(1 - \eta)$ of the well deficit does not vanish, but approaches an extensive constant:
> $$\lim_{N \to \infty} \frac{r_\eta(N)}{q(N)} = c_\eta \in (0, 1).$$
> For $\eta = 0.05$, the finite-$N$ data indicate $c_{0.05} \approx 0.90$.

This means the well deficit is **genuinely extensive** in the continuum truncation dimension. The collective mechanism is high-dimensional; it cannot be replaced by a low-dimensional Galerkin sub-model.

### 4.3 The Physical Mechanism: Dual-Reservoir Architecture
Rather than a low-dimensional active subspace, the coupled ground state $v(1)$ operates via a **dual-reservoir architecture**:
1. **The Probability Reservoir (Top Cluster $j \le 3$):**  
   The top eigenvalues are exponentially compressed ($\delta\nu_1 \sim 10^{-8}, \delta\nu_2 \sim 10^{-5}, \delta\nu_3 \sim 10^{-4}$). The ground state places nearly half its mass ($49.29\%$) here, allowing it to satisfy boundary constraints and optimize kinetic overlap virtually free of well-potential cost ($1.087 \times 10^{-5}$).
2. **The Energy-Paying Reservoir (Deep Bulk $j \ge 4$):**  
   The well deficit $\Delta W(1) = 0.5507$ is paid broadly across the entire continuum spectrum. Modes $j \ge 20$ contribute two-thirds ($67.58\%$) of the energy sacrifice despite carrying only $22.49\%$ of the probability mass.

---

## 5. Forward Strategic Direction: Cell 143 (Spectral-Measure Scaling)

Because $r_\eta(N) / q(N) \approx 0.90$ scales proportionally to the dimension, further discrete eigenvalue counting (e.g. extending to 30 or 40 modes) is redundant.  

**Cell 143 will investigate the continuum limit of the energy-deficit distribution:**

1. **The Energy-Deficit Spectral Measure:**  
   Define the empirical probability measure $\mu_N = \sum_{j=0}^{q-1} P_W^{(N)}(j) \delta_{\delta\nu_j^{(N)}}$, with $\Delta W_N = \int x \, d\mu_N(x)$.  
   Define the normalized **energy-deficit measure**:
   $$d\lambda_N(x) \equiv \frac{x \, d\mu_N(x)}{\Delta W_N}.$$
2. **Spectral Quantiles:**  
   Evaluate the spectral quantiles $Q_{50}(N), Q_{75}(N), Q_{90}(N), Q_{95}(N)$ on the continuous energy axis $x = \delta\nu$, replacing discrete mode indices $j$.
3. **Scaling Collapse Test:**  
   Test whether $\lambda_N(x)$ collapses toward a universal limiting continuum distribution $\lambda_\infty(x)$ as $N \to \infty$, determining whether the extensive deficit reflects a stationary continuum spectral density or a moving edge.
