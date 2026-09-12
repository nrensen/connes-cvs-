# Cell 117 Analytical Note: Boundary-Row Modal Profiling & Boundary Layer Scaling Collapse Diagnostic

**Companion Computational Script:** [`cell117.py`](file:///c:/data/github/connes-cvs-/cell117.py) | **Verification Log:** `cell117.out` (pending compute node execution)  
**Status:** Pre-Flight Authoring (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell116.md`](file:///c:/data/github/connes-cvs-/cell116.md); [`cell116.out`](file:///c:/data/github/connes-cvs-/cell116.out); [`cell115.md`](file:///c:/data/github/connes-cvs-/cell115.md); [`cell112a.md`](file:///c:/data/github/connes-cvs-/cell112a.md); [`cell111.md`](file:///c:/data/github/connes-cvs-/cell111.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Mathematical Architecture

Following the analytical reconnaissance in Cell 116, which identified an empirical proportionality $|\alpha_N| / |T_v(0)| \approx 284$ and showed compatibility with an empirical $N^{-1/3}$ scaling envelope, Cell 117 provides the rigorous diagnostic required to test whether the $N^{-1/3}$ behavior is driven by a true discrete boundary-layer mechanism (the Airy hypothesis) or is merely a heuristic coincidence.

### Core Mathematical Objectives

1. **Exact Analytical & Numerical Audit of $H_{0, N}$ and $a_N$:**
   Derive the closed-form identities:
   $$H_{0, N} = \sqrt{2} \frac{\psi(N)}{N}, \qquad a_N = 2 N \psi(N).$$
   Decompose $\psi(N) = \psi_{\mathrm{prime}}(N) + \psi_{\mathrm{pole}}(N) + \psi_{\mathrm{arch}}(N)$ to explain analytically why $a_N / \sqrt{2}$ oscillates wildly in sign (from $+24.17$ to $-145.38$) and grows linearly in envelope as $\mathcal{O}(N)$ due to the almost-periodic, non-decaying prime sum $\psi_{\mathrm{prime}}(N)$.
2. **Modal Boundary-Row Profiling ($j = N - k$):**
   Evaluate the individual modal contributions to the boundary flux:
   $$F_{N, k} \equiv H_{Nk} k^2 v_{N, k}$$
   re-indexed by distance from the boundary $j \equiv N - k \in \{0, 1, 2, \dots, 24\}$ across five dimensions $N \in \{64, 96, 128, 160, 192\}$.
3. **Boundary Layer Scaling Collapse Diagnostic:**
   Directly test whether the modal profile obeys a universal similarity collapse:
   $$F_{N, N-j} \stackrel{?}{\approx} N^{-\mu} f\left(\frac{j}{N^\theta}\right)$$
   for candidate exponents:
   - $\theta = 0$: Discrete fixed-depth boundary layer (independent of $N$).
   - $\theta = 1/3$: Airy Dirichlet boundary layer in modal space ($\Delta j \sim N^{1/3}$, corresponding to $\Delta t \sim N^{-2/3}$ in physical coordinate space).
   - $\theta = 1/2$: Diffusive scaling ($\Delta j \sim \sqrt{N}$).
   - $\theta = 1$: Macroscopic / bulk scaling ($\Delta j \sim N$).
4. **Binned Flux Decomposition & Multi-Order Cancellation Tracking:**
   Quantify the exact destructive cancellation between the core modes ($k \le 24$) and intermediate/boundary modes by computing the peak partial sum $M_{\mathrm{peak}}(N) = \max_M |S_{\mathrm{core}}(M)|$ and the cancellation factor:
   $$\mathcal{C}_{\mathrm{cancel}}(N) \equiv \frac{\max_M |\sum_{k=1}^M F_{N, k}|}{|(H u_N)_N|}.$$

---

## 1. Analytical Dissection of $H_{0, N}$ and the Oscillation of $a_N$

### 1.1 The Exact Operator Identity

In the $(2N+1)$-dimensional exponential basis $e_k(t) = \exp(2\pi i k t / L)$, the Connes–CvS Galerkin matrix $Q$ has entries given by divided differences of the symbol $\psi$:
$$Q_{m, n} = \begin{cases} \dfrac{\psi(m) - \psi(n)}{m - n}, & m \ne n, \\[1ex] \psi'(n), & m = n. \end{cases}$$
The symbol $\psi$ satisfies the exact odd parity identity $\psi(-n) = -\psi(n)$, with $\psi(0) = 0$.

Evaluating row $m = 0$ against column $n = k > 0$:
$$Q_{0, k} = \frac{\psi(0) - \psi(k)}{0 - k} = \frac{-\psi(k)}{-k} = \frac{\psi(k)}{k}.$$

In the canonical $(N+1)$-dimensional $v$-basis, the first row is given by parity projection (cf. `cell.py`, `extract_canonical_H`):
$$H_{0, k} = \sqrt{2} Q_{0, k} = \sqrt{2} \frac{\psi(k)}{k}.$$

Recall that the commutator vector $a_k$ (satisfying $[K^2, H]_{jk} = a_j - a_k$) is defined by:
$$a_k \equiv \sqrt{2} k^2 H_{0, k} = \sqrt{2} k^2 \left( \sqrt{2} \frac{\psi(k)}{k} \right) = 2 k \psi(k).$$
Therefore:
$$\boxed{a_N = 2 N \psi(N)} \qquad \Longleftrightarrow \qquad \boxed{\frac{a_N}{\sqrt{2}} = \sqrt{2} N \psi(N)}.$$

### 1.2 The Structure of the Symbol $\psi(N)$

The full Weil symbol decomposes into three distinct arithmetic contributions:
$$\psi(N) = \psi_{\mathrm{prime}}(N) + \psi_{\mathrm{pole}}(N) + \psi_{\mathrm{arch}}(N).$$

1. **The Prime Piece $\psi_{\mathrm{prime}}(N)$:**
   $$\psi_{\mathrm{prime}}(N) = -\frac{1}{\pi} \sum_{p^k \le c} \frac{\Lambda(p^k)}{\sqrt{p^k}} \sin\left(2\pi N \left(1 - \frac{\log(p^k)}{L}\right)\right).$$
   Using $\sin(2\pi N - \theta) = -\sin\theta$ for integer $N$:
   $$\psi_{\mathrm{prime}}(N) = \frac{1}{\pi} \sum_{p^k \le c} \frac{\Lambda(p^k)}{\sqrt{p^k}} \sin\left(2\pi N \frac{\log(p^k)}{L}\right).$$
   Because the frequencies $\omega_{p^k} \equiv 2\pi \frac{\log(p^k)}{L}$ are incommensurate (as logarithms of distinct prime powers are linearly independent over $\mathbb{Q}$), $\psi_{\mathrm{prime}}(N)$ is a **quasi-periodic, almost-periodic trigonometric polynomial** in $N$.
   - It does **NOT** decay as $N \to \infty$.
   - Its envelope is strictly $\mathcal{O}(1)$.
   - It oscillates indefinitely between positive and negative values.

2. **The Pole Piece $\psi_{\mathrm{pole}}(N)$:**
   $$\psi_{\mathrm{pole}}(N) = \frac{1}{\pi} \int_0^L \sin\left(2\pi N \left(1 - \frac{y}{L}\right)\right) 2 \cosh(y/2) \, dy.$$
   Integrating by parts against the smooth hyperbolic cosine:
   $$\psi_{\mathrm{pole}}(N) = \mathcal{O}\left(\frac{1}{N}\right) \quad \text{as } N \to \infty.$$
   The pole contribution vanishes asymptotically.

3. **The Archimedean Piece $\psi_{\mathrm{arch}}(N)$:**
   $$\psi_{\mathrm{arch}}(N) = \frac{1}{2\pi^2} \int_{-T}^T h_+(\tau) \operatorname{Re}\big(\hat{S}_N(\tau)\big) \, d\tau.$$
   For fixed $T$, this represents a smooth, slowly varying Mellin multiplier baseline.

### 1.3 Resolution of the $a_N$ Oscillation Puzzle

Because $\psi_{\mathrm{prime}}(N) = \mathcal{O}(1)$ oscillates indefinitely without decaying:
$$\psi(N) \sim \mathcal{O}(1) \times \text{oscillatory factor}.$$
Multiplying by $2N$:
$$a_N = 2 N \psi(N) \sim \mathcal{O}(N) \times \text{oscillatory factor}.$$

This completely resolves the observation in `cell116.out`:
- For $N \in \{64, 80, 96, 112, 128, 144, 160, 192\}$, $a_N / \sqrt{2}$ took values:
  $$+5.27, -4.60, -17.45, -13.04, -44.68, -43.07, +24.17, -145.38.$$
- The ratio $\frac{a_N / \sqrt{2}}{N} = \sqrt{2} \psi(N)$ fluctuates between $-0.76$ and $+0.15$, perfectly matching the bounded oscillations of $\psi_{\mathrm{prime}}(N)$.
- **Consequence for Equipartition:** The contact term in Theorem 1 is $\frac{a_N}{\sqrt{2}} T_v(0) = \sqrt{2} N \psi(N) T_v(0)$. Because $a_N$ oscillates with $N$, its ratio to $\alpha_N = \sum_{k=1}^N a_k v_k$ is not a universal constant. The numerical equipartition $(H u_N)_N \approx \frac{1}{2}\alpha_N \approx \frac{a_N}{\sqrt{2}} T_v(0)$ observed at $N = 192$ occurred because $a_{192}/\sqrt{2} \approx -145.38 \approx -284/2$; at $N = 64$, $a_{64}/\sqrt{2} = +5.27$ contributes less than $2\%$ of $\alpha_{64}$. The $1/2$ ratio at $N = 192$ was an accidental feature of $N = 192$, not an asymptotic law.

---

## 2. Formulation of the Boundary Layer Scaling Collapse Diagnostic

### 2.1 The Modal Flux Representation

The exact boundary flux of the localized solitary wave is:
$$(H u_N)_N = \sum_{k=1}^N H_{Nk} k^2 v_{N, k} = \sum_{j=0}^{N-1} F_{N, N-j},$$
where $j \equiv N - k$ measures the distance from the boundary cutoff:
$$F_{N, N-j} \equiv H_{N, N-j} (N - j)^2 v_{N, N-j}.$$
- $j = 0$ is the diagonal endpoint term: $F_{N, N} = H_{NN} N^2 v_{N, N}$.
- $j = 1$ is the first off-diagonal boundary mode: $F_{N, N-1} = H_{N, N-1} (N-1)^2 v_{N, N-1}$.
- $j \ge 1$ modes represent the high-frequency edge of the Galerkin domain.

### 2.2 Candidate Scaling Hypotheses

The general scaling collapse ansatz posits:
$$F_{N, N-j} \approx N^{-\mu} f\left(\frac{j}{N^\theta}\right).$$
In Cell 116, the total flux $(H u_N)_N$ scaled roughly as $N^{-1/3}$, suggesting $\mu = 1/3$.

Cell 117 tests four distinct physical regimes for the scaling exponent $\theta$:

1. **Hypothesis $\theta = 0$ (Discrete Fixed-Depth Boundary Layer):**
   The boundary layer is purely lattice-discrete and independent of $N$:
   $$F_{N, N-j} \approx N^{-1/3} f_0(j).$$
   *Prediction:* At any fixed discrete distance $j \in \{0, 1, 2, \dots\}$, the quantity $N^{1/3} F_{N, N-j}$ is constant across $N$.
2. **Hypothesis $\theta = 1/3$ (Airy Boundary Layer):**
   The boundary layer width in modal frequency space scales as $\Delta j \sim N^{1/3}$:
   $$F_{N, N-j} \approx N^{-1/3} f_{\mathrm{Airy}}\left(\frac{j}{N^{1/3}}\right).$$
   *Prediction:* Points with matched dimensionless coordinate $\eta \equiv j / N^{1/3}$ collapse onto a universal profile $f_{\mathrm{Airy}}(\eta)$.
3. **Hypothesis $\theta = 1/2$ (Diffusive Boundary Layer):**
   The boundary layer width scales as $\Delta j \sim \sqrt{N}$.
4. **Hypothesis $\theta = 1$ (Macroscopic Bulk Scaling):**
   The boundary terms reflect global domain rescaling $\eta = j / N$, with no localized boundary layer.

---

## 3. Multi-Order Destructive Cancellation Tracking

The total boundary flux $(H u_N)_N = -8.35 \times 10^{-23}$ at $N = 192$ is ten orders of magnitude smaller than the partial sum over the low-mode core ($k \le 24$: $-5.03 \times 10^{-13}$).

To quantify whether this cancellation is a systematic feature of the operator or an isolated anomaly, Cell 117 evaluates:
1. **The Peak Partial Sum:**
   $$M_{\mathrm{peak}}(N) \equiv \max_{1 \le M \le N} \left| \sum_{k=1}^M F_{N, k} \right|.$$
2. **The Cancellation Factor:**
   $$\mathcal{C}_{\mathrm{cancel}}(N) \equiv \frac{M_{\mathrm{peak}}(N)}{|(H u_N)_N|}.$$

### Decision Tree for Cancellation Scaling

- **If $\mathcal{C}_{\mathrm{cancel}}(N) \sim \mathcal{O}(1)$:** The cancellation is mild; boundary-row terms directly govern the boundary flux without significant destructive interference.
- **If $\mathcal{C}_{\mathrm{cancel}}(N) \sim N^p$ (Power-Law Growth):** Systematic algebraic destructive phase cancellation quenches the core flux into the boundary layer.
- **If $\mathcal{C}_{\mathrm{cancel}}(N) \sim e^{c N}$ (Exponential Growth):** The boundary flux is governed by an exponential mode-by-mode cancellation mechanism, explaining why macroscopic sector estimates fail to capture the net scaling.

---

## 4. Falsification Criteria for Cell 117

| Quantitative Metric | Outcome A (Airy Hypothesis Confirmed) | Outcome B (Discrete $\theta=0$ Boundary Layer) | Outcome C (Airy Hypothesis Refuted / Phase Quenching) |
| :--- | :--- | :--- | :--- |
| **Scaling Collapse Test** | $N^{1/3} F_{N, N-j}$ collapses onto universal curve vs $\eta = j / N^{1/3}$ (spread $< 15\%$). | $N^{1/3} F_{N, N-j}$ is stationary at fixed discrete $j$ (spread $< 15\%$). | Neither collapses; modal terms oscillate in sign or fluctuate with prime frequencies (spread $> 50\%$). |
| **$a_N$ Decomposition** | Confirmed: $a_N = 2 N \psi(N)$, oscillations driven by $\psi_{\mathrm{prime}}(N)$. | Confirmed: $a_N = 2 N \psi(N)$, oscillations driven by $\psi_{\mathrm{prime}}(N)$. | Confirmed: $a_N = 2 N \psi(N)$, oscillations driven by $\psi_{\mathrm{prime}}(N)$. |
| **Cancellation Ratio $\mathcal{C}_{\mathrm{cancel}}$** | Moderately growing or stable. | Moderately growing or stable. | Rapidly rising across $N \in [64, 192]$, confirming phase cancellation dominance. |
| **Architectural Verdict** | Airy boundary layer confirmed as genuine physical mechanism. Promote to NR2 with rigorous model. | Discrete boundary layer confirmed; continuum limit is non-Airy. | Airy hypothesis abandoned; pivot to phase-cancellation and coordinate transport models. |

---

## 5. Summary of Pre-Flight Script Configuration (`cell117.py`)

- **Dimensions Tested:** $N \in \{64, 96, 128, 160, 192\}$.
- **Parameters:** $c = 13$, $T = 600$, precision $\mathrm{dps} = 70$.
- **Sub-Matrix Extraction:** Re-uses cached $N = 192$ Hamiltonian; total runtime estimated at $\approx 300 - 450\text{ s}$ on compute node.
- **Execution Sentinel:** Clean 3-line termination sentinel present.
