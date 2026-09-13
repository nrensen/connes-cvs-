# Cell 121 Analytical Note: Core-Size ($L$) Invariant Product Mapping & The Bound-State-to-Continuum Transition for Gate 1

**Companion Computational Script:** [`cell121.py`](file:///c:/data/github/connes-cvs-/cell121.py) | **Verification Log:** [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out)  
**Status:** Pre-Flight Design / Active Diagnostic Probe (Gate 1 / Milestone M-G1.0 / Core-Size Scaling)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell120.md`](file:///c:/data/github/connes-cvs-/cell120.md); [`cell120.out`](file:///c:/data/github/connes-cvs-/cell120.out); [`cell96.md`](file:///c:/data/github/connes-cvs-/cell96.md); [`cell97.md`](file:///c:/data/github/connes-cvs-/cell97.md); [`cell60.md`](file:///c:/data/github/connes-cvs-/cell60.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §8.26–8.27  
**Date:** September 2026  

---

## 1. Context & Lessons from Cell 120

Cell 120 provided three decisive insights into the Gate 1 architecture:

### 1.1 Route 1B is Decisively Dead for Fixed Small $L$
The pre-flight conjecture that $R_{\mathrm{spec}}(N, L) = \mathcal{O}(1)$ for fixed $L = 4$ is completely refuted. While the operator norm $\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}$ remains bounded ($3.32 \to 6.47$), the lower Ritz gap collapses exponentially:
$$g_{2, 4}(16) = 6.43 \times 10^{-12} \quad \longrightarrow \quad g_{2, 4}(192) = 2.75 \times 10^{-27}.$$
Consequently, the spectral denominator factor explodes by 31 orders of magnitude:
$$R_{\mathrm{spec}}(16, 4) \approx 8.02 \times 10^{22} \quad \longrightarrow \quad R_{\mathrm{spec}}(192, 4) \approx 8.56 \times 10^{53}.$$
The remote spectrum factor is **not benign** when the tail boundary is placed inside the low-energy well.

### 1.2 Physical Cause: Collective Bound-State Clustering Near Zero
The gap $g_{2, 4}(N) \equiv E_5^{(N)} - E_3^{(N)}$ does not collapse because of a pathological isolated eigenvalue; rather, **the entire low-lying positive spectrum collapses toward zero together**:
$$\text{At } N = 192: \qquad E_0 \sim -10^{-51}, \quad E_2 \approx 6.60 \times 10^{-45}, \quad E_3 \approx 1.47 \times 10^{-38}, \quad E_5 \approx 2.75 \times 10^{-27}.$$
As established in Paper NR2 §8, the Galerkin potential well possesses a bound-state capacity of $\bar{N}_{\mathrm{bound}} \approx 11$ states beneath the barrier top. Setting $L = 4$ forces the tail index $L+1 = 5$ to fall **inside this bound-state cluster**, where level spacing is exponentially suppressed.

### 1.3 The Core-Size Discovery: Tunneling Overwhelms Spectral Crowding for $L \ge 8$
In sharp contrast to $L = 4$ (where the Gate 1 product $\Delta_2 R_{\mathrm{spec}}$ reaches $5 \times 10^{12}$), for $L = 8$ the product is **astonishingly small and rapidly decaying**:
$$\mathcal{P}_2(N, L=8): \quad 3.54 \times 10^{-19} \quad (N=16) \;\longrightarrow\; 3.88 \times 10^{-23} \quad (N=64),$$
even though $R_{\mathrm{spec}}(64, 8) \approx 1.53 \times 10^{18}$ is still large! Tunneling suppression $\Delta_2(N)$ dominates residual spectral crowding by over 20 decades once $L$ is large enough.

### 1.4 Category Error in the Tridiagonal Surrogate (Route 1A)
The Galerkin Hamiltonian $H$ is a **dense matrix** with long-range couplings $H_{mn} = \frac{2(m\psi(m) - n\psi(n))}{m^2 - n^2} \ne 0$ for all $|m - n| > 1$. Truncating to nearest-neighbor hopping $H_{m, m+1}$ and applying a scalar Jacobi Agmon formula is an unproven nearest-neighbor surrogate. Furthermore, $V_{\mathrm{eff}}(m) = H_{mm}$ oscillates between $0.26$ and $2.51$, producing 13 alternating turning points rather than a single classical barrier. A genuine Route 1A proof requires dense-matrix Agmon conjugation ($e^{\theta \Phi} H e^{-\theta \Phi}$), continuum barrier transfer, or solitary quasimodes.

### 1.5 Numerical Precision Floor at 70 dps
At 70 dps, the ground doublet splitting $\Delta_0(N)$ hits the eigensolver precision floor near $10^{-50}$ at $N \ge 40$ ($4.89 \times 10^{-51} \to 2.78 \times 10^{-50}$), producing an unphysical sign flip in the estimated finite-difference slope $\sigma_0$. Reliable numerical rates must be tracked strictly within the clean numerical window ($N \le 40$ for $j=0, 1$; $N \le 48-64$ for $j=2$ where $\Delta_2 \ge 10^{-41} \gg 10^{-70}$).

---

## 2. Mathematical Target & Formulation of Cell 121

Cell 121 addresses the primary question formulated by the analyst:
> **Is the apparent success at $L = 8$ robust as $L$ increases, or is it another finite-$N$ coincidence?**

### 2.1 The Gate 1 Central Proposition
Gate 1 does not require $R_{\mathrm{spec}}(N, L) = \mathcal{O}(1)$ for every fixed small $L$. The governing statement is the **joint limit**:
$$\boxed{\lim_{L \to \infty} \limsup_{N \to \infty} \Delta_j(N) R_{\mathrm{spec}}(N, L) = 0 \quad \Longrightarrow \quad \lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1.}$$

### 2.2 Definitions on the 2D $(N, L)$ Grid
For each dimension $N$ and core size $L$ (with $N \ge L + 2$):
1. **Target Bound States:** Doublets $j \in \{0, 1, 2\}$, with even eigenvalues $E_j^{(N)}, E_{j+1}^{(N)}$ and odd eigenvalue $E_{\mathrm{odd}, j}^{(N)}$.
2. **Parity Doublet Splitting:**
   $$\Delta_j(N) \equiv |E_{\mathrm{odd}, j}^{(N)} - E_{\mathrm{even}, j}^{(N)}|.$$
3. **Lower Ritz Gap:**
   $$g_{j, L}(N) \equiv E_{L+1}^{(N)} - E_{j+1}^{(N)}.$$
4. **Spectral Remote-Tail Factor:**
   $$R_{\mathrm{spec}, j}(N, L) \equiv \frac{\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)})}.$$
5. **Gate 1 Invariant Product:**
   $$\mathcal{P}_j(N, L) \equiv \Delta_j(N) R_{\mathrm{spec}, j}(N, L).$$
6. **The Supremum Envelope:**
   $$\mathcal{S}_j(L) \equiv \sup_{N \in [N_{\min}(L), N_{\max}^{\mathrm{rel}}]} \mathcal{P}_j(N, L).$$

### 2.3 The Core-Size Transition Hypothesis
- **Inside the Well ($L < \bar{N}_{\mathrm{bound}} \approx 11$):**
  $E_{L+1}$ is trapped inside the bound-state cluster. The gap $g_{j, L}(N)$ collapses exponentially with $N$, causing $R_{\mathrm{spec}}$ to grow as $e^{c_L N}$.
- **In the Continuum ($L \ge \bar{N}_{\mathrm{bound}} \approx 11$):**
  $E_{L+1}$ enters the delocalized scattering continuum. By Paper NR2 §8.26 and Cell 90, the boundary gap is macroscopic:
  $$g_{11} \approx 0.42 - 0.57 > 0 \quad \text{uniformly in } N.$$
  Consequently, for $L \ge 12$, the denominator $(E_{L+1} - E_j)(E_{L+1} - E_{j+1}) \ge g_{11}^2 \approx 0.18 > 0$ is **strictly macroscopic**, yielding:
  $$R_{\mathrm{spec}, j}(N, L) = \Theta(1) \quad \text{for } L \ge 12.$$
  In this continuum regime, the Gate 1 product collapses at the **pure tunneling rate**:
  $$\mathcal{P}_j(N, L) \asymp \Delta_j(N) \le C_j e^{-\sigma_j N} \longrightarrow 0.$$

---

## 3. Computational Design of `cell121.py`

### 3.1 Parameter Grid
- **Precision:** `mp.mp.dps = 70` (matching cached Hamiltonian at $c = 13, T = 600, N_{\max} = 192$).
- **Core sizes:** $L \in \{4, 6, 8, 10, 12, 14, 16\}$.
  - $L = 4, 6$: Deep inside the bound-state cluster.
  - $L = 8, 10$: Transition zone near the well barrier top.
  - $L = 12, 14, 16$: Continuum scattering band (above the macroscopic gap $g_{11} \approx 0.42$).
- **Dimensions:** $N \in \{16, 20, 24, 28, 32, 36, 40, 48, 64\}$.
  - Reliable range for $j = 2$: $N \in [16, 64]$ (all $\Delta_2 \ge 2.5 \times 10^{-41} \gg 10^{-70}$).
  - Reliable range for $j = 1$: $N \in [16, 40]$ ($\Delta_1 \ge 9.3 \times 10^{-44}$).
  - Reliable range for $j = 0$: $N \in [16, 32]$ ($\Delta_0 \ge 5.4 \times 10^{-46}$).

### 3.2 Key Quantities to Output
1. **Table 1: Lower Ritz Gaps $g_{2, L}(N)$ Across $(N, L)$:** Tracking the transition from exponential gap collapse ($L = 4$) to macroscopic stabilization ($L \ge 12$).
2. **Table 2: Spectral Factor $R_{\mathrm{spec}, 2}(N, L)$ Across $(N, L)$:** Demonstrating the drop from $10^{53}$ to $\mathcal{O}(1)$.
3. **Table 3: The Gate 1 Product $\mathcal{P}_2(N, L)$ Across $(N, L)$:** Demonstrating suppression across the full grid.
4. **Table 4: Envelope $\mathcal{S}_2(L) = \sup_N \mathcal{P}_2(N, L)$ vs $L$:** Mapping the supremum envelope to confirm that $\mathcal{S}_2(L) \to 0$ monotonically.
5. **Table 5: Doublet Comparison ($j = 0, 1, 2$):** Comparing $\mathcal{P}_j(N, L)$ across the three lowest doublets at representative dimensions ($N = 24, 32$).

---

## 4. Pre-Flight Verification Criteria

| Outcome | Quantitative Criterion | Theoretical Implication |
| :--- | :---: | :--- |
| **Continuum Gap Stabilization** | $g_{2, L}(N) \ge 0.20$ for all tested $N$ when $L \ge 12$ | Confirms that placing $L$ past the bound-state well eliminates small denominators unconditionally. |
| **Spectral Factor Saturation** | $R_{\mathrm{spec}, 2}(N, L) \le 50$ for $L \ge 12$ across all $N$ | Re-establishes Route 1B's bounded $R_{\mathrm{spec}}$ mechanism in the proper continuum regime ($L \ge 12$). |
| **Robust Gate 1 Product Extinction** | $\mathcal{P}_2(N, L) \le 10^{-20}$ for all $N \ge 24$ when $L \ge 8$ | Empirically confirms that tunneling suppression completely dominates residual spectral crowding for $L \ge 8$. |
| **Monotone Envelope Collapse** | $\mathcal{S}_2(L) \equiv \sup_{N \ge 20} \mathcal{P}_2(N, L)$ strictly decreases with $L$ | Provides decisive empirical evidence for $\lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{P}_2(N, L) = 0$. |
