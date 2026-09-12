# Cell 114 Analytical Note: Two-State Spectral Reordering Anatomy & Localization Invariants

**Companion Computational Script:** [`cell114.py`](file:///c:/data/github/connes-cvs-/cell114.py) | **Verification Log:** `cell114.out`  
**Status:** Working Research Note (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell113.md`](file:///c:/data/github/connes-cvs-/cell113.md); [`cell113.out`](file:///c:/data/github/connes-cvs-/cell113.out); [`cell112a.md`](file:///c:/data/github/connes-cvs-/cell112a.md); [`cell112a.out`](file:///c:/data/github/connes-cvs-/cell112a.out); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary

Cell 113 tested the hypothesis that the lowest finite-rank Galerkin eigenstate ceases to coincide with the continuum solitary wave for $N \ge 64$. While Cell 113 confirmed the reality of the spectral reordering numerically (eigensolver residuals $\sim 10^{-70}$), its single-vector overlap continuation algorithm:
$$k_{\mathrm{sol}}(N) = \operatorname{argmax}_k \big| \langle v_{N_{\mathrm{prev}}}^{\mathrm{sol}}, v_N^{(k)} \rangle \big|$$
remained locked at $k_{\mathrm{sol}} = 0$ across the entire sweep $N \in [24, 96]$.

### The Methodological Obstruction in Cell 113
The reason for this lock-in is revealed by the detailed overlap progression across the transition region:
- At $N = 44$: $\mathcal{O}_0 = 0.9999, \quad v_0 = 0.5350$
- At $N = 48$: $\mathcal{O}_0 = 0.9844, \quad v_0 = 0.4558$
- At $N = 52$: $\mathcal{O}_0 = 0.8760, \quad v_0 = 0.1629$
- At $N = 56$: $\mathcal{O}_0 = 0.9935, \quad v_0 = 0.0881$
- At $N = 60$: $\mathcal{O}_0 = 0.9998, \quad v_0 = 0.0764$

Because the lowest eigenvector rotates gradually through the two-dimensional eigenspace across $N \in [48, 56]$, its overlap with the preceding state remained larger than the overlap with the excited state ($\mathcal{O}_0 \approx 0.876 > \mathcal{O}_1$). Consequently, the continuation algorithm followed the state into the delocalized edge branch ($v_0 \to 0.064$), and Part C computed extinction on the exact same state as Cell 111, reproducing the $2 \times 10^{-21}$ plateau identically ($P_\alpha^{\mathrm{sol}} / P_\alpha^{(0)} = 1.0$).

### Strategic Mandate for Cell 114
Cell 114 resolves this methodological limitation by introducing:
1. **Two-Dimensional Subspace SVD & Principal Angles:** Tracking the invariant 2D subspace $\mathcal{V}_2 = \operatorname{span}\{v^{(0)}, v^{(1)}\}$ across a fine grid $N \in [40, 60]$ ($\Delta N = 2$), measuring principal singular values $\sigma_1, \sigma_2$ and the continuous rotation angle $\phi(N) = \arctan(|M_{01}| / |M_{00}|)$.
2. **Operator-Independent Physical Localization Invariants:** Characterizing eigenstates by physical observables rather than discrete eigenvalue indices:
   - Core mass concentration: $L_{24}(v) \equiv \sum_{m=0}^{24} v_m^2$.
   - Kinetic Sobolev moment: $\mathcal{K}_2(v) \equiv \sum_{m=1}^N m^4 v_m^2$.
   - Central amplitude: $v_0$.
3. **High-$N$ Branch Extinction Audit:** Evaluating both states ($k=0$ and $k=1$) at $N \in \{64, 80, 96, 128, 192\}$ to measure boundary observables ($T_v(0), \alpha_N, P_\alpha(N)$) on the state that genuinely maintains solitary localization ($L_{24} \approx 1, v_0 \approx 0.666$).
4. **Finite-$T$ Structure at Fixed $N = 48$:** Re-evaluating both states across cached cutoffs $T \in \{400, 500, 600\}$.

---

## 1. Two-Dimensional Subspace Tracking (Module 1)

### 1.1 Mathematical Formulation
Let $V_A = [v_A^{(0)}, v_A^{(1)}] \in \mathbb{R}^{(N_A+1) \times 2}$ and $V_B = [v_B^{(0)}, v_B^{(1)}] \in \mathbb{R}^{(N_B+1) \times 2}$ denote the orthonormal bases for the two lowest eigenspaces at adjacent dimensions $N_A < N_B$. We embed $V_A$ into $\mathbb{R}^{N_B+1}$ by zero-padding higher modes.

The $2 \times 2$ inter-dimensional overlap matrix is:
$$M = V_A^T V_B = \begin{pmatrix} \langle v_A^{(0)}, v_B^{(0)} \rangle & \langle v_A^{(0)}, v_B^{(1)} \rangle \\ \langle v_A^{(1)}, v_B^{(0)} \rangle & \langle v_A^{(1)}, v_B^{(1)} \rangle \end{pmatrix}.$$

The singular values $\sigma_1 \ge \sigma_2 \ge 0$ of $M$ are the square roots of the eigenvalues of $M^T M$.
The principal angles $\theta_1, \theta_2 \in [0, \pi/2]$ between the two subspaces satisfy:
$$\cos \theta_1 = \sigma_1, \qquad \cos \theta_2 = \sigma_2.$$

### 1.2 Diagnostic Hypotheses
- **Subspace Invariance Hypothesis:** If the low-energy physics is governed by an isolated 2-state manifold, then:
  $$\sigma_1 \approx 1, \qquad \sigma_2 \approx 1, \qquad \det(M) \approx 1.$$
- **Planar Rotation Dynamics:** The internal rotation angle between the two basis states across the transition is:
  $$\phi(N) \equiv \arctan\left(\frac{|M_{01}|}{|M_{00}|}\right).$$
  - Before the transition ($N \le 44$): $M_{00} \approx 1, M_{01} \approx 0 \implies \phi \approx 0^\circ$.
  - At the transition midpoint ($N \approx 50-52$): $M_{00} \approx M_{01} \approx 1/\sqrt{2} \implies \phi \approx 45^\circ$ (maximum hybridization).
  - After the transition ($N \ge 58$): $M_{00} \approx 0, M_{01} \approx 1 \implies \phi \approx 90^\circ$ (complete basis exchange $v^{(0)} \leftrightarrow v^{(1)}$).

---

## 2. Operator-Independent Physical Localization Invariants (Module 2)

To avoid circularity or reliance on eigenvalue sorting, Cell 114 defines physical localization metrics:

1. **Core Mass Concentration $L_M(v)$:**
   $$L_M(v) \equiv \sum_{m=0}^M v_m^2 \quad (M \in \{8, 16, 24\}).$$
   - For a localized solitary wave: $L_{24}(v) \approx 0.999+$ (virtually all $L^2$ weight concentrated in low Fourier modes).
   - For a delocalized edge/background mode: $L_{24}(v)$ is significantly lower as modal weight diffuses into higher modes.

2. **Kinetic Sobolev Moment $\mathcal{K}_2(v)$:**
   $$\mathcal{K}_2(v) \equiv \sum_{m=1}^N m^4 v_m^2.$$
   - For the continuum solitary wave: $\mathcal{K}_2 \approx 83.1$ (stabilizes to a finite constant independent of $N$).
   - For a background mode: $\mathcal{K}_2$ scales with higher powers of $N$ due to boundary and high-frequency oscillations.

3. **Branch Identification Criterion:**
   $$k_{\mathrm{loc}}(N) \equiv \operatorname{argmax}_{k \in \{0, 1\}} L_{24}\big(v_N^{(k)}\big).$$
   This provides a robust, operator-independent definition of the solitary wave branch.

---

## 3. High-$N$ Branch Comparison & Extinction Audit (Module 3)

For $N \in \{64, 80, 96, 128, 192\}$, Cell 114 evaluates both states $k=0$ and $k=1$ side-by-side:
- Energy $E_k$
- Central amplitude $v_0(k)$
- Core mass $L_{24}(k)$
- Kinetic moment $\mathcal{K}_2(k)$
- Boundary contact defect $T_{v_k}(0) = v_{k, 0} + \sqrt{2}\sum_{m=1}^N v_{k, m}$
- Boundary coupling scalar $\alpha_N(k) = \sum_{m=1}^N a_m v_{k, m}$
- Extinction product $P_\alpha(k) = |\alpha_N(k)| \sqrt{N}$

### Testable Hypotheses for High $N$:
- **Hypothesis 3A (State Bifurcation):** State $k=1$ maintains $v_0 \approx 0.666$ and $L_{24} \approx 1$, identifying it as the solitary wave continuation, while state $k=0$ maintains $v_0 \approx 0.064$ and $E_0 \approx -10^{-51}$.
- **Hypothesis 3B (Extinction on $k=1$):** Does the boundary defect $P_\alpha(k=1)$ decay toward zero as $N$ increases, or does it also exhibit a plateau?

---

## 4. Finite-$T$ Structure at Fixed $N = 48$ (Module 4)

Using the cached Galerkin matrices at $T \in \{400, 500, 600\}$ (retrieved with zero additional quadrature overhead), Cell 114 evaluates $E_k, v_0, L_{24}, \mathcal{K}_2, T_v(0), \alpha_{48}$ for both states $k \in \{0, 1\}$ to characterize the spectral sensitivity to the finite-$T$ boundary truncation.

---

## 5. Pre-Flight Specification for `cell114.py`

| Parameter / Module | Specification | Purpose |
| :--- | :--- | :--- |
| **Precision** | `mp.mp.dps = 70`, `GROUND_DPS = 70` | Instant cache hit for $N \le 192$ at $T = 600$ |
| **Fine Sweep Grid** | $N \in \{40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60\}$ | Resolves transition anatomy with $\Delta N = 2$ |
| **Subspace SVD** | 2D overlap matrix $M = V_{\mathrm{old}}^T V_{\mathrm{new}}$ | Measures principal angles $\theta_1, \theta_2$ and rotation $\phi$ |
| **Physical Invariants** | $L_8, L_{16}, L_{24}, \mathcal{K}_2, S_2, v_0$ | Operator-independent localized branch identification |
| **High-$N$ Span** | $N \in \{64, 80, 96, 128, 192\}$ | Side-by-side extinction audit on both $k=0$ and $k=1$ |
| **Fixed-$N$ Sweep** | $T \in \{400, 500, 600\}$ at $N = 48$ | Instant cached finite-$T$ sensitivity test |
| **Output Style** | Dispassionate, objective numerical tables | Adheres to repository dry-output standard |
