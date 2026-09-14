# CELL 141 — SPECTRAL DEGENERACY OF THE WELL OPERATOR: TOP-SPECTRUM SPLITTINGS, CLUSTER MANIFOLDS, AND GROUND-STATE COUPLING

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.6 / Variational Lower Bound & Coupled Operator Geometry)  
**Target Propositions & Tested Hypotheses:**
1. **Theorem 141.1 (Top-Spectrum Splitting Hierarchy & Scenario Classification):**
   On the continuum constraint subspace $\mathcal{B}_{11}^\perp$ (dimension $q = N - 10$), let $\nu_0(N) > \nu_1(N) \ge \nu_2(N) \ge \dots$ denote the descending eigenvalues of the projected step-potential well operator $W_\perp \equiv U_{\mathrm{cont}}^T \widetilde{W} U_{\mathrm{cont}}$.
   Define the $j$-th spectral splitting:
   $$\delta \nu_j(N) \equiv \nu_0(N) - \nu_j(N), \qquad j = 1, 2, \dots, 9.$$
   Classify the asymptotic structure into three mutually exclusive scenarios:
   - **Scenario A (Isolated Exponential Doublet):**
     $$\lim_{N \to \infty} \delta \nu_1(N) = 0, \qquad \text{while} \quad \liminf_{N \to \infty} \delta \nu_2(N) = c_2 > 0.$$
     The well possesses an isolated near-degenerate doublet $\{y_0, y_1\}$, separated by a macroscopic gap $c_2$ from the rest of the spectrum.
   - **Scenario B (Growing Near-Degenerate Cluster):**
     $$\lim_{N \to \infty} \delta \nu_j(N) = 0 \quad \text{for multiple } j \in \{1, \dots, r-1\}.$$
     The top of the well develops an $r$-dimensional near-degenerate band $\mathcal{Y}_{\mathrm{top}}$.
   - **Scenario C (Numerical / Truncation Artefact):**
     The splitting collapse does not persist under precision increases or exhibits non-monotonic instability.
2. **Theorem 141.2 (Effective Degeneracy Dimension & Manifold Counter):**
   For tolerance thresholds $\epsilon \in \{10^{-1}, 10^{-2}, 10^{-3}, 10^{-4}, 10^{-6}\}$, define the degeneracy counter:
   $$r(\epsilon; N) \equiv 1 + \#\{j \ge 1 : \delta \nu_j(N) < \epsilon\}.$$
   Track the effective dimension $r(\epsilon; N)$ across $N \in [24, 28, 32, 40, 48, 64]$ to determine the scaling of the degenerate subspace $\mathcal{Y}_{\mathrm{top}}(\epsilon) \equiv \operatorname{span}\{y_0, \dots, y_{r-1}\}$.
3. **Theorem 141.3 (Physical Ground-State Subspace Coupling & Harvest Distribution):**
   Let $v(1)$ denote the physical coupled ground state of $H(1) = K_{\mathrm{rest}} - W_\perp$ at $\gamma = 1.00$.
   Define the modal mass distribution across well eigenstates:
   $$P_W(j) \equiv |\langle y_j, v(1) \rangle|^2, \qquad \Pi_{\mathrm{top}}(r) \equiv \sum_{j=0}^{r-1} P_W(j).$$
   The total well deficit decomposes into top-sector and bulk contributions:
   $$\Delta W(1) \equiv \nu_0 - \langle v(1), W_\perp v(1) \rangle = \sum_{j=1}^{r-1} \delta \nu_j P_W(j) + \sum_{j \ge r} \delta \nu_j P_W(j).$$
   Test whether the well sacrifice $\Delta W(1) \approx 0.5507$ is paid inside the near-degenerate manifold $\mathcal{Y}_{\mathrm{top}}$ or comes from deep non-degenerate bulk modes.
4. **Diagnostic 141.4 (Hard Pre-Flight Regression Audit at $N=64$):**
   Exact numerical verification of $\omega_0 = 2.9315260463$, $\nu_0 = 4.2604953336$, and $\mu_0 = -0.4869792210$ to machine precision prior to executing the multi-$N$ sweep.

**Companion Computational Script:** [`cell141.py`](file:///c:/data/github/connes-cvs-/cell141.py)  
**Execution Standard:** Self-contained 50-dps verification suite ready for external compute node execution.

---

## 1. Executive Context & Structural Motivation

### 1.1 The Breakthrough and Diagnosis from Cell 140
In Cell 140, an attempt to normalize the Pareto frontier coordinates by the spectral gaps $(u, v) = (\Delta W / \Delta \nu, \Delta K / \Delta \omega)$ and test an analytical 2-level model in $\operatorname{span}\{x_0, x_1\}$ revealed two profound facts:
1. **The First Excited Well Gap Collapses Exponentially:**
   $$\Delta \nu(N) = \nu_0(N) - \nu_1(N): \quad 0.490 \to 0.168 \to 0.0463 \to 0.00229 \to 3.1 \times 10^{-5} \to < 10^{-7}.$$
   This sudden collapse caused the naive 2-level model $W_2 = \nu_1 I_2 + \Delta \nu y_0 y_0^T$ to collapse to a multiple of the identity ($\nu_0 I_2$), trivially producing zero deficits $\Delta K_{\mathrm{2lvl}} = 0, \Delta W_{\mathrm{2lvl}} = 0$.
2. **The Local Pareto Geometry Decouples from $\Delta \nu$:**
   While $\Delta \nu(N) \to 0$ over six orders of magnitude, the physical ground-state Pareto curvature stabilized cleanly:
   $$E''(1): \quad -0.4267 \to -0.4045 \to -0.4007 \qquad (\kappa(1) \to 2.4958 \approx 2.50).$$
   Furthermore, the well-side spectral moments rapidly stabilized ($M_W^{(1)} \approx 1.236, \sigma_W \approx 0.797$), even while the kinetic moments broadened.

### 1.2 The Central Mission of Cell 141
The collapse of $\nu_0 - \nu_1$ is not a failure of the project; it is the discovery of an unexpected **spectral degeneracy in the potential well operator $W_\perp$**.
Cell 141 investigates:
- Is this degeneracy confined to a single pair of states (Scenario A: isolated exponential doublet), or does a macroscopic cluster of eigenvalues coalesce at the top of the well (Scenario B)?
- What is the effective dimension $r(\epsilon; N)$ of the top manifold $\mathcal{Y}_{\mathrm{top}}$?
- Does the physical coupled ground state $v(1)$ live largely inside $\mathcal{Y}_{\mathrm{top}}$, and how much of the well depth does $\mathcal{Y}_{\mathrm{top}}$ capture?

---

## 2. Theoretical Formulation & Hypotheses

### 2.1 The Top 10 Eigenvalue Splittings
Let $W_\perp = U_{\mathrm{cont}}^T \widetilde{W} U_{\mathrm{cont}}$ denote the projected step-potential matrix on $\mathcal{B}_{11}^\perp$ (dimension $q = N - 10$).
Let its sorted eigenvalues be $\nu_0 > \nu_1 \ge \nu_2 \ge \dots \ge \nu_{q-1} > 0$, with corresponding orthonormal eigenvectors $y_0, y_1, y_2, \dots, y_{q-1}$.
Define the splittings relative to the ground state:
$$\delta \nu_j(N) \equiv \nu_0(N) - \nu_j(N), \qquad j = 1, 2, \dots, 9.$$

### 2.2 Classification of the Three Scenarios

#### Scenario A: Isolated Exponential Doublet
If the step potential $W(t)$ and the boundary constraint $\mathcal{B}_{11}$ induce an effective reflection symmetry in the highest bound states of the well, the top two states will form an isolated doublet:
$$\delta \nu_1(N) \sim C_1 e^{-\alpha N} \longrightarrow 0, \qquad \delta \nu_2(N) \longrightarrow c_2 > 0.$$
- **Consequence:** The top sector is strictly two-dimensional: $\mathcal{Y}_{\mathrm{top}} = \operatorname{span}\{y_0, y_1\}$.
- **Significance for Reduced Modeling:** The correct minimal model is a **3-state Hamiltonian** $\mathcal{H}_3 = \operatorname{span}\{x_0, y_0, y_1\}$ rather than the collapsed 2-state model. Because $y_0$ and $y_1$ are degenerate, the coupled state can rotate freely in $\operatorname{span}\{y_0, y_1\}$ without paying potential energy.

#### Scenario B: Growing Near-Degenerate Cluster
If multiple states coalesce at $\nu_0$:
$$\delta \nu_j(N) \longrightarrow 0 \quad \text{for } j = 1, 2, \dots, r-1,$$
then the top of $W_\perp$ forms a macroscopic quasi-continuum or degenerate band.
- **Consequence:** The well admits an $r$-dimensional flat plateau in energy space, explaining why the well-side moments $M_W^{(1)}$ and $\sigma_W$ stabilize rapidly.
- **Significance for Reduced Modeling:** The effective well sector must be modeled by the full projector $P_{\mathrm{top}} = \sum_{j=0}^{r-1} y_j y_j^T$.

#### Scenario C: Truncation Artefact
If the splittings do not scale smoothly with $N$ or the apparent collapse is reversed upon altering the Galerkin cutoff or quadrature precision, the phenomenon is a numerical artefact.

---

### 2.3 Ground-State Projection into $\mathcal{Y}_{\mathrm{top}}$
Let $v(1)$ be the coupled ground state:
$$(K_{\mathrm{rest}} - W_\perp) v(1) = \mu_0 v(1), \qquad \|v(1)\|_2 = 1.$$
Expanding $v(1)$ in the eigenbasis $\{y_j\}$ of $W_\perp$:
$$v(1) = \sum_{j=0}^{q-1} c_j y_j, \qquad c_j = \langle y_j, v(1) \rangle, \qquad P_W(j) = |c_j|^2.$$
The cumulative mass captured by the first $r$ well eigenstates is:
$$\Pi_{\mathrm{top}}(r) \equiv \sum_{j=0}^{r-1} P_W(j).$$
The total well harvest sacrifice $\Delta W(1) \equiv \nu_0 - \langle v(1), W_\perp v(1) \rangle$ satisfies:
$$\Delta W(1) = \sum_{j=1}^{q-1} \delta \nu_j P_W(j) = \underbrace{\sum_{j=1}^{r-1} \delta \nu_j P_W(j)}_{\text{Top-Manifold Deficit}} + \underbrace{\sum_{j=r}^{q-1} \delta \nu_j P_W(j)}_{\text{Bulk Deficit}}.$$
- If $\delta \nu_j \approx 0$ for $j < r$, the top manifold pays **negligible potential deficit**: $\sum_{j=1}^{r-1} \delta \nu_j P_W(j) \approx 0$.
- Therefore, the observed deficit $\Delta W(1) = 0.550712$ must be paid entirely by the projection of $v(1)$ onto non-degenerate bulk modes $j \ge r$.
- Computing $P_W(j)$ tests quantitatively whether $v(1)$ concentrates its mass inside $\mathcal{Y}_{\mathrm{top}}$ or is distributed deep into the well bulk.

---

## 3. Experimental Suite for Cell 141

### 3.1 Dimensions & Parameters
- **Cutoffs:** $c = 13, L = \log 13 \approx 2.56494935746, T = 600, \text{dps} = 50$.
- **Constraint Subspace:** $N_{\mathrm{bound}} = 11$, continuum dimension $q = N - 10$.
- **Dimension Sweep:** $N \in [24, 28, 32, 40, 48, 64]$ ($q = 14, 18, 22, 30, 38, 54$).

### 3.2 Output Tables
1. **Table 1 (Top-10 Eigenvalue Splittings):** Reports $\delta \nu_j(N) = \nu_0(N) - \nu_j(N)$ for $j = 1, \dots, 9$ across all 6 dimensions.
2. **Table 2 (Degeneracy Counter):** Reports $r(\epsilon; N) \equiv 1 + \#\{j \ge 1 : \delta \nu_j(N) < \epsilon\}$ for $\epsilon \in \{10^{-1}, 10^{-2}, 10^{-3}, 10^{-4}, 10^{-6}\}$.
3. **Table 3 (Coupled Ground State Well Distribution):** Reports $P_W(j) = |\langle y_j, v(1) \rangle|^2$ for $j = 0, \dots, 9$ and cumulative mass $\Pi_{\mathrm{top}}(r)$ at $N = 64$.
4. **Table 4 (Scenario Diagnostic & Ratio Tests):** Reports consecutive gap ratios $\delta \nu_2 / \delta \nu_1$, $\delta \nu_3 / \delta \nu_2$, and local exponential decay exponents to classify Scenario A vs B vs C.
