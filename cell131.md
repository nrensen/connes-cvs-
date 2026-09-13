# CELL 131 — ARCHIMEDEAN OFF-DIAGONAL CONTROL & PROJECTED NEGATIVE-POTENTIAL COMPRESSION ON $\Phi^\perp$

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.5)  
**Target Proposition:** Subspace Coercivity & Operator Domination on the Continuum Subspace:
$$\widehat{Q}_{\mathrm{even}} \equiv U_{\mathrm{cont}}^T Q_{\mathrm{even}} U_{\mathrm{cont}} = \widehat{\Omega} + \widehat{\Delta Q}_{\mathrm{arch}} - \widehat{\mathcal{K}}_{\mathrm{neg}} + \widehat{Q}_{\mathrm{pole}} \succeq E_{11} I \approx 0.58 I > 0$$
on the codimension-11 continuum spectral subspace $\Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}}) \subset \mathcal{H}_{\mathrm{even}}(N)$.  
**Verification / Falsification Criteria:**  
1. **Probe 1 (Archimedean Off-Diagonal Control):** Measure the compressed operator norm $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} = \|U_{\mathrm{cont}}^T [Q_{\mathrm{arch}}^{\mathrm{even}} - \operatorname{diag}(h_+(a_m))] U_{\mathrm{cont}}\|_{\mathrm{op}}$ across $N \in \{16, 20, 24, 28, 32, 40, 48, 64\}$. If $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} \gg c_0 \approx 0.1567$, off-diagonal divided differences dominate the diagonal multiplier and the scalar multiplier cannot serve as the autonomous lower bound.
2. **Probe 2 (Projected Negative Potential Control):** Measure $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) = \lambda_{\max}(U_{\mathrm{cont}}^T (\widetilde{W} - \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}})$. If $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) < \lambda_{\min}(\widehat{\Omega})$, direct operator form domination $\widehat{\Omega} \succ \widehat{\mathcal{K}}_{\mathrm{neg}}$ holds on $\Phi^\perp$. If $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) > \lambda_{\min}(\widehat{\Omega})$, coercivity requires localized phase cancellation within the coupled Friedrichs form rather than separate form domination.  
**Companion Computational Script:** [`cell131.py`](file:///c:/data/github/connes-cvs-/cell131.py)  

---

## 1. Executive Summary & Epistemic Progression

### 1.1 The Context of the Investigation
The progression through Cells 127–130 successfully opened the black box of the prime operator:
1. **Cell 127:** Discovered the exact physical-space shifted autocorrelation representation:
   $$\langle v, Q_{\mathrm{prime}}^{\mathrm{even}} v \rangle = -\frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt.$$
2. **Cell 128:** Proved the Identical Vanishing Theorem $\mathcal{B}_q[v] \equiv 0$ for boundary mass mismatch via midpoint reflection symmetry $T_v(L - t) \equiv T_v(t)$.
3. **Cell 129:** Proved the **Two-Regime Multiplier Theorem**, establishing that the effective diagonal multiplier $\Omega(m) \equiv h_+(a_m) + 4 M(m)$ is unconditionally strictly positive across all modes:
   $$\inf_{m \ge 1} \Omega(m) \ge c_0 \equiv h_+(a_3) \approx +0.1567 > 0,$$
   completely overcoming the Diophantine incommensurability obstruction $\inf_{m \ge 1} M(m) = 0$.
4. **Cell 130:** Proved and certified to 50-digit precision the **Exact Component Decomposition Identity**:
   $$Q_{\mathrm{prime}}^{\mathrm{even}} \equiv -\widetilde{W} + \widetilde{\mathcal{D}}^{\mathrm{per}} + \Delta\widetilde{\mathcal{D}},$$
   with $\widetilde{W} \succeq 0$, $\widetilde{\mathcal{D}}^{\mathrm{per}} \succeq 0$, and $\Delta\widetilde{\mathcal{D}} \preceq 0$.

### 1.2 The Analytical Mission of Cell 131
With the exact tripartite decomposition certified, the central question of Gate 1 has been distilled down to its fundamental core:

$$\textbf{Can the negative operator } \mathcal{K}_{\mathrm{neg}} \equiv \widetilde{W} - \Delta\widetilde{\mathcal{D}} \textbf{ be controlled relative to the positive backbone on } \Phi^\perp?$$

To answer this, Cell 131 investigates two decoupled, precise analytical questions:
1. **Archimedean Off-Diagonal Control:** How large is the off-diagonal divided-difference remainder $\Delta Q_{\mathrm{arch}} = Q_{\mathrm{arch}}^{\mathrm{even}} - \operatorname{diag}(h_+(a_m))$ when projected onto $\Phi^\perp$?
2. **Projected Negative-Potential Control:** Does orthogonal projection onto the continuum spectral subspace $\Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}})$ suppress the negative potential $\mathcal{K}_{\mathrm{neg}} = \widetilde{W} - \Delta\widetilde{\mathcal{D}}$ sufficiently to permit positive form domination?

---

## 2. Mathematical Formulation of the Continuum Compression

### 2.1 The Codimension-11 Spectral Subspace $\Phi^\perp$
On the discrete even Galerkin space $\mathbb{R}^{N+1}$ ($N \ge 16$), the full Friedrichs matrix is:
$$Q_{\mathrm{even}} = Q_{\mathrm{arch}}^{\mathrm{even}} + Q_{\mathrm{prime}}^{\mathrm{even}} + Q_{\mathrm{pole}}^{\mathrm{even}}.$$
Let $Q_{\mathrm{even}} u_k = E_k u_k$ be its complete orthonormal eigendecomposition, sorted in ascending energy order:
$$E_0 \le E_1 \le \dots \le E_{10} < E_{11} \le \dots \le E_N.$$
For $c = 13$, the potential well carries $\bar{N}_{\mathrm{bound}} \approx 11$ bound states ($E_0, \dots, E_{10} \ll 1$), separated from the continuum by the macroscopic boundary gap:
$$g_{11} \equiv E_{11} - E_{10} \approx 0.58 - 10^{-12} \approx 0.58 > 0.$$
We define the bound-state projector:
$$P_{\mathrm{bound}} \equiv \sum_{k=0}^{10} u_k u_k^T \in \mathbb{R}^{(N+1) \times (N+1)},$$
and the continuum spectral projector:
$$P_{\mathrm{cont}} \equiv I - P_{\mathrm{bound}} = \sum_{k=11}^N u_k u_k^T \in \mathbb{R}^{(N+1) \times (N+1)}.$$
The continuum subspace $\Phi^\perp = \operatorname{Ran}(P_{\mathrm{cont}})$ has dimension $q \equiv N - 10$.  
Let $U_{\mathrm{cont}} \in \mathbb{R}^{(N+1) \times q}$ be the orthonormal column matrix:
$$U_{\mathrm{cont}} \equiv \big[ u_{11}, u_{12}, \dots, u_N \big].$$
Because $U_{\mathrm{cont}}^T U_{\mathrm{cont}} = I_q$, $U_{\mathrm{cont}}$ is an exact isometry from $\mathbb{R}^q$ into $\mathbb{R}^{N+1}$.

### 2.2 Compressed Operators on $\Phi^\perp$
For any symmetric operator $A$ on $\mathcal{H}_{\mathrm{even}}(N)$, its compression to $\Phi^\perp$ is the $q \times q$ symmetric matrix:
$$\widehat{A} \equiv U_{\mathrm{cont}}^T A U_{\mathrm{cont}}.$$
By construction, the full Friedrichs operator compresses to a strictly diagonal matrix:
$$\widehat{Q}_{\mathrm{even}} = U_{\mathrm{cont}}^T Q_{\mathrm{even}} U_{\mathrm{cont}} = \operatorname{diag}\big(E_{11}, E_{12}, \dots, E_N\big).$$
Consequently:
$$\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) = E_{11} \approx 0.58 > 0.$$
The continuum spectral floor is an exact, unconditional identity on $\Phi^\perp$. The question is how its constituent physical pieces balance.

---

## 3. The Constituent Operators and the Exact Form Balance

From the Cell 130 exact decomposition, $Q_{\mathrm{even}}$ expands into:
$$Q_{\mathrm{even}} = Q_{\mathrm{arch}}^{\mathrm{even}} + \widetilde{\mathcal{D}}^{\mathrm{per}} - \widetilde{W} + \Delta\widetilde{\mathcal{D}} + Q_{\mathrm{pole}}^{\mathrm{even}}.$$
We define the constituent operators as follows:

### 3.1 The Diagonal Positive Multiplier Backbone ($\Omega_{\mathrm{diag}}$)
The continuous Fourier multiplier symbol induces the diagonal matrix:
$$D_{\mathrm{mult}} \equiv \operatorname{diag}\big(h_+(0), h_+(a_1), \dots, h_+(a_N)\big).$$
Adding the diagonal periodic translation-defect matrix $\widetilde{\mathcal{D}}^{\mathrm{per}} = \operatorname{diag}(0, 4M(1), \dots, 4M(N))$:
$$\Omega_{\mathrm{diag}} \equiv D_{\mathrm{mult}} + \widetilde{\mathcal{D}}^{\mathrm{per}} = \operatorname{diag}\big(h_+(0), \Omega(1), \Omega(2), \dots, \Omega(N)\big),$$
where $\Omega(m) = h_+(a_m) + 4M(m)$.  
- For $m \ge 1$: $\Omega(m) \ge c_0 \approx +0.1567 > 0$ by the Two-Regime Multiplier Theorem.
- For $m = 0$: $\Omega(0) = h_+(0) \approx -5.3722 < 0$.

#### Proposition 3.1 (Zero-Mode Isolation on $\Phi^\perp$)
For any unit vector $w \in \mathbb{R}^q$ ($\|w\|_2 = 1$), the corresponding physical state $v = U_{\mathrm{cont}} w \in \Phi^\perp$ satisfies:
$$\langle w, \widehat{\Omega} w \rangle = h_+(0) |v_0|^2 + \sum_{m=1}^N \Omega(m) v_m^2 \ge h_+(0) |v_0|^2 + c_0 (1 - |v_0|^2).$$
Because the bound states $\{u_0, \dots, u_{10}\}$ capture the spatial average of the potential well, the zero-mode leakage into $\Phi^\perp$:
$$\kappa_0(N) \equiv \|P_{\mathrm{cont}} e_0\|_2^2 = \sum_{k=11}^N |(u_k)_0|^2 = 1 - \sum_{k=0}^{10} |(u_k)_0|^2$$
governs the extent to which the negative $h_+(0)$ penetrates the continuum subspace. If $\kappa_0(N) \ll 1$, $\Phi^\perp$ is insulated from $h_+(0)$.

### 3.2 The Archimedean Off-Diagonal Defect ($\Delta Q_{\mathrm{arch}}$)
The divided-difference matrix $Q_{\mathrm{arch}}^{\mathrm{even}}$ deviates from its diagonal multiplier by:
$$\Delta Q_{\mathrm{arch}} \equiv Q_{\mathrm{arch}}^{\mathrm{even}} - D_{\mathrm{mult}}.$$
Its entries are:
$$(\Delta Q_{\mathrm{arch}})_{mn} = \begin{cases} (Q_{\mathrm{arch}}^{\mathrm{even}})_{mn}, & m \ne n, \\ \psi_{\mathrm{arch}}'(m) - h_+(a_m), & m = n. \end{cases}$$
Its compression onto $\Phi^\perp$ is:
$$\widehat{\Delta Q}_{\mathrm{arch}} \equiv U_{\mathrm{cont}}^T \Delta Q_{\mathrm{arch}} U_{\mathrm{cont}}.$$

### 3.3 The Total Negative Potential ($\mathcal{K}_{\mathrm{neg}}$)
The negative contributions to $Q_{\mathrm{prime}}^{\mathrm{even}}$ combine into:
$$\mathcal{K}_{\mathrm{neg}} \equiv \widetilde{W} - \Delta\widetilde{\mathcal{D}} \succeq 0.$$
Because $\widetilde{W} \succeq 0$ and $\Delta\widetilde{\mathcal{D}} \preceq 0$, $\mathcal{K}_{\mathrm{neg}}$ is strictly **positive semi-definite**. Its compression onto $\Phi^\perp$ is:
$$\widehat{\mathcal{K}}_{\mathrm{neg}} \equiv U_{\mathrm{cont}}^T (\widetilde{W} - \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}} \succeq 0.$$

### 3.4 The Exact Subspace Decomposition Identity
Projecting all components onto $\Phi^\perp$ yields the exact algebraic balance:
$$\boxed{\widehat{Q}_{\mathrm{even}} \equiv \widehat{\Omega} + \widehat{\Delta Q}_{\mathrm{arch}} - \widehat{\mathcal{K}}_{\mathrm{neg}} + \widehat{Q}_{\mathrm{pole}},}$$
where $\widehat{Q}_{\mathrm{pole}} \equiv U_{\mathrm{cont}}^T Q_{\mathrm{pole}}^{\mathrm{even}} U_{\mathrm{cont}}$.  
Taking the trace:
$$\operatorname{Tr}(\widehat{Q}_{\mathrm{even}}) \equiv \operatorname{Tr}(\widehat{\Omega}) + \operatorname{Tr}(\widehat{\Delta Q}_{\mathrm{arch}}) - \operatorname{Tr}(\widehat{\mathcal{K}}_{\mathrm{neg}}) + \operatorname{Tr}(\widehat{Q}_{\mathrm{pole}}).$$

---

## 4. Analytical Hypotheses Under Test

### Hypothesis H-OffDiag (Archimedean Off-Diagonal Decoupling)
The off-diagonal divided-difference coupling decays like $\mathcal{O}(|m - n|^{-1})$ away from the diagonal. On the high-frequency subspace $\Phi^\perp$, destructive phase interference suppresses the off-diagonal operator norm:
$$\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} < \delta_{\mathrm{arch}} < \infty.$$
*Test:* Measure $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}}$ across discrete dimensions $N \in [16, 64]$ and compare with the diagonal floor $c_0 \approx 0.1567$.

### Hypothesis H-Supp (Negative-Potential Subspace Suppression)
On the unprojected Galerkin space, the step potential has maximum eigenvalue $\lambda_{\max}(\widetilde{W}) \approx 9.945$, which is large.  
However, the 11 bound states $u_0, \dots, u_{10}$ localize strongly inside the potential well. Orthogonal projection onto $\Phi^\perp$ forces test functions to be orthogonal to the potential well eigenstates:
$$\int_0^L T_v(t) T_{u_k}(t) \, dt = 0 \qquad (k = 0, \dots, 10).$$
*Hypothesis:* The compressed operator $\widehat{\mathcal{K}}_{\mathrm{neg}}$ satisfies:
$$\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) \ll \lambda_{\max}(\mathcal{K}_{\mathrm{neg}}) \approx 9.945.$$
*Quantitative Metric:* The suppression factor:
$$\mathcal{S}_{\mathrm{supp}}(N) \equiv \frac{\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}})}{\lambda_{\max}(\mathcal{K}_{\mathrm{neg}})} \ll 1.$$

### Hypothesis H-Dom (Coupled Backbone Form Domination)
Define the coupled competition operator:
$$\widehat{\mathcal{Q}}_{\mathrm{comp}} \equiv \widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}.$$
*Hypothesis:* The positive translation lift $4M(m)$ and Archimedean growth $h_+(a_m)$ dominate the projected negative potential:
$$\lambda_{\min}\big(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}\big) > 0.$$

---

## 5. Pre-Flight Computational Protocol (`cell131.py`)

The companion Python script [`cell131.py`](file:///c:/data/github/connes-cvs-/cell131.py) implements this audit under strict numerical standards:

1. **Parameters & Precision:**
   - Precision: `mp.mp.dps = 50`
   - Primary Cutoff: $c = 13$, $L = \log 13 \approx 2.56494935746$
   - Parameter $T$: $T = 600$
   - Dimension Grid: $N \in [16, 20, 24, 28, 32, 40, 48, 64]$
2. **Persistent Matrix Retrieval:**
   - Uses `get_galerkin_matrix` from [`cell.py`](file:///c:/data/github/connes-cvs-/cell.py) with persistent JSON caching to assemble $Q_{\mathrm{even}}$ without redundant quadrature evaluations.
3. **Subspace Isometry Construction:**
   - Computes eigendecomposition $Q_{\mathrm{even}} = V \operatorname{diag}(E) V^T$.
   - Sorts eigenvalues via `symmetric_eigenvalues` to guarantee exact ordering.
   - Extracts $U_{\mathrm{cont}} = V[:, 11:N+1]$.
4. **Audited Quantities per Dimension $N$:**
   - Continuum floor check: $\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) \stackrel{?}{=} E_{11}$.
   - Zero-mode leakage: $\kappa_0(N) = \|P_{\mathrm{cont}} e_0\|^2$.
   - Backbone compression: $\lambda_{\min}(\widehat{\Omega})$, $\lambda_{\max}(\widehat{\Omega})$, $\|\widehat{\Omega}\|_F$.
   - Archimedean off-diagonal defect: $\lambda_{\min}(\widehat{\Delta Q}_{\mathrm{arch}})$, $\lambda_{\max}(\widehat{\Delta Q}_{\mathrm{arch}})$, $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}}$.
   - Negative potential compression: $\lambda_{\min}(\widehat{\mathcal{K}}_{\mathrm{neg}})$, $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}})$, $\|\widehat{\mathcal{K}}_{\mathrm{neg}}\|_{\mathrm{op}}$, and suppression ratio $\mathcal{S}_{\mathrm{supp}}$.
   - Coupled competition spectrum: $\lambda_{\min}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}})$, $\lambda_{\max}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}})$.
   - Pole contribution: $\lambda_{\min}(\widehat{Q}_{\mathrm{pole}})$, $\lambda_{\max}(\widehat{Q}_{\mathrm{pole}})$.
   - Trace identity residual: $|\operatorname{Tr}(\widehat{Q}_{\mathrm{even}}) - [\operatorname{Tr}(\widehat{\Omega}) + \operatorname{Tr}(\widehat{\Delta Q}_{\mathrm{arch}}) - \operatorname{Tr}(\widehat{\mathcal{K}}_{\mathrm{neg}}) + \operatorname{Tr}(\widehat{Q}_{\mathrm{pole}})]|$.

---

## 6. Epistemic Labeling & Standards Adherence

In accordance with [`AGENTS.md`](file:///c:/data/github/connes-cvs-/AGENTS.md):
- **Dispassionate Output:** All script printouts report raw numerical norms, eigenvalues, and ratios without speculative narrative.
- **No Local Python Execution:** [`cell131.py`](file:///c:/data/github/connes-cvs-/cell131.py) is authored as a self-contained script ready for external compute node execution.
- **Epistemic Discipline:** Numerical suppression ratios across $N \in [16, 64]$ will be described strictly as *empirical evidence for potential well shielding*, never as *analytical proof of asymptotic extinction*.

---

## References
- [`cell130.md`](file:///c:/data/github/connes-cvs-/cell130.md) — Exact Component Decomposition of $Q_{\mathrm{prime}}^{\mathrm{even}}$
- [`cell129.md`](file:///c:/data/github/connes-cvs-/cell129.md) — Two-Regime Multiplier Theorem & Monotonicity of $h_+(r)$
- [`cell128.md`](file:///c:/data/github/connes-cvs-/cell128.md) — Identical Vanishing of Boundary Mismatch & Step Potential
- [`cell127.md`](file:///c:/data/github/connes-cvs-/cell127.md) — Shifted Autocorrelation Representation of $Q_{\mathrm{prime}}$
- [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md) — Core Submatrix Spectra & Cauchy Interlacing Bound
- [`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md) — Strategic Roadmap (Gate 1 Pipeline)
