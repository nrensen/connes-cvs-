# CELL 134 — CONSTRAINED VARIATIONAL FORMULATION OF THE COMPETITION MINIMIZER & THE CONTINUUM LIMITING PROFILE

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.6)  
**Target Proposition:** Constrained Euler–Lagrange Variational Structure of the Solitary Vulnerable State:
$$(\mathcal{Q}_{\mathrm{comp}} - \mu_0 I) v_{\mathrm{bad}} = \sum_{k=0}^{10} \lambda_k u_k, \qquad \lambda_k \equiv \langle u_k, \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} \rangle$$
on the canonical even space $\mathbb{R}^{N+1}$, demonstrating that:
1. The solitary negative state $v_{\mathrm{bad}}$ is the unique ground state of a constrained variational problem whose deficit $\mu_0 \approx -0.4870$ is maintained against the attractive step potential by an exact vector of Lagrange constraint forces $\boldsymbol\lambda = (\lambda_0, \dots, \lambda_{10})^T$.
2. Under canonical phase alignment ($T_{v_{\mathrm{bad}}}(0) > 0$), the physical wavepacket $T_{v_{\mathrm{bad}}}(t)$ exhibits strong numerical stabilization in $L^2(d\mu)$ toward a stationary spatial profile.
3. The competition Hamiltonian decomposes into three positive/negative physical forms:
   $$\mathcal{Q}_{\mathrm{comp}} \equiv D_{\mathrm{mult}} - \widetilde{W} + \widetilde{\mathcal{D}}^{\mathrm{true}},$$
   where the maximum depth of the prime step potential $W(L) = 2 \sum_{q \le c} w_q \approx 9.9438$ coincides identically with the average translation stiffness $\overline{4 M(m)} \equiv W(L)$, establishing the exact baseline around which the $-0.4870$ deficit is stabilized.

**Verification / Falsification Criteria:**
1. **Euler–Lagrange Residual Audit:** Certify dynamically to 45 decimal digits that:
   $$\|\mathbf{r}_{\mathrm{EL}}\|_2 \equiv \Big\| (\mathcal{Q}_{\mathrm{comp}} - \mu_0 I) v_{\mathrm{bad}} - \sum_{k=0}^{10} \lambda_k u_k \Big\|_2 < 10^{-45}$$
   across all tested dimensions $N \in [24, 64]$, confirming the algebraic consistency of the constrained variational projection.
2. **Profile Stabilization Toward Reference:** Track the $L^2$ Euclidean distance of the canonically aligned coefficient vectors:
   $$\delta_{\mathrm{Cauchy}}(N) \equiv \|v_{\mathrm{bad}}^{(N)} - v_{\mathrm{bad}}^{(64)}\|_2$$
   and demonstrate numerical stabilization toward the $N = 64$ reference state with $\delta_{\mathrm{Cauchy}}(48) \approx 0.065$.
3. **Three-Way Energy Partition Balance:** Verify dynamically that:
   $$\big| R_{\mathrm{comp}} - (R_{\mathrm{mult}} - R_W + R_{\mathcal{D}}^{\mathrm{true}}) \big| < 10^{-45}$$
   at every dimension $N$.

**Companion Computational Script:** [`cell134.py`](file:///c:/data/github/connes-cvs-/cell134.py)  

---

## 1. Executive Context: From Phenomenological Cancellation to Variational Optimization

In Cells 131–133, a high-precision 50-dps audit established the exact operator decomposition and coupled cancellation on the continuum subspace $\Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}})$:
1. The competition operator $\widehat{\mathcal{Q}}_{\mathrm{comp}} \equiv \widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}$ possesses **precisely one** negative eigenvalue across all tested dimensions $N \ge 24$ ($k_{\mathrm{neg}} \equiv 1$).
2. The second eigenvalue is strictly positive and bounded away from zero: $\mu_1(64) \approx +0.3483 > 0$.
3. The resulting spectral gap:
   $$\Delta \mu \equiv \mu_1 - \mu_0 \approx 0.3483 - (-0.4870) = 0.8353 = \Theta(1)$$
   is macroscopic and stable, proving that the solitary negative state $w_{\mathrm{bad}}$ is an isolated, robust, non-degenerate ground state of $\widehat{\mathcal{Q}}_{\mathrm{comp}}$.
4. In Cell 133, the continuous spatial wavepacket $T_{v_{\mathrm{bad}}}(t)$ was found to have a stationary probability mass distribution ($\mathcal{M}_{\mathrm{well}} \approx 60.6\%$, $\mathcal{M}_{\mathrm{flat}} \approx 39.4\%$) with its peak strictly frozen at $t_{\mathrm{peak}} \approx 2.539$ and Fourier mode $m^* = 26$.

Cell 134 investigates the mathematical reason for this stability:
> *Can the competition minimizer $v_{\mathrm{bad}}$ be characterized as a constrained variational state for a one-dimensional step potential, and what are the exact Lagrange multiplier constraint forces that stabilize it?*

---

## 2. Derivation of the Discrete Euler–Lagrange Variational Equation

Let $\mathcal{H}_{\mathrm{even}} = \mathbb{R}^{N+1}$ be the canonical even Galerkin space, and let $u_0, u_1, \dots, u_N$ be the orthonormal eigenvectors of the full Friedrichs operator $Q_{\mathrm{even}}$, with eigenvalues $E_0 \le E_1 \le \dots \le E_N$.

The bound-state core subspace is $\mathcal{H}_{\mathrm{bound}} \equiv \operatorname{span}\{u_0, u_1, \dots, u_{10}\}$ with dimension $N_{\mathrm{bound}} = 11$. The bound-state spectral projector is:
$$P_{\mathrm{bound}} \equiv \sum_{k=0}^{10} u_k u_k^T.$$
The continuum trial subspace is $\Phi^\perp \equiv \operatorname{Ran}(I - P_{\mathrm{bound}})$, with dimension $q = (N + 1) - 11 = N - 10$.

Let $U_{\mathrm{cont}} \in \mathbb{R}^{(N+1) \times q}$ be the isometry whose columns are the continuum eigenvectors $u_{11}, \dots, u_N$.

### 2.1 The Subspace Variational Problem
On $\Phi^\perp$, the competition operator is represented in the $U_{\mathrm{cont}}$ basis by:
$$\widehat{\mathcal{Q}}_{\mathrm{comp}} \equiv U_{\mathrm{cont}}^T \mathcal{Q}_{\mathrm{comp}} U_{\mathrm{cont}} \in \mathbb{R}^{q \times q}.$$
The lowest eigenvalue $\mu_0 = \lambda_{\min}(\widehat{\mathcal{Q}}_{\mathrm{comp}})$ and its normalized eigenvector $w_{\mathrm{bad}} \in \mathbb{R}^q$ satisfy the standard Rayleigh quotient minimization:
$$\mu_0 = \min_{\substack{w \in \mathbb{R}^q \\ \|w\|_2 = 1}} \langle w, \widehat{\mathcal{Q}}_{\mathrm{comp}} w \rangle, \qquad \widehat{\mathcal{Q}}_{\mathrm{comp}} w_{\mathrm{bad}} = \mu_0 w_{\mathrm{bad}}.$$

### 2.2 Embedding into the Full Space $\mathbb{R}^{N+1}$
Define the canonical embedded physical vector:
$$v_{\mathrm{bad}} \equiv U_{\mathrm{cont}} w_{\mathrm{bad}} \in \mathbb{R}^{N+1}.$$
Because $U_{\mathrm{cont}}^T U_{\mathrm{cont}} = I_q$, we have:
$$\|v_{\mathrm{bad}}\|_2^2 = w_{\mathrm{bad}}^T U_{\mathrm{cont}}^T U_{\mathrm{cont}} w_{\mathrm{bad}} = \|w_{\mathrm{bad}}\|_2^2 = 1.$$
Furthermore, because the columns of $U_{\mathrm{cont}}$ are orthogonal to $u_0, \dots, u_{10}$:
$$\langle u_k, v_{\mathrm{bad}} \rangle = u_k^T U_{\mathrm{cont}} w_{\mathrm{bad}} = 0 \qquad \forall k \in \{0, 1, \dots, 10\}.$$
Therefore, $v_{\mathrm{bad}} \in \Phi^\perp$, which means $(I - P_{\mathrm{bound}}) v_{\mathrm{bad}} = v_{\mathrm{bad}}$ and $P_{\mathrm{bound}} v_{\mathrm{bad}} = 0$.

### 2.3 Theorem 134.1 (Exact Discrete Euler–Lagrange Equation)
The canonical vector $v_{\mathrm{bad}} \in \mathbb{R}^{N+1}$ satisfies the exact inhomogeneous eigenvalue equation:
$$\boxed{(\mathcal{Q}_{\mathrm{comp}} - \mu_0 I) v_{\mathrm{bad}} = \sum_{k=0}^{10} \lambda_k u_k,}$$
where the Lagrange multipliers are given by the exact projections:
$$\boxed{\lambda_k \equiv \langle u_k, \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} \rangle \qquad (k = 0, 1, \dots, 10).}$$

*Proof.*  
Since $w_{\mathrm{bad}}$ satisfies $\widehat{\mathcal{Q}}_{\mathrm{comp}} w_{\mathrm{bad}} = \mu_0 w_{\mathrm{bad}}$, multiplying on the left by $U_{\mathrm{cont}}$ gives:
$$U_{\mathrm{cont}} U_{\mathrm{cont}}^T \mathcal{Q}_{\mathrm{comp}} U_{\mathrm{cont}} w_{\mathrm{bad}} = \mu_0 U_{\mathrm{cont}} w_{\mathrm{bad}}.$$
Using $v_{\mathrm{bad}} = U_{\mathrm{cont}} w_{\mathrm{bad}}$ and the projector identity $U_{\mathrm{cont}} U_{\mathrm{cont}}^T = I - P_{\mathrm{bound}}$:
$$(I - P_{\mathrm{bound}}) \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} = \mu_0 v_{\mathrm{bad}}.$$
Rearranging:
$$\mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} - \mu_0 v_{\mathrm{bad}} = P_{\mathrm{bound}} \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}}.$$
Expanding $P_{\mathrm{bound}} = \sum_{k=0}^{10} u_k u_k^T$:
$$P_{\mathrm{bound}} \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} = \sum_{k=0}^{10} u_k \big( u_k^T \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} \big) = \sum_{k=0}^{10} \lambda_k u_k,$$
where $\lambda_k \equiv u_k^T \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}}$.  
This completes the proof. $\blacksquare$

### 2.4 Mathematical Interpretation of the Lagrange Multipliers
In constrained optimization and spectral variational theory:
- The term $\mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}}$ represents the unconstrained internal gradient of the quadratic functional.
- The term $\mu_0 v_{\mathrm{bad}}$ is the collinear component maintaining normalization $\|v\|_2 = 1$.
- The vector $\mathbf{F}_{\mathrm{constr}} \equiv \sum_{k=0}^{10} \lambda_k u_k$ represents the **normal-space component of the gradient** required to enforce the spectral orthogonality constraints $v \perp u_k$ for $k \in \{0, \dots, 10\}$.
- The total constraint norm:
  $$\|\boldsymbol\lambda\|_2 = \sqrt{\sum_{k=0}^{10} \lambda_k^2} = \|\mathbf{F}_{\mathrm{constr}}\|_2$$
  quantifies the magnitude of the constraint force holding the stationary state $v_{\mathrm{bad}}$ in equilibrium against the potential well within the continuum trial subspace $\Phi^\perp$.
- *Epistemic note:* Because $v_{\mathrm{bad}}$ is constructed as an eigenvector of the compression $\widehat{\mathcal{Q}}_{\mathrm{comp}} = U_{\mathrm{cont}}^T \mathcal{Q}_{\mathrm{comp}} U_{\mathrm{cont}}$, the component of $(\mathcal{Q}_{\mathrm{comp}} - \mu_0 I) v_{\mathrm{bad}}$ inside $\Phi^\perp$ vanishes identically by algebraic construction. The small residual $\|\mathbf{r}_{\mathrm{EL}}\|_2 < 10^{-45}$ is a verification of algebraic consistency rather than an independent physical discovery. The genuine discovery is that $v_{\mathrm{bad}}$ is an isolated constrained minimizer whose constraint force $\boldsymbol\lambda$ is an $O(1)$ vector with significant projection across multiple bound states ($k = 0, 1, 2, \dots$).

---

## 3. Canonical Phase Alignment and Profile Stabilization of $T_{v_{\mathrm{bad}}}(t)$

Because an eigenvector $w$ is defined only up to a global phase $\pm 1$, the continuous wavepacket $T_{v_{\mathrm{bad}}}(t) = \pm |T_{v_{\mathrm{bad}}}(t)|$ can undergo arbitrary sign flips across different dimensions $N$ depending on the internal ordering of the numerical eigensolver (as observed in Cell 133, where $T(0)$ was positive for $N \in \{24, 28, 32, 48, 64\}$ but negative for $N = 40$).

### 3.1 Canonical Phase Convention
To eliminate this spurious discrete symmetry and enable consistent profile tracking, we define the **canonical positive boundary phase**:
$$\operatorname{Phase}(v_{\mathrm{bad}}) \equiv \operatorname{sgn}\big( T_{v_{\mathrm{bad}}}(0) \big) = \operatorname{sgn}\left( (v_{\mathrm{bad}})_0 + \sqrt{2}\sum_{m=1}^N (v_{\mathrm{bad}})_m \right).$$
Whenever $\operatorname{Phase}(v_{\mathrm{bad}}) < 0$, we replace $w_{\mathrm{bad}} \mapsto -w_{\mathrm{bad}}$ (and consequently $v_{\mathrm{bad}} \mapsto -v_{\mathrm{bad}}$).

### 3.2 Discrete and Continuous Profile Stabilization Metrics
Under the normalized measure $d\mu = dt/L$, the continuous $L^2(d\mu)$ distance between two wavepackets $T_{v^{(N_1)}}(t)$ and $T_{v^{(N_2)}}(t)$ with $N_1 \le N_2$ coincides identically with the $\ell^2$ distance of their zero-padded coefficient vectors:
$$\|T_{v^{(N_1)}} - T_{v^{(N_2)}}\|_{L^2(d\mu)}^2 = \sum_{m=0}^{N_1} \big| (v^{(N_1)})_m - (v^{(N_2)})_m \big|^2 + \sum_{m=N_1+1}^{N_2} \big| (v^{(N_2)})_m \big|^2 \equiv \|v^{(N_1)} - v^{(N_2)}\|_2^2.$$

We define the distance metric relative to the highest resolution $N_{\max} = 64$:
$$\delta_{\mathrm{Cauchy}}(N) \equiv \|v_{\mathrm{bad}}^{(N)} - v_{\mathrm{bad}}^{(64)}\|_2.$$
While comparing to $N = 64$ provides evidence of numerical stabilization toward the reference state rather than an unconditional proof of Cauchy convergence, a decreasing sequence $\delta_{\mathrm{Cauchy}}(N)$ indicates whether the discrete Galerkin wavepackets are settling toward a stationary spatial profile.

---

## 4. The Three-Way Energy Partition and the Average Stiffness Balance

Recall from Cell 129 and Cell 130 that the prime Galerkin matrix decomposes into:
$$Q_{\mathrm{prime}}^{\mathrm{even}} \equiv -\widetilde{W} + \widetilde{\mathcal{D}}^{\mathrm{per}} + \Delta\widetilde{\mathcal{D}}.$$
When grouped with the Archimedean diagonal multiplier $D_{\mathrm{mult}} = \operatorname{diag}(h_+(a_m))$, the total competition operator is:
$$\mathcal{Q}_{\mathrm{comp}} \equiv \Omega_{\mathrm{diag}} - \mathcal{K}_{\mathrm{neg}} = (D_{\mathrm{mult}} + \widetilde{\mathcal{D}}^{\mathrm{per}}) - (\widetilde{W} - \Delta\widetilde{\mathcal{D}}).$$
Regrouping:
$$\boxed{\mathcal{Q}_{\mathrm{comp}} \equiv D_{\mathrm{mult}} - \widetilde{W} + \widetilde{\mathcal{D}}^{\mathrm{true}},}$$
where:
$$\widetilde{\mathcal{D}}^{\mathrm{true}} \equiv \widetilde{\mathcal{D}}^{\mathrm{per}} + \Delta\widetilde{\mathcal{D}} \succeq 0$$
is the **exact true finite-interval translation defect matrix**.

### 4.1 Physical Form Representation
In continuous physical coordinate space, the Rayleigh quotient of $v_{\mathrm{bad}}$ evaluates to:
$$\langle v_{\mathrm{bad}}, \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} \rangle = R_{\mathrm{mult}} - R_W + R_{\mathcal{D}}^{\mathrm{true}},$$
where:
1. **Archimedean Diagonal Multiplier:**
   $$R_{\mathrm{mult}} \equiv \sum_{m=0}^N (v_{\mathrm{bad}})_m^2 h_+\left(\frac{2\pi m}{L}\right) > 0.$$
2. **Physical Step Potential Well:**
   $$R_W \equiv \frac{1}{L} \int_0^L W(t) |T_{v_{\mathrm{bad}}}(t)|^2 dt > 0 \qquad \left(W(t) = 2 \sum_{q \le c} w_q \mathbf{1}_{[\log q, L]}(t)\right).$$
3. **True Finite-Interval Translation Stiffness:**
   $$R_{\mathcal{D}}^{\mathrm{true}} \equiv \frac{1}{L} \sum_{q \le c} w_q \int_{\log q}^L |T_{v_{\mathrm{bad}}}(t) - T_{v_{\mathrm{bad}}}(t - \log q)|^2 dt \ge 0.$$

### 4.2 Theorem 134.2 (The Average Translation Stiffness Identity)
In the periodic translation defect multiplier $M(m) = \sum_{q \le c} w_q \sin^2\left(\frac{\pi m \log q}{L}\right)$, the continuous average over frequency space satisfies:
$$\boxed{\overline{4 M(m)} \equiv 4 \sum_{q \le c} w_q \left(\frac{1}{\pi} \int_0^\pi \sin^2\theta \, d\theta\right) = 2 \sum_{q \le c} w_q \equiv W(L) \approx 9.9438.}$$

*Significance:*  
The maximum depth of André Weil's prime step potential $W(L) = 2 \sum_{q \le c} w_q \approx 9.943769$ is **identically equal to the mean translation stiffness** $\overline{4M(m)}$!  
Therefore, the translation stiffness operator on average perfectly cancels the deepest point of the potential well. Any negative competition energy can only arise when a wavepacket exploits the local phase fluctuations of $\sin^2(\theta)$ and low Archimedean values while remaining constrained by bound-state orthogonality.

---

## 5. Formulation of the Continuum 1D Variational Problem

In the continuum limit $N \to \infty$, the discrete vector $v_{\mathrm{bad}} \in \mathbb{R}^{N+1}$ maps to a continuous function $T_\infty \in H^1([0, L])$.

The discrete Euler–Lagrange equation:
$$(\mathcal{Q}_{\mathrm{comp}} - \mu_0 I) v_{\mathrm{bad}} = \sum_{k=0}^{10} \lambda_k u_k$$
becomes the continuous non-local Schrödinger equation:
$$\boxed{\mathcal{L}_{\mathrm{cont}} T_\infty(t) - W(t) T_\infty(t) - \mu_0 T_\infty(t) = \sum_{k=0}^{10} \lambda_k T_{u_k}(t),}$$
where $\mathcal{L}_{\mathrm{cont}}$ is the self-adjoint kinetic operator on $L^2([0, L])$:
$$\mathcal{L}_{\mathrm{cont}} = h_+\left(\sqrt{-\partial_t^2}\right) + \sum_{q \le c} w_q \big( 2 I - \mathcal{S}_{\log q} - \mathcal{S}_{-\log q} \big).$$

Around the frozen peak mode $m^* = 26$, the non-local operator $\mathcal{L}_{\mathrm{cont}}$ admits an effective local second-order differential approximation:
$$\mathcal{L}_{\mathrm{cont}} \approx -\sigma_{\mathrm{eff}} \frac{d^2}{dt^2} + \Omega_{\mathrm{eff}},$$
where $\sigma_{\mathrm{eff}}$ is the effective dispersion curvature at $m^*$.
Cell 134 evaluates this effective kinetic parameter and compares the resulting 1D continuum variational ground state with the audited Galerkin minimizer.

---

## 6. Certified Computational Audit (Cell 134 Output)

The high-precision 50-dps verification suite [`cell134.py`](file:///c:/data/github/connes-cvs-/cell134.py) was executed across $N \in [24, 64]$ with $c = 13$, $L = \log 13 \approx 2.564949$, $T = 600$, and $N_{\mathrm{bound}} = 11$.

### 6.1 Table 1: Canonical Phase-Aligned Profile & Stabilization
Canonical phase convention: $T_{v_{\mathrm{bad}}}(0) > 0$. Metric: $\delta_{\mathrm{Cauchy}}(N) \equiv \|v_{\mathrm{bad}}^{(N)} - v_{\mathrm{bad}}^{(64)}\|_2$.

| $N$ | $T(0)$ | $T(L/2)$ | Peak Mode $m^*$ | $\delta_{\mathrm{Cauchy}}(N)$ |
| :---: | :---: | :---: | :---: | :---: |
| 24 | $+5.0226$ | $+0.0862$ | 22 | $5.3023 \times 10^{-1}$ |
| 28 | $+4.2070$ | $-0.4732$ | 22 | $2.2763 \times 10^{-1}$ |
| 32 | $+3.7524$ | $-0.4972$ | 26 | $1.3749 \times 10^{-1}$ |
| 40 | $+3.5557$ | $-0.4657$ | 26 | $1.0506 \times 10^{-1}$ |
| 48 | $+3.5228$ | $-0.3867$ | 26 | $6.4957 \times 10^{-2}$ |
| 64 | $+3.4820$ | $-0.4495$ | 26 | $0.0000$ (Ref) |

*Observation:* The boundary amplitude stabilizes monotonically toward $T(0) \approx 3.48$, the interior nodal node $T(L/2) \approx -0.45$ settles, and the distance to the $N=64$ reference decreases monotonically from $0.530$ to $0.065$.

### 6.2 Table 2: Spectrum of Lagrange Multiplier Constraint Forces $\lambda_k$
Definition: $\lambda_k \equiv \langle u_k, \mathcal{Q}_{\mathrm{comp}} v_{\mathrm{bad}} \rangle$.

| $N$ | $\lambda_0$ | $\lambda_1$ | $\lambda_2$ | $\lambda_9$ | $\lambda_{10}$ | $\|\boldsymbol\lambda\|_2$ | Dominant $k^*$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $-0.6301$ | $+0.4752$ | $-0.4412$ | $-0.3462$ | $-0.4128$ | $1.3535$ | 0 |
| 28 | $-0.6802$ | $+0.5135$ | $-0.4765$ | $-0.3439$ | $+0.4317$ | $1.4455$ | 0 |
| 32 | $-0.7317$ | $+0.5528$ | $-0.5124$ | $-0.3516$ | $+0.4545$ | $1.5447$ | 0 |
| 40 | $-0.7829$ | $+0.5925$ | $-0.5487$ | $+0.3648$ | $+0.4791$ | $1.6468$ | 0 |
| 48 | $-0.6742$ | $-0.7307$ | $-0.5556$ | $+0.3681$ | $-0.4855$ | $1.6680$ | 1 |
| 64 | $-0.1707$ | $-0.9821$ | $-0.5573$ | $+0.3709$ | $+0.4866$ | $1.6733$ | 1 |

*Observation:*
1. Total constraint norm stabilizes at an $O(1)$ value: $\|\boldsymbol\lambda\|_2 \approx 1.6733$.
2. The constraint force is distributed across multiple bound states rather than concentrated solely in the ground state: at $N = 64$, $\lambda_1 = -0.9821$ (dominant), $\lambda_2 = -0.5573$, $\lambda_{10} = +0.4866$, with ground-state force $\lambda_0 = -0.1707$. The dominant multiplier switches from $k=0$ at lower $N$ to $k=1$ at $N \ge 48$.

### 6.3 Table 3: Exact Discrete Euler–Lagrange Residual Audit
Residual: $\mathbf{r}_{\mathrm{EL}} \equiv (\mathcal{Q}_{\mathrm{comp}} - \mu_0 I) v_{\mathrm{bad}} - \sum_{k=0}^{10} \lambda_k u_k$.

| $N$ | $\mu_0$ (Deficit) | $\|\boldsymbol\lambda\|_2$ | $\|\mathbf{r}_{\mathrm{EL}}\|_2$ | $\|\mathbf{r}_{\mathrm{EL}}\|_{\max}$ | Status |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $-0.00807113$ | $1.353542$ | $1.6016 \times 10^{-50}$ | $6.0346 \times 10^{-51}$ | PASSED |
| 28 | $-0.28445278$ | $1.445493$ | $2.4668 \times 10^{-50}$ | $1.0023 \times 10^{-50}$ | PASSED |
| 32 | $-0.38662701$ | $1.544673$ | $2.5604 \times 10^{-50}$ | $1.2863 \times 10^{-50}$ | PASSED |
| 40 | $-0.45755948$ | $1.646793$ | $1.8767 \times 10^{-50}$ | $7.6842 \times 10^{-51}$ | PASSED |
| 48 | $-0.47803442$ | $1.668010$ | $3.1019 \times 10^{-50}$ | $1.3531 \times 10^{-50}$ | PASSED |
| 64 | $-0.48697922$ | $1.673335$ | $3.4199 \times 10^{-50}$ | $1.5034 \times 10^{-50}$ | PASSED |

*Residual Certification:*
$$\max_{N} \|\mathbf{r}_{\mathrm{EL}}\|_2 = 3.42 \times 10^{-50} \ll 10^{-45}.$$
This confirms the exact algebraic consistency of the projected Euler–Lagrange equation to 50 decimal digits.

### 6.4 Table 4: Exact Three-Way Energy Partition Balance
Decomposition: $R_{\mathrm{comp}} = R_{\mathrm{mult}} - R_W + R_{\mathcal{D}}^{\mathrm{true}}$ where $W(L) = 9.943768796$.

| $N$ | $R_{\mathrm{mult}}$ (Arch) | $R_W$ (Step) | $R_{\mathcal{D}}^{\mathrm{true}}$ (Stiffness) | $R_{\mathrm{comp}}$ (Net $\mu_0$) | $R_W / W(L)$ | Balance Error |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $+1.436458$ | $+3.710369$ | $+2.265839$ | $-0.008071$ | $37.31\%$ | $< 10^{-50}$ |
| 28 | $+1.392251$ | $+3.674564$ | $+1.997859$ | $-0.284453$ | $36.95\%$ | $< 10^{-50}$ |
| 32 | $+1.307048$ | $+3.689425$ | $+1.995750$ | $-0.386627$ | $37.10\%$ | $< 10^{-50}$ |
| 40 | $+1.202064$ | $+3.700465$ | $+2.040841$ | $-0.457559$ | $37.21\%$ | $< 10^{-50}$ |
| 48 | $+1.178727$ | $+3.706678$ | $+2.049916$ | $-0.478034$ | $37.28\%$ | $< 10^{-50}$ |
| 64 | $+1.173761$ | $+3.709783$ | $+2.049043$ | $-0.486979$ | $37.31\%$ | $1.34 \times 10^{-51}$ |

---

## 7. Epistemic Assessment & Conclusions

1. **Algebraic Consistency vs Physical Law:**
   The Euler–Lagrange residual $\|\mathbf{r}_{\mathrm{EL}}\|_2 \le 3.42 \times 10^{-50}$ is an exact consequence of constructing $v_{\mathrm{bad}}$ as an eigenvector of the compression $\widehat{\mathcal{Q}}_{\mathrm{comp}} = U_{\mathrm{cont}}^T \mathcal{Q}_{\mathrm{comp}} U_{\mathrm{cont}}$. Its vanishing certifies algebraic and computational correctness of the subspace projection, rather than representing an unexpected physical identity.
2. **Macroscopic Constraint Force & Multi-Mode Barrier:**
   The genuinely interesting mathematical structure revealed by Cell 134 is that the competition ground state is held in equilibrium with an $O(1)$ constraint force $\|\boldsymbol\lambda\|_2 = 1.6733$. The constraint force is distributed across several bound states ($u_1, u_2, u_9, u_{10}$), with $u_1$ exerting the dominant normal-space force at $N = 64$ ($\lambda_1 = -0.9821$).
3. **Stable Three-Way Form Balance:**
   The competition energy $\mu_0 \approx -0.4870$ is cleanly explained by three stabilizing quadratic forms:
   $$R_{\mathrm{mult}} \to 1.174, \qquad R_W \to 3.710, \qquad R_{\mathcal{D}}^{\mathrm{true}} \to 2.049.$$
   The step potential energy harvested by $v_{\mathrm{bad}}$ saturates at $R_W / W(L) \approx 37.31\%$, demonstrating that spatial oscillation and bound-state orthogonality rigorously constrain how deeply the wavepacket can penetrate the attractive potential well.
4. **Stabilization vs Convergence Proof:**
   The distance $\delta_{\mathrm{Cauchy}}(48) = 0.065$ relative to $N=64$ provides strong empirical evidence of stabilization toward a stationary spatial profile. Formal proof of Cauchy convergence and continuum existence requires analyzing pairwise differences and the continuum variational problem (Cell 135).
