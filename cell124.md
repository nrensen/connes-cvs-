# Cell 124 Analytical Note: Spectral-Subspace Projection, Non-Circularity Resolution, and Variational Min-Max Lower Bounds for the Continuum Threshold

**Companion Document:** Analytical Research Note (Gate 1 / Milestone M-G1.5)  
**Status:** Working Analytical Note (Theoretical Architecture)  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.5: Spectral-Subspace & Min-Max Continuum Threshold](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md); [`cell123.out`](file:///c:/data/github/connes-cvs-/cell123.out); [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md); [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out); [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py)  
**Date:** September 2026  

---

## 1. Executive Summary & Epistemological Resolution of Cell 123

Cell 123 provided a decisive, high-precision empirical audit of the tripartite operator decomposition $Q_{\mathrm{even}} = Q_{\mathrm{arch}} + Q_{\mathrm{prime}} + Q_{\mathrm{pole}}$ and tested the hypothesis that coordinate principal submatrices $C_M = Q_{\mathrm{even}}[M..N, M..N]$ could provide an analytical proof of the continuum spectral threshold $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0$ for Gate 1.

The findings from [`cell123.out`](file:///c:/data/github/connes-cvs-/cell123.out) established two unequivocal falsifications and one major structural discovery:

1. **Definitive Falsification of $M = 3$ Coercivity (Route C):**
   $$\lambda_{\min}(C_3) \approx 4.64 \times 10^{-39} \le E_3 = 2.98 \times 10^{-38} \quad (\text{at } N = 64).$$
   Because the potential well contains $\bar{N}_{\mathrm{bound}} \approx 11$ confined bound states ($E_0, \dots, E_{10} \ll 1$), deleting only 3 coordinate modes leaves the remaining $\sim 8$ bound states supported inside $C_3$. By Cauchy interlacing, $\lambda_{\min}(C_3) \le E_3 \approx 10^{-38}$, rendering $C_3 \succeq 0.386 I$ mathematically impossible.
2. **Definitive Falsification of $M = 12$ Coordinate Coercivity:**
   $$\lambda_{\min}(C_{12}) \text{ collapses by four orders of magnitude: } 3.83 \times 10^{-2} \; (N=16) \longrightarrow 7.87 \times 10^{-6} \; (N=64).$$
   Coordinate truncation to $m \ge 12$ does **not** yield a uniform macroscopic spectral floor $C_{12} \succeq c_{12} I > 0$. The minimum eigenvalue decays toward zero as dimension increases.
3. **Component Indefiniteness & Emergent Cancellation:**
   On high coordinate modes ($M = 12, N = 64$), the Archimedean sector is coercive ($C_{\mathrm{arch}} \succeq 1.553 I$), but the prime sector is strongly indefinite ($\lambda_{\min}(C_{\mathrm{prime}}) \approx -2.373$). The near-positivity of $C_{12}$ ($7.87 \times 10^{-6}$) is produced by near-perfect destructive cancellation between Archimedean and prime distributions, not by componentwise positivity.

### The Foundational Epistemic Lesson
$$\boxed{\text{coordinate-mode truncation } (\operatorname{span}\{e_M, \dots, e_N\}) \ne \text{spectral-subspace projection } (P_{\mathrm{cont}} = I - P_{\mathrm{bound}}).}$$

In a dense Galerkin representation, physical bound-state eigenfunctions $u_j$ have non-compact Fourier frequency support: their Fourier coefficients $u_{j, m}$ decay rapidly (exponentially or super-polynomially) but remain strictly non-zero for all $m \ge M$. Consequently:
$$\operatorname{span}\{e_M, \dots, e_N\} \not\perp \operatorname{span}\{u_0, \dots, u_{10}\}.$$
Projecting a bound state onto the high coordinate subspace leaves a residual vector $w_j = \sum_{m=M}^N u_{j, m} e_m$ whose Rayleigh quotient is of the order of the bound-state energy plus tail leakage:
$$\frac{\langle w_j, Q_{\mathrm{even}} w_j \rangle}{\|w_j\|^2} \approx E_j + \mathcal{O}(\text{tail mass}) \ll 1.$$
This creates an artificial near-zero Ritz eigenvalue in $C_M$. Cauchy interlacing on coordinate blocks is therefore fundamentally incapable of separating the bound-state sector from the continuum. All coordinate-submatrix $M$-sweeps are permanently retired.

### The Solid Foundation (Cell 121)
The failure of coordinate submatrix coercivity does **not** affect the empirical existence of the physical continuum gap discovered in Cell 121:
$$g_{2, 12}(N) \equiv E_{13}^{(N)} - E_3^{(N)} \approx 0.582 \quad (\text{at } N = 64).$$
Gate 1 relies exclusively on the eigenvalues of the **full operator**, not on the properties of its coordinate submatrices. While empirical across tested discrete dimensions, the continuum gap is an observed property of the full spectrum.

---

## 2. The Spectral-Subspace Projection Formalism

To isolate the continuum without suffering from coordinate tail leakage, the projection must be performed onto the **spectral subspace** spanned by the physical eigenfunctions.

### Definition 2.1 (Bound-State & Continuum Spectral Projectors)
Let $Q_{\mathrm{even}}^{(N)} \in \mathbb{R}^{(N+1) \times (N+1)}$ be the canonical even Galerkin matrix, with ordered orthonormal eigenbasis $\{u_j^{(N)}\}_{j=0}^N$ and eigenvalues $E_0^{(N)} \le E_1^{(N)} \le \dots \le E_N^{(N)}$.

For cutoff parameter $c = 13, T = 600$, current Galerkin calculations show an observed low-energy cluster of approximately $\bar{N}_{\mathrm{bound}} \approx 11$ states with energies $E_0, \dots, E_{10} \le 10^{-5}$, separated from higher states with $E_{11} \approx 0.575$ and $E_{12} \approx 0.582$ at $N = 64$. (Note: 11 is the empirical count from finite-dimensional calculations, not yet an analytically proved theorem).

1. The **bound-state spectral projector** is defined as:
   $$P_{\mathrm{bound}}^{(N)} \equiv \sum_{j=0}^{\bar{N}_{\mathrm{bound}}-1} u_j^{(N)} (u_j^{(N)})^T.$$
2. The **continuum spectral projector** is the orthogonal complement:
   $$P_{\mathrm{cont}}^{(N)} \equiv I - P_{\mathrm{bound}}^{(N)} = \sum_{j=\bar{N}_{\mathrm{bound}}}^N u_j^{(N)} (u_j^{(N)})^T.$$
3. The **continuum-projected operator** is:
   $$Q_{\mathrm{cont}}^{(N)} \equiv P_{\mathrm{cont}}^{(N)} Q_{\mathrm{even}}^{(N)} P_{\mathrm{cont}}^{(N)} = \sum_{j=\bar{N}_{\mathrm{bound}}}^N E_j^{(N)} u_j^{(N)} (u_j^{(N)})^T.$$

### Proposition 2.1 (Exact Spectral Coercivity on the Continuum Subspace)
Let $\mathcal{H}_{\mathrm{cont}}^{(N)} \equiv \operatorname{Ran}(P_{\mathrm{cont}}^{(N)}) = \{v \in \mathbb{R}^{N+1} : \langle u_j^{(N)}, v \rangle = 0, \;\; \forall j < \bar{N}_{\mathrm{bound}}\}$.
Then for every unit vector $v \in \mathcal{H}_{\mathrm{cont}}^{(N)}$ ($\|v\|_2 = 1$):
$$\langle v, Q_{\mathrm{even}}^{(N)} v \rangle = \langle v, Q_{\mathrm{cont}}^{(N)} v \rangle \ge E_{\bar{N}_{\mathrm{bound}}}^{(N)}.$$
Consequently, on its invariant subspace $\mathcal{H}_{\mathrm{cont}}^{(N)}$, $Q_{\mathrm{cont}}^{(N)}$ satisfies:
$$\lambda_{\min}\left(Q_{\mathrm{cont}}^{(N)}\big|_{\mathcal{H}_{\mathrm{cont}}^{(N)}}\right) \equiv E_{\bar{N}_{\mathrm{bound}}}^{(N)}.$$

*Proof.* Expanding $v = \sum_{j=\bar{N}_{\mathrm{bound}}}^N c_j u_j^{(N)}$ with $\sum |c_j|^2 = 1$:
$$\langle v, Q_{\mathrm{even}}^{(N)} v \rangle = \sum_{j=\bar{N}_{\mathrm{bound}}}^N E_j^{(N)} |c_j|^2 \ge E_{\bar{N}_{\mathrm{bound}}}^{(N)} \sum_{j=\bar{N}_{\mathrm{bound}}}^N |c_j|^2 = E_{\bar{N}_{\mathrm{bound}}}^{(N)}.$$
Because every bound-state component has been projected out orthogonally ($\langle u_j, v \rangle \equiv 0$ for $j < \bar{N}_{\mathrm{bound}}$), there is zero bound-state tail contamination. $\blacksquare$

*Epistemic Remark on Proposition 2.1:*
Proposition 2.1 is mathematically exact, but as recognized in §3, it is essentially a restatement of the spectral theorem. It identifies the exact mathematical object whose positivity is needed, but it does **not** constitute an independent proof of continuum positivity, because $P_{\mathrm{cont}}^{(N)}$ is defined using the exact eigenvectors $u_j^{(N)}$.

---

## 3. The Non-Circularity Dilemma

While Proposition 2.1 is mathematically exact, it exposes the central analytical challenge of Gate 1:

> **The Non-Circularity Question:**  
> The projector $P_{\mathrm{cont}}^{(N)}$ is constructed from the exact eigenvectors $u_0^{(N)}, \dots, u_{10}^{(N)}$ of $Q_{\mathrm{even}}^{(N)}$. If one must compute or assume the spectral decomposition of $Q_{\mathrm{even}}^{(N)}$ in order to define $P_{\mathrm{cont}}^{(N)}$, one cannot use coercivity of $P_{\mathrm{cont}} Q P_{\mathrm{cont}}$ as an *independent proof* that $E_{\bar{N}_{\mathrm{bound}}}^{(N)} \ge c_* > 0$.

To establish an audit-proof analytical proof of the continuum threshold, we must characterize the continuum subspace and lower-bound the energy $E_{\bar{N}_{\mathrm{bound}}}^{(N)}$ **directly from the quadratic form**, without presupposing the full spectral theorem for the discrete matrix.

We formulate three constructive analytical pathways to resolve this dilemma:
- **Route 1:** Courant–Fischer Min-Max with Adapted Quasimode Codimension Subspaces.
- **Route 2:** Dunford–Schwartz Resolvent Contour Integral Projection.
- **Route 3:** Continuum Spatial Delocalization and Riemann Phase Cancellation against the Prime Dirac Comb.

---

## 4. Route 1: Courant–Fischer Min-Max with Adapted Quasimodes

The Courant–Fischer (min-max) characterization of eigenvalues provides a variational mechanism to lower-bound the continuum threshold without computing exact eigenvectors.

### Theorem 4.1 (Courant–Fischer Codimension Lower Bound)
Let $Q \in \mathbb{R}^{(N+1) \times (N+1)}$ be a real symmetric matrix with eigenvalues $E_0 \le E_1 \le \dots \le E_N$.
For any integer $K \ge 1$, the $(K+1)$-th eigenvalue satisfies:
$$E_K = \max_{\operatorname{codim}(V) = K} \min_{v \in V, \|v\|=1} \langle v, Q v \rangle = \max_{w_0, \dots, w_{K-1} \in \mathbb{R}^{N+1}} \min_{\substack{v \perp w_0, \dots, w_{K-1} \\ \|v\|=1}} \langle v, Q v \rangle.$$
Consequently, for **any choice** of $K$ trial vectors $\{\phi_0, \dots, \phi_{K-1}\}$:
$$E_K \ge \min_{\substack{v \perp \phi_0, \dots, \phi_{K-1} \\ \|v\|=1}} \langle v, Q v \rangle.$$

*Significance:* This inequality is **unconditional**. We do not need $\{\phi_k\}$ to be exact eigenvectors. If we can find a family of analytically constructed trial vectors $\{\phi_k\}$ such that the quadratic form $\langle v, Q v \rangle$ is bounded below by $c_* > 0$ on the orthogonal complement $V_\Phi = \{\phi_0, \dots, \phi_{K-1}\}^\perp$, then:
$$E_K \ge c_* > 0$$
follows rigorously and unconditionally.

### Construction 4.1 (Analytic Well Quasimodes)
In the semiclassical analysis of Paper NR2 (§4–§5), the potential well $V_{\mathrm{conf}}(t)$ confined between $t = 0$ and $t = L$ creates bound states that are well-approximated by:
1. **The Ground State:** The solitary wave profile $T_\infty(t) = C t (L - t) e^{-\dots}$ with dual Dirichlet vanishing $T(0) = T(L) = 0$.
2. **Excited Bound States:** Parabolic cylinder / Weber functions (or localized Hermite–Gauss wavepackets) with $k$ nodes inside the well $[0, L]$:
   $$\phi_k(t) = \mathcal{N}_k H_k\left(\sqrt{\omega}\left(t - \frac{L}{2}\right)\right) e^{-\frac{\omega}{2}(t - L/2)^2}, \qquad k = 0, \dots, K-1.$$
3. **Galerkin Projection of Quasimodes:**
   The Fourier coefficients of these analytic functions form discrete trial vectors $\boldsymbol\phi_k \in \mathbb{R}^{N+1}$:
   $$(\boldsymbol\phi_k)_m = \int_0^L \phi_k(t) \sqrt{\frac{2}{L}} \cos\left(\frac{2\pi m t}{L}\right) dt.$$

### Heuristic Lemma 4.1 (Quasimode Perturbation Conjecture — Target for Cell 125)
Let $\Phi = \operatorname{span}\{\boldsymbol\phi_0, \dots, \boldsymbol\phi_{K-1}\}$ be an orthonormal family of analytic quasimodes satisfying the approximate eigenvalue equation:
$$\|Q_{\mathrm{even}} \boldsymbol\phi_k - \lambda_k \boldsymbol\phi_k\|_2 \le \varepsilon_k \ll 1, \qquad 0 \le \lambda_k \le \lambda_{\max}^{\mathrm{bound}} \ll 1.$$
Let $V_\Phi = \Phi^\perp$. It is tempting to write a Temple/Davis–Kahan-type estimate:
$$\langle v, Q_{\mathrm{even}} v \rangle \ge E_K - \sum_{k=0}^{K-1} \frac{\varepsilon_k^2}{E_K - \lambda_k} \quad (\forall v \in V_\Phi, \|v\|=1).$$
*Critical Epistemic Gap & Cell 125 Mandate:*  
The formula above is an **unverified heuristic assertion**. The mere fact that $\|Q \boldsymbol\phi_k - \lambda_k \boldsymbol\phi_k\| \le \varepsilon_k$ does not automatically imply this lower bound on every vector orthogonal to all $\boldsymbol\phi_k$ without explicit hypotheses on mutual orthogonality, residual projection, and spectral gaps.  
Deriving a mathematically watertight variational lower bound on $\Phi^\perp$ directly from first principles—or replacing it with invariant-subspace separation theorems ($\sin\Theta$ machinery)—is the primary analytical mandate of **Cell 125**.

---

## 5. Route 2: Dunford–Schwartz Resolvent Contour Projection

Instead of constructing eigenvectors individually, the spectral projection onto the bound-state sector can be expressed as a global Cauchy contour integral in the complex plane.

### Theorem 5.1 (Riesz Resolvent Projection)
Let $Q \equiv Q_{\mathrm{even}}^{(N)}$. Let $\Gamma \subset \mathbb{C}$ be a positively oriented closed Jordan curve enclosing the bound-state spectrum $\{E_0, \dots, E_{K-1}\}$ while enclosing no other points of $\sigma(Q)$.
Then the bound-state spectral projector is given exactly by:
$$P_{\mathrm{bound}} = \frac{1}{2\pi i} \oint_\Gamma (z I - Q)^{-1} dz.$$
The continuum projector is:
$$P_{\mathrm{cont}} = I - \frac{1}{2\pi i} \oint_\Gamma (z I - Q)^{-1} dz.$$

### Structure of the Spectral Gap in the Complex Plane
From the numerical audits of Cells 121 and 123:
- The bound-state cluster satisfies: $\sigma_{\mathrm{bound}}(Q) \subset [0, 10^{-5}]$ (observed).
- The continuum spectrum satisfies: $\sigma_{\mathrm{cont}}(Q) \subset [0.575, \infty)$ (observed).
- The spectral gap is macroscopic:
  $$\mathrm{dist}(\sigma_{\mathrm{bound}}, \sigma_{\mathrm{cont}}) = E_{11} - E_{10} \approx 0.57557 > 0.57.$$

Choose the circular contour:
$$\Gamma = \{z \in \mathbb{C} : |z| = r_0\}, \qquad r_0 \equiv 0.20.$$
Along this contour:
1. Distance to the bound states: $|z - E_j| \ge r_0 - 10^{-5} \approx 0.20$.
2. Distance to the continuum: $|z - E_k| \ge 0.575 - 0.20 = 0.375$.
3. Resolvent norm bound:
   $$\|(z I - Q)^{-1}\|_{\mathrm{op}} = \frac{1}{\mathrm{dist}(z, \sigma(Q))} \le \frac{1}{0.20} = 5.0 \quad (\forall z \in \Gamma).$$

### Operational Role: Auxiliary Stability, Not Primary Gap Proof
Route 2 is best regarded as a **stability and perturbation tool once a spectral gap is independently established**, rather than the primary mechanism that creates the gap. Approximating $P_{\mathrm{bound}}$ via resolvent identities requires independent control over the resolvent of the unperturbed or limiting operator $\|(z I - Q_0)^{-1}\|$ on the same contour $\Gamma$.

---

## 6. Route 3: Speculative Heuristic — Prime Phase Interference

Cell 123 revealed that high-mode coercivity fails in coordinate space because the prime operator $Q_{\mathrm{prime}}$ is **strongly indefinite** ($\lambda_{\min}(C_{\mathrm{prime}}) \approx -2.373$), canceling almost completely with the Archimedean piece ($C_{\mathrm{arch}} \ge 1.553$).

### Informal Physical Heuristic (Outside Active Proof Pipeline)
In coordinate space, the prime operator is an explicit sum over prime powers with discrete spikes at $\log(p^k)$.  
It is physically plausible that delocalized scattering states oscillating at high wavenumber $k \gg 1$ experience destructive phase cancellation against these prime spikes:
$$\sum_{p^k \le e^L} \frac{\Lambda(p^k)}{p^{k/2}} \cos(k_j \log(p^k) + \theta) \stackrel{?}{=} \mathcal{O}\left(k_j^{-1/2} \log k_j\right) \longrightarrow 0.$$
*Epistemic Warning:*  
This is a **speculative number-theoretic heuristic**, not a proven analytical bound. Asserting sub-unitary bounds on prime exponential sums without proof violates repository standards. Route 3 is explicitly **excluded from the active proof pipeline** and serves only as an intuitive qualitative picture.

---

## 7. Strategic Conclusions & Forward Path for Gate 1

### Summary of Cell 124 Analytical Findings
1. **Coordinate Submatrix Coercivity is Officially Retired:** Deleting coordinate modes $e_0, \dots, e_{M-1}$ leaves wavepacket tail leakage $\sim \mathcal{O}(\text{tail mass})$, producing near-zero Ritz eigenvalues ($\lambda_{\min}(C_{12}) \approx 7.87 \times 10^{-6}$).
2. **Spectral Projection is Exactly Coercive:** On $\mathcal{H}_{\mathrm{cont}} = \operatorname{Ran}(I - P_{\mathrm{bound}})$, $Q_{\mathrm{cont}}$ has spectral floor identically equal to $E_{\bar{N}_{\mathrm{bound}}}$ by the spectral theorem. This identifies the target subspace but does not prove positivity independently.
3. **Variational Min-Max is the True Analytical Route:** Through Courant–Fischer (Theorem 4.1), any trial subspace $\Phi$ provides an unconditional lower bound $E_K \ge \min_{v \in \Phi^\perp, \|v\|=1} \langle v, Q v \rangle$.
4. **The Target for Cell 125:**
   $$\boxed{\text{We now need to prove a variational lower bound, not define a spectral projector.}}$$

### Gate 1 Target Proposition 1.0 Pipeline
Gate 1 separates into two decoupled problems:
- **Problem A (Low-Energy Sector Removal):** Construct an analytic $K$-dimensional subspace $\Phi$ and prove $\langle v, Q_{\mathrm{even}} v \rangle \ge c_* > 0$ on $\Phi^\perp$ (Cell 125).
- **Problem B (Tunneling Decay):** Establish $\Delta_j(N) \le C_j e^{-\sigma_j N}$ via dense-operator Agmon metrics (Route 1A).
Importantly, Problem A does not depend on the boundary-defect decay rate $\alpha_N$, insulating Gate 1 from boundary-defect collapse.

---

## References
- [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md) — Operator Decomposition Audit & High-Mode Coercivity
- [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md) — Bound-State-to-Continuum Transition & Core-Size Invariant Mapping
- [`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md) — Strategic Roadmap & Gate 1 Milestones M-G1.4 & M-G1.5
- [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) — Semiclassical Barrier Mechanics & Asymptotic Weil Positivity
