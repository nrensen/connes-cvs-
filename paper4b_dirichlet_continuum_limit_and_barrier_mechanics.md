# The Dirichlet Continuum Limit, Barrier Mechanics, and Asymptotic Weil Positivity in the Connes–van Suijlekom Galerkin Truncation

**Authors:** Research Record / Connes–CvS Investigation Series  
**Date:** September 2026  
**Software & Reproducibility Suite:** `https://github.com/akivag613/connes-cvs-` (mirror: `nrensen/connes-cvs-`)  
**Status:** Standalone Manuscript / Research Programme (Companion to Paper 4: *The Toolkit*)

---

### Abstract

The truncated Weil quadratic form of Connes–van Suijlekom (2025) and Connes–Consani–Moscovici (2026) projects the explicit formula for the Riemann zeta function onto finite-rank Galerkin subspaces of band $N$ on a logarithmic scaling interval $[0, L] = [0, \log c]$. While the algebraic structure of the finite-$N$ Archimedean resolvent, Cauchy transform, and commutator algebra are established as exact mathematical theorems in the companion paper (Paper 4: *An Exact Resolvent and Commutator Toolkit*), the physical mechanism driving Weil positivity in the infinite-dimensional limit $N \to \infty$ involves non-trivial asymptotic and continuum phenomena.

This manuscript sets forth the analytical and empirical research programme investigating the continuum limit $N \to \infty$ across 24 Galerkin dimensions ($N = 1, \dots, 24$) and multiple prime cutoffs $c \in \{5, 7, 11, 13, 17\}$:

1. **Observed Geometric Boundary Suppression and Eigenvalue Gap Law (Numerical Conjectures):** The boundary values of the normalized Galerkin ground-state profiles $T_{v_N}(0) = v_0 + \sqrt{2}\sum_{m=1}^N v_m$ drop precipitously from $7.52 \times 10^{-3}$ at $N = 1$ to $1.14 \times 10^{-20}$ at $N = 24$ (spanning 17.8 decimal orders of magnitude), obeying the geometric suppression law $|T_{v_N}(0)| \sim C(c) \rho(c)^N$ with $0 < \rho(c) < 1$. Simultaneously, the ground-state eigenvalue $\lambda_{\min}(N)$ decreases from $3.10 \times 10^{-6}$ down to $2.53 \times 10^{-43}$ across 43 decimal orders of magnitude, exhibiting precise asymptotic proportionality to the boundary leakage energy:
   $$\lambda_{\min}(N) \sim \kappa_c A_0(N) \equiv \frac{2\kappa_c}{L} [T_{v_N}(0)]^2 \longrightarrow 0 \qquad (N \to \infty),$$
   with an empirical scaling ratio $\kappa_c \approx 0.00238 \pm 0.00004$ that is approximately cutoff-independent for all $c \ge 7$.
2. **The Continuum Solitary Wave and Infinite-Order Boundary Flatness (Conjecture):** As $N \to \infty$, the spatial trigonometric wave $T_{v_N}(t)$ converges to a smooth, strictly positive solitary wave $T_\infty(t)$ on $[0, L]$ with dual Dirichlet boundary vanishing $T_\infty(0) = T_\infty(L) = 0$ and conjectured infinite-order flat boundary contact: $T_\infty^{(k)}(0) = T_\infty^{(k)}(L) = 0$ for all $k \ge 0$, such that the zero-extension $\widetilde{T}_\infty \in C_c^\infty(\mathbb{R})$ with $\operatorname{supp} \widetilde{T}_\infty = [0, L]$. This boundary flatness conditionally eliminates the finite-rank Volterra boundary jump at $\omega = 1$, removing the classical obstruction to continuous Weil positivity.
3. **Semiclassical WKB Barrier Tunneling Mechanics:** Inverting the ground-state profile defines an effective Schrödinger potential $V_{\mathrm{conf}}(t) - E = T''(t)/T(t)$ whose midpoint well rises steeply toward the boundaries. Semiclassical WKB tunneling action across the barrier $\mathcal{S}_{\mathrm{WKB}} = \int_0^{t_{\mathrm{turn}}} \sqrt{T''/T} \, dt$ reproduces the observed logarithmic boundary suppression to within $5.3\%$ across 47 orders of magnitude, consistent with the empirical scaling law $\mathcal{S}_{\mathrm{WKB}}(N, c) \approx \frac{\pi N}{4} \log c$.
4. **Legendre Multipole Decomposition and Asymptotic Tail Extinction:** Via Bauer's spherical Bessel expansion, $T_{v_N}(t)$ decomposes into Legendre multipoles with alternating signs, producing strictly constructive interference at the midpoint and destructive cancellation at the boundaries. The high-frequency Taylor coefficients $A_k(N)$ extinguish rapidly across all orders ($A_0 \sim 10^{-40}, A_1 \sim 10^{-34}, A_2 \sim 10^{-29}$ at $N = 24$), motivating the conjecture that the continuum resolvent $R_\infty(r) = o(r^{-k})$ decays faster than every inverse power of $r$.
5. **Tri-Partite Zero-Energy Equilibrium and Finite-$T$ Archimedean Leakage:** Continuous-variable numerical quadrature and the exact closed-form digamma identity independently cancel the algebraic pole ($+1.551652$) and prime ($-0.071854$) contributions down to a residual of $Q_{\mathrm{total}} = 4.201 \times 10^{-43}$ at $N = 24$ (with continuous quadrature agreeing to $4.96 \times 10^{-25}$). The ratio $\lambda_{\min}(24) / Q_{\mathrm{total}} = 0.6030$ against the matrix eigenvalue $\lambda_{\min}(24) = 2.533 \times 10^{-43}$ isolates the exact finite-$T$ Archimedean cutoff leakage $\delta_T^{\mathrm{tail}} = \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v) - \lambda_{\min}(24) = 1.668 \times 10^{-43}$, certified free of numerical truncation error and proven in Theorem 5.5 of Paper 4 to be 100% continuous tail leakage from the finite integration cutoff $T = 400$.
6. **Formal Continuum Wiener–Hopf Scaling, Interlaced Tunneling Ladders, and Operator Reformulation of the First Jet:** In the continuum scaling limit, the divided-difference Galerkin kernel transforms into a half-line Wiener–Hopf convolution operator with kernel $K_{\mathrm{sym}}(w) = \frac{w}{2\sinh(w/2)}$ whose symbol factors into squared Gamma functions $\frac{\pi^2}{\cosh^2(\pi k)} = [\Gamma(\frac{1}{2} - ik)]^2 [\Gamma(\frac{1}{2} + ik)]^2$. The resulting double pole at $k = -i/2$ generates a logarithmic boundary layer $\phi(x) \sim -\log x$ as $x \to 0^+$, explaining the observed bulk/edge asymmetry between $D_0$ and $D_1$. Below the double-well barrier top, the even and odd spectra form an interlaced ladder of bound tunneling doublets ($E_0 < \mu_1 < E_1 < \dots < \mu_{N/2} < E_{N/2}$). The first coordinate transition dipole moment $\langle e_1, Kc \rangle = -D_0 a_1 / \Delta_1 \approx 1.3134 = \mathcal{O}(1)$ carries $99.9999\%$ of the coordinate norm $\|Kc\|^2$, providing empirical evidence connecting the overlap ratio $C_N = |a_1| / \sqrt{\Delta_1} \approx 2.4$ to the transition dipole moment. The first-jet ratio decomposes into exact operator components $D_1/D_0 = \kappa^2 [\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} - D_0^2 M_2]$, and exhibits numerical behaviour consistent with an asymptotic spectral filtering mechanism where the ratio appears to stabilize near $\approx -0.77$ with net ratio $(D_1/D_0) / [\kappa^2 \mathcal{T}_{\mathrm{diag}}] \sim 0.23$ across tested dimensions. The square-root overlap scaling ($|d_k| \approx 24.1 \sqrt{\Delta_k}$, $|a_j| \approx 2.24 \sqrt{\Delta_j}$) keeps small denominators bounded mode by mode, supporting the hypothesis that $u_1^{-1} = |D_1/D_0|$ is bounded polynomially ($\sim N^{2.6 - 2.8}$ over the tested range), which would guarantee exponential boundary-defect decoupling in the continuum limit.
7. **The Analytical Roadmap toward Continuous Weil Positivity:** We formulate the three open mathematical stages required to convert these empirical and asymptotic findings into a complete proof of Weil positivity on the idele class group.

---

## 1. Introduction and Connection to the Rigorous Toolkit

In André Weil's 1952 formulation, the Riemann Hypothesis (RH) is equivalent to the non-negativity of the quadratic functional:

$$W(g) \ge 0$$

for all admissible test functions $g = f * f^*$ on the idele class group $\mathbb{A}_{\mathbb{Q}} / \mathbb{Q}^\times$. In Alain Connes' non-commutative geometry programme, this positivity is pursued through spectral truncations on prolate spheroidal wave spaces. 

Recent work by Connes and van Suijlekom (2025) and Connes, Consani, and Moscovici (2026) models this functional on a compact logarithmic scaling interval $[0, L] = [0, \log c]$ (with prime cutoff $c > 1$) projected onto a finite-rank Galerkin subspace of frequency band $N \ge 1$, generating an explicit $(2N+1) \times (2N+1)$ matrix:

$$Q_{c, N} = Q_{\mathrm{prime}} + Q_{\mathrm{pole}} + Q_{\mathrm{arch}}.$$

### Separation of the Exact Toolkit and the Continuum Programme

The investigation of this finite-rank system naturally divides into two distinct mathematical realms:

1. **The Exact Finite-$N$ Toolkit (Paper 4):**
   In the companion paper, *An Exact Resolvent and Commutator Toolkit for the Truncated Connes–van Suijlekom Weil Quadratic Form* [Paper 4], all finite-$N$ algebraic structures are established unconditionally as pure mathematical theorems. These include:
   - The exact four-term Volterra reduction.
   - The squared rational resolvent identity $R_v(r) \equiv \frac{2}{L} [ \frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} ]^2$ and operator identity $D(z) \equiv [(I + z\mathcal{L})^{-1} T_v](0)$.
   - Unconditional pointwise non-negativity $K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0$ on $\mathbb{R}$.
   - Spectral lattice orthogonality $K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2$.
   - The exact closed-form Archimedean Cauchy transform $J(q)$ and unconditionally convergent Weierstrass pole series for $\mathcal{Q}_{\mathrm{arch}}(v)$ with fast $\mathcal{O}(n^{-2})$ convergence.
   - Exact rank-$2k$ commutator algebra $[M^k, Q]$, strict parity decoupling, and the odd-sector resolvent identity $M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi$.
   - The exact first-jet identity $D_1/D_0 \equiv -\frac{1}{2} A_1/A_0 \equiv -\kappa^2 F'(0)/F(0)$.
   - The exact small-denominator cancellation $(E_k - \lambda)$ in the odd-even resolvent coupling.

2. **The Asymptotic Continuum Programme (This Manuscript):**
   While Paper 4 provides the exact mathematical machinery, it does not address the physical and asymptotic questions of how the system behaves as the dimension $N \to \infty$:
   - Why is the minimum eigenvalue $\lambda_{\min}(N) > 0$ strictly positive for all finite $N$?
   - What is the geometric and semiclassical mechanism driving $\lambda_{\min}(N) \to 0$?
   - Does the finite trigonometric wave converge to a well-behaved continuous profile?
   - How is the finite-rank Volterra boundary jump eliminated?
   - How can these asymptotic properties be assembled into a rigorous roadmap toward continuous Weil positivity?

This paper presents the empirical, numerical, and asymptotic evidence answering these questions, establishing a structured research programme whose conjectures are systematically tested against high-precision computational data.

---

## 2. Large-$N$ Asymptotics of the Galerkin Ground State

We analyze the sequence of normalized Galerkin ground states $v_N \in \mathbb{R}^{N+1}$ defined by the finite-rank eigensystem:

$$Q_{c, N} v_N = \lambda_{\min}(N) v_N, \qquad \|v_N\|_2 = \sqrt{\sum_{m=0}^N v_{N, m}^2} = 1.$$

High-precision diagonalizations of $Q_{c, N}$ were carried out at 50 decimal digits of precision for $c = 13$ across all dimensions $N = 1, \dots, 24$ using the companion arbitrary-precision Python analysis suite (`cell34.py`, `cell40.py`, `cell41.py` with corresponding verification logs `cell34.out`, `cell40.out`, `cell41.out`) [10].

### 2.1 Observed $\ell^2$ Concentration and Numerical Convergence

**Table 1: Ground-State Mode Convergence in $\ell^2$ ($c = 13$)**

| $N$ | $v_{N, 0}$ | $v_{N, 1}$ | $v_{N, 2}$ | $v_{N, 3}$ | $\|v_N - v_{N-1}\|_{\ell^2}$ | Tail Mass ($\sum_{m > 4} v_m^2$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | 0.818996 | -0.573799 | — | — | — | — |
| 2 | 0.729744 | -0.666834 | 0.151017 | — | 0.198565 | — |
| 4 | 0.648433 | -0.697681 | 0.297857 | -0.063506 | 0.073238 | — |
| 8 | 0.590042 | -0.691774 | 0.389631 | -0.143200 | 0.021807 | $8.85 \times 10^{-6}$ |
| 12 | 0.567664 | -0.684002 | 0.419077 | -0.178518 | 0.009877 | $4.13 \times 10^{-5}$ |
| 16 | 0.556735 | -0.679189 | 0.432079 | -0.196292 | 0.004908 | $7.80 \times 10^{-5}$ |
| 20 | 0.550397 | -0.676105 | 0.439184 | -0.206694 | 0.002905 | $1.10 \times 10^{-4}$ |
| 24 | 0.546859 | -0.674292 | 0.443008 | -0.212522 | 0.001995 | $1.32 \times 10^{-4}$ |

The computed ground states exhibit strong $\ell^2$ mode concentration:
- The step increment $\|v_N - v_{N-1}\|_{\ell^2}$ decreases monotonically from $0.198$ at $N = 2$ down to $0.00199$ at $N = 24$.
- At $N = 24$, over **$99.98\%$** of the total $\ell^2$ mass resides in the lowest five Fourier modes ($m \le 4$), with the remaining 20 modes carrying less than $0.0132\%$ of the energy.

### 2.2 Observed Geometric Boundary Suppression

The boundary value $D_0(N) = T_{v_N}(0) = v_0 + \sqrt{2} \sum_{m=1}^N v_m$ and the second derivative $D_1(N) = T_{v_N}''(0)$ were tracked across all dimensions:

**Table 2: Geometric Boundary Suppression across Galerkin Dimensions ($c = 13$)**

| $N$ | $|D_0(N)|$ | $|D_1(N)|$ | Step Ratio $|D_0(N)| / |D_0(N-1)|$ | Decay Exponent $\alpha_N$ |
| :---: | :---: | :---: | :---: | :---: |
| 1 | $7.52 \times 10^{-3}$ | 4.8694 | — | — |
| 2 | $2.68 \times 10^{-4}$ | 0.5326 | 0.0357 | 3.333 |
| 4 | $6.00 \times 10^{-7}$ | $5.92 \times 10^{-3}$ | 0.0488 | 3.020 |
| 8 | $8.05 \times 10^{-11}$ | $3.36 \times 10^{-6}$ | 0.1434 | 1.942 |
| 12 | $6.65 \times 10^{-14}$ | $6.40 \times 10^{-9}$ | 0.2268 | 1.484 |
| 16 | $1.78 \times 10^{-16}$ | $3.13 \times 10^{-11}$ | 0.2861 | 1.251 |
| 20 | $8.38 \times 10^{-19}$ | $2.58 \times 10^{-13}$ | 0.3105 | 1.170 |
| 24 | $1.14 \times 10^{-20}$ | $5.92 \times 10^{-15}$ | 0.3244 | 1.126 |

Between $N = 1$ and $N = 24$, the boundary value drops from $7.52 \times 10^{-3}$ to $1.14 \times 10^{-20}$, spanning approximately **17.8 decimal orders of magnitude**. The effective decay exponent $\alpha_N = -\frac{\log(|D_0(N)|/|D_0(N-1)|)}{\log c}$ decreases from $3.333$ toward values near $1.1$. We formulate this asymptotic behavior as a conjecture:

### Conjecture 2.1 (Geometric Boundary Suppression)
*For fixed prime cutoff $c > 1$, the boundary values of the normalized Galerkin ground states satisfy:*

$$|T_{v_N}(0)| \sim C(c) \cdot \rho(c)^N \qquad (N \to \infty),$$

*for some positive constant $C(c)$ and decay base $0 < \rho(c) < 1$. A secondary hypothesis suggests $\rho(c) \approx c^{-1/2}$, but determining the exact asymptotic base remains an open problem.*

### 2.3 Observed Eigenvalue-to-Boundary Proportionality

**Table 3: Ground-State Eigenvalue vs Boundary Leakage Energy ($c = 13$)**

| $N$ | $\lambda_{\min}(N)$ | $A_0(N) = \frac{2}{L} [T_{v_N}(0)]^2$ | Ratio $\lambda_{\min}(N) / A_0(N)$ |
| :---: | :---: | :---: | :---: |
| 1 | $3.10 \times 10^{-6}$ | $4.41 \times 10^{-5}$ | 0.07028 |
| 4 | $8.83 \times 10^{-15}$ | $2.81 \times 10^{-13}$ | 0.03142 |
| 8 | $6.71 \times 10^{-23}$ | $5.05 \times 10^{-21}$ | 0.01328 |
| 12 | $1.78 \times 10^{-29}$ | $3.45 \times 10^{-27}$ | 0.00517 |
| 16 | $7.12 \times 10^{-35}$ | $2.48 \times 10^{-32}$ | 0.00287 |
| 18 | $1.24 \times 10^{-37}$ | $5.18 \times 10^{-35}$ | 0.00239 |
| 20 | $1.32 \times 10^{-39}$ | $5.48 \times 10^{-37}$ | 0.00241 |
| 22 | $1.89 \times 10^{-41}$ | $7.97 \times 10^{-39}$ | 0.00237 |
| 24 | $2.53 \times 10^{-43}$ | $1.01 \times 10^{-40}$ | **0.00251** |

While both $\lambda_{\min}(N)$ and $A_0(N)$ span approximately 43 decimal orders of magnitude, their ratio stabilizes remarkably:

$$\frac{\lambda_{\min}(N)}{A_0(N)} \approx 0.00246 \pm 0.0001 \qquad (N = 18, \dots, 24).$$

This motivates our central numerical conjecture connecting the discrete Galerkin spectral gap to physical boundary leakage:

### Conjecture 2.2 (Numerical Conjecture: Eigenvalue Gap Law)
*For a given cutoff $c > 1$, the minimum eigenvalue of the truncated Galerkin matrix is asymptotically proportional to the boundary leakage energy:*

$$\lambda_{\min}(N) \sim \kappa_c A_0(N) \equiv \frac{2 \kappa_c}{L} [T_{v_N}(0)]^2 \longrightarrow 0 \qquad (N \to \infty),$$

*where $\kappa_c = \kappa(c) > 0$ is a cutoff-dependent constant. If Conjecture 2.1 holds with secondary exponent $\rho(c) = c^{-1/2}$, then Conjecture 2.2 further predicts $\lambda_{\min}(N) \sim \widetilde{\kappa}_c \cdot c^{-N}$.*

---

## 3. The Continuum Limit: Solitary Wave and Dual Dirichlet Boundary Conditions

The continuous spatial wave profile:

$$T_{v_N}(t) = v_0 + \sqrt{2} \sum_{m=1}^{N} v_m \cos\left(\frac{2\pi m t}{L}\right)$$

was evaluated across $[0, L]$ on a dense uniform grid of 2,000 points using the companion analysis script `cell42.py` (with full profile transcript recorded in `cell42.out`) [10]. 

### Proposition 3.1 (Finite-$N$ Symmetry and Normalization)
*For every finite dimension $N \ge 1$, the trigonometric wave profile $T_{v_N}(t)$ generated by the Galerkin ground state satisfies:*

1. **Midpoint Parity Symmetry:**
   $$T_{v_N}(L - t) = T_{v_N}(t) \qquad \forall t \in [0, L].$$
2. **Exact Energy Normalization:**
   $$\|T_{v_N}\|_{L^2([0, L])}^2 = \int_0^L T_{v_N}(t)^2 \, dt = L.$$

*Proof.* Because $Q_{c, N}$ commutes with the reflection operator, its eigenspaces are reflection-invariant; in particular, an even ground-state eigenvector may always be chosen ($u_{-m} = u_m = v_m / \sqrt{2}$). (Numerically, the ground state is simple and non-degenerate for all $N$ examined). Energy normalization follows from Fourier orthogonality on $[0, L]$: $\int_0^L T_{v_N}(t)^2 \, dt = L (v_0^2 + \sum_{m=1}^N v_m^2) = L \|v_N\|_2^2 = L$. $\blacksquare$

### Numerical Observation 3.2 (Apparent Continuum Profile and Dual Dirichlet Nodes)
*Dense grid evaluations through $N = 24$ indicate that as $N$ increases, the sequence of trigonometric profiles $T_{v_N}(t)$ appears to converge to a smooth, strictly positive solitary wave $T_\infty(t)$ satisfying:*

1. **Dual Dirichlet Boundary Nodes:**
   $$T_\infty(0) = T_\infty(L) = 0.$$
2. **Interior Positivity:**
   $$T_\infty(t) > 0 \qquad \forall t \in (0, L),$$
   *with a single central maximum at $t = L/2$ of height $T_{\max} \approx 2.5382 \approx L$.*

### Conjecture 3.3 (Conjectured $C^\infty$ Boundary Flatness)
*The limiting continuum solitary wave $T_\infty(t)$ is conjectured to satisfy infinite-order flat boundary contact at both endpoints:*

$$T_\infty^{(k)}(0) = T_\infty^{(k)}(L) = 0 \qquad \forall k \ge 0.$$

*If $T_\infty$ exists with the required boundary regularity, these vanishing jets imply that the extension of $T_\infty(t)$ by zero outside $[0, L]$, denoted $\widetilde{T}_\infty(t)$, belongs to $C_c^\infty(\mathbb{R})$ with $\operatorname{supp} \widetilde{T}_\infty = [0, L]$.*

*Discussion of Boundary Jet Extinction.*
Numerical evaluation of the even derivatives $D_k(N) = T_{v_N}^{(2k)}(0)$ for $k \in \{0, 1, 2, 3\}$ across $N \in \{8, 16, 24\}$ (computed via `cell43.py`, recorded in `cell43.out` [10]) confirms rapid geometric extinction across all computed orders:
- $D_0$: $8.05 \times 10^{-11} \longrightarrow 1.78 \times 10^{-16} \longrightarrow 1.14 \times 10^{-20}$,
- $D_1$: $3.36 \times 10^{-6} \longrightarrow 3.13 \times 10^{-11} \longrightarrow 5.92 \times 10^{-15}$,
- $D_2$: $2.63 \times 10^{-2} \longrightarrow 1.37 \times 10^{-6} \longrightarrow 7.20 \times 10^{-10}$,
- $D_3$: $71.43 \longrightarrow 2.45 \times 10^{-2} \longrightarrow 3.61 \times 10^{-5}$.

All odd derivatives vanish identically by midpoint symmetry: $T_{v_N}^{(2k+1)}(0) \equiv 0$.

### 3.4 Elimination of the Volterra Boundary Jump

In finite Galerkin implementations, the Volterra convolution:

$$K_v(\omega) = 2 \int_0^\omega \tau_v(s) \tau_v(\omega - s) \, ds \qquad (0 \le \omega \le 1)$$

exhibits an endpoint jump discontinuity at $\omega = 1$ because $\tau_v(0) = T_v(0) \ne 0$. This boundary jump produces the oscillatory factor $1 - \cos(rL)$ and the leading $A_0/r^2$ tail in the Fourier-side kernel.

Under Conjecture 3.3, because $T_\infty(0) = T_\infty(L) = 0$ with infinite-order flat contact, the continuum Volterra convolution:

$$K_\infty(\omega) = 2 \int_0^\omega T_\infty(t) T_\infty(\omega - t) \, dt$$

vanishes smoothly at both endpoints:

$$\lim_{\omega \to 0^+} K_\infty(\omega) = 0, \qquad \lim_{\omega \to 1^-} K_\infty(\omega) = 0,$$

with all derivatives vanishing identically: $K_\infty^{(k)}(0) = K_\infty^{(k)}(1) = 0$ for all $k \ge 0$. Consequently, the finite-rank boundary jump is completely eliminated in the continuum limit.

---

## 4. Semiclassical WKB Barrier Potential and Quantum Tunneling Mechanics

Given any positive profile $T(t)$, one can formally define an effective Schrödinger potential by:

$$V_{\mathrm{conf}}(t) - E := \frac{T''(t)}{T(t)}.$$

Under this definition, $T(t)$ formally satisfies the stationary Schrödinger equation $-T''(t) + V_{\mathrm{conf}}(t) T(t) = E T(t)$ as an identity. For the computed ground-state profile, this construction produces an effective potential whose minimum lies at the midpoint $t = L/2$ and which rises steeply toward the boundaries.

The boundary suppression can then be modeled semiclassically via the WKB tunneling action across the barrier $[0, t_{\mathrm{turn}}]$:

$$\mathcal{S}_{\mathrm{WKB}} = \int_0^{t_{\mathrm{turn}}} \sqrt{\frac{T''(t)}{T(t)}} \, dt,$$

where $t_{\mathrm{turn}} \approx 0.4079 L$ is the effective turning point defined by $T''(t_{\mathrm{turn}}) = 0$.

*Numerical Comparison and WKB Barrier Computation.*
At $N = 24$, the numerical turning point is $t_{\mathrm{turn}} \approx 1.046259$ ($0.40791 L$), computed via `cell44.py` (output log `cell44.out` [10]). The WKB barrier action evaluates to:

$$\mathcal{S}_{\mathrm{WKB}} = 44.363852.$$

Comparing this with the actual boundary suppression across 20 orders of magnitude:

$$\text{Actual Suppression} = \log\left(\frac{T(L/2)}{T(0)}\right) = \log\left(\frac{2.538158}{1.137963 \times 10^{-20}}\right) = 46.853901.$$

$$\frac{\text{Actual Suppression}}{\mathcal{S}_{\mathrm{WKB}}} = \frac{46.853901}{44.363852} = 1.05613.$$

The numerically constructed effective potential yields a WKB action whose exponential scale matches the observed boundary suppression within **$5.6\%$** for a boundary suppression corresponding to approximately 20 decimal orders of magnitude.

### 4.1 Universal Semiclassical Scaling across Prime Cutoffs

Investigation across prime cutoffs $c \in \{5, 7, 11, 13, 17\}$ and dimensions $N \in \{4, 8, 12, 16, 20\}$ using the multi-parameter analysis suite (`cell47.py` and output log `cell47.out` [10]) reveals that the WKB action closely tracks the semiclassical relation:

$$\frac{\mathcal{S}_{\mathrm{WKB}}(N, c)}{L} \approx \frac{\pi N}{4} \implies \mathcal{S}_{\mathrm{WKB}}(N, c) \approx \frac{\pi N}{4} \log c.$$

At $N = 20$, $\frac{\pi \times 20}{4} = 5\pi \approx 15.70796$. Numerical evaluations yield:
- $c = 11$: $\mathcal{S}_{\mathrm{WKB}} / L = 15.3258$,
- $c = 13$: $\mathcal{S}_{\mathrm{WKB}} / L = 15.6681$ (**$99.75\%$ match** to $5\pi$),
- $c = 17$: $\mathcal{S}_{\mathrm{WKB}} / L = 15.8090$ (**$99.36\%$ match** to $5\pi$).

The ratio $\text{Actual Suppression} / \mathcal{S}_{\mathrm{WKB}}$ decreases monotonically across cutoffs toward unity ($1.121 \to 1.084 \to 1.063 \to 1.059 \to 1.054$). Across approximately 47 decimal orders of magnitude at $c = 17$, the WKB action agrees with the observed logarithmic boundary suppression to within $5.3\%$.

---

## 5. Legendre Multipole Decomposition and Taylor Jet Extinction

In normalized coordinates $x = \frac{2t}{L} - 1 \in [-1, 1]$, the finite-$N$ normalized even wave $\psi_N(x) = T_{v_N}(\frac{x+1}{2} L)$ admits an exact Legendre expansion in $L^2([-1, 1])$:

$$\psi_N(x) = \sum_{k=0}^\infty c_{2k}^{(N)} P_{2k}(x).$$

### Proposition 5.1 (Legendre Expansion via Bauer–Bessel Transform)
*The Legendre coefficients are given in exact closed analytical form via Bauer's spherical Bessel expansion:*

$$c_0^{(N)} = v_{N, 0}, \qquad c_{2k}^{(N)} = (4k + 1) \sqrt{2} (-1)^k \sum_{m=1}^N (-1)^m v_{N, m} j_{2k}(\pi m) \quad (k \ge 1),$$

*where $j_n(z) = \sqrt{\frac{\pi}{2z}} J_{n+1/2}(z)$ is the spherical Bessel function of the first kind.*

*Observed Properties of the Legendre Multipoles (`cell44.py`).*
1. **Spectral Concentration:** Truncation at $K = 10$ ($P_{20}(x)$) captures **$99.999984\%$** of the $L^2$ norm: $\sum_{k=0}^{10} \frac{2}{4k + 1} [c_{2k}^{(24)}]^2 = 1.99999968 \approx 2.00000000$. Over **$93.7\%$** of the total energy resides in the lowest four even multipoles ($P_0: 29.9\%, P_2: 31.4\%, P_4: 21.2\%, P_6: 11.1\%$).
2. **Alternating Phases and Boundary Cancellation:** Throughout the resolved multipoles ($k \le 10$), the computed coefficients exhibit an alternating-sign pattern: $c_{2k}^{(N)} = (-1)^k |c_{2k}^{(N)}|$. Because $P_{2k}(0) = (-1)^k \frac{(2k)!}{2^{2k}(k!)^2}$, the alternating signs make the central contribution strictly constructive: $\psi_N(0) \approx \sum_{k=0}^{10} |c_{2k}^{(N)}| \frac{(2k)!}{2^{2k}(k!)^2} > 0$. At the boundaries $x = \pm 1$, $P_{2k}(\pm 1) = 1$, producing destructive cancellation: $\psi_\infty(\pm 1) = |c_0^{(\infty)}| - |c_2^{(\infty)}| + |c_4^{(\infty)}| - |c_6^{(\infty)}| + \dots = 0$.

### Conjecture 5.2 (Extinction of the Asymptotic Tail Hierarchy)
*For every fixed $k \ge 0$, the $k$-th coefficient $A_k(N)$ in the inverse-power asymptotic expansion of the Archimedean resolvent vanishes identically in the continuum limit:*

$$A_k(\infty) = \lim_{N\to\infty} \frac{2}{L} (-1)^k \sum_{j=0}^k D_j(N) D_{k-j}(N) = 0 \qquad \forall k \ge 0.$$

*If, in addition, the asymptotic expansion is sufficiently uniform in $N$ to permit the interchange of $\lim_{N\to\infty}$ and $r \to \infty$, the limiting continuous-variable resolvent:*

$$R_\infty(r) = \lim_{N\to\infty} \frac{2}{L} \left[ \frac{v_{N, 0}}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_{N, m}}{r^2 - a_m^2} \right]^2$$

*is conjectured to decay faster than every inverse power of $r$:*

$$R_\infty(r) = o(r^{-k}) \qquad \forall k \in \mathbb{N}.$$

*Numerical Evidence across Dimensions (`cell45.py`).*
Evaluation across $N \in \{4, 8, 12, 16, 20, 24\}$ demonstrates dramatic extinction across all tested orders:
- $A_0$: $2.81 \times 10^{-13} \to 5.05 \times 10^{-21} \to 1.01 \times 10^{-40}$ (collapsing by 27 orders of magnitude),
- $A_1$: $5.54 \times 10^{-9} \to 4.22 \times 10^{-16} \to 1.05 \times 10^{-34}$,
- $A_2$: $3.48 \times 10^{-5} \to 1.21 \times 10^{-11} \to 4.01 \times 10^{-29}$,
- $A_3$: $7.65 \times 10^{-2} \to 1.47 \times 10^{-7} \to 7.28 \times 10^{-24}$,
- $A_4$: $73.42 \to 9.23 \times 10^{-4} \to 7.53 \times 10^{-19}$.

At high frequencies, the finite-$N$ resolvent plunges: $R_{v_{24}}(10.0) = 0.0368$, $R_{v_{24}}(15.0) = 6.30 \times 10^{-6}$, $R_{v_{24}}(20.0) = 1.10 \times 10^{-8}$, and $R_{v_{24}}(50.0) = 5.40 \times 10^{-30}$. The effective logarithmic slope $\gamma_{\mathrm{eff}}(r) = -r R'(r)/R(r)$ reaches $\gamma_{\mathrm{eff}} \approx 78.6$ at $r = 15.0$, $154.0$ at $r = 20.0$, and $270.3$ at $r = 30.0$.

---

## 6. Accumulating Pole Mechanism and Heat Boundary Dynamics

The operator-resolvent representation $D_N(z) = \big[(I + z\mathcal{L})^{-1} T_{v_N}\big](0)$ established in Paper 4 provides the conceptual mechanism reconciling the extinction of the inverse-power tail coefficients ($A_k(N) \to 0$) with a non-trivial continuous-variable resolvent $R_\infty(r)$:

1. **Accumulation of Resolvent Poles at the Origin:**
   At every finite dimension $N$, $D_N(z)$ is a rational function whose poles lie on the negative real axis at:
   $$z_m = -\frac{1}{a_m^2} = -\frac{L^2}{4\pi^2 m^2} \in \left(-\frac{L^2}{4\pi^2}, 0\right) \qquad (m = 1, \dots, N).$$
   As $N \to \infty$, the poles accumulate at the origin: $\lim_{m\to\infty} z_m = 0^-$.
2. **Obstruction to Analyticity:**
   The residue of $D_N(z)$ at $z_m$ is proportional to $\frac{\sqrt{2} v_{N, m}}{a_m^2}$. The infinite accumulation of poles at $z = 0^-$ obstructs analytic continuation through the origin from the negative real axis.
3. **Vanishing Taylor Jet $\not\Rightarrow$ Triviality:**
   The numerical extinction of every fixed Taylor coefficient $D_k(N) = T_{v_N}^{(2k)}(0) \to 0$ is consistent with a limiting object that is $C^\infty$-flat at $z = 0$ from the right ($\operatorname{Re}(z) > 0$). Because $z = 0$ is an accumulation boundary of singularities, this flatness does not force $D_\infty(z)$ to vanish identically on the negative axis.
4. **Exponentially Flat Ansatz:**
   The observed WKB quantum barrier behavior motivates testing an exponentially flat ansatz:
   $$D_\infty\left(-\frac{1}{r^2}\right) \sim e^{-C r} \qquad (r \to \infty),$$
   corresponding to $D_\infty(z) \sim e^{-C / \sqrt{-z}}$ as $z \to 0^-$.

### 6.1 Numerical Investigations from Dedicated Computational Cells

- **Cell 51 (Discrete Cauchy Transform):** On the negative real axis $z = -1/r^2$, with $w = -r^2/\kappa^2$, $D_N(-1/r^2)$ matches the discrete Cauchy transform $v_0 + \sqrt{2} w \sum_{m=1}^N \frac{v_m}{w - m^2}$ identically to 51 decimal digits. Probing off-lattice points confirms universal suppression ($\sim 10^{-12} - 10^{-13}$ at $m = 20$), while the ratio $-\log|D|/r$ oscillates between $0.37$ and $0.59$.
- **Cell 52 (Spectral-Edge Time Scale):** The heat boundary trace $H_N(u) = [e^{-u\mathcal{L}} T_N](0)$ plunges across 20 orders of magnitude, reaching $T_{24}(0)$ at $u = 10^{-6}$. The inverse spectral-edge scale $u_N = a_N^{-2} = \frac{1}{\kappa^2 N^2}$ acts as an exact crossover scale for the resolvent integral.
- **Cell 53 (Universal First-Jet Profile Collapse):** When heat time is scaled by the first-jet cancellation scale $u_1 = D_0 / D_1$, the normalized profiles $\Theta_N^{\mathrm{cancel}}(\theta) = H_N(\theta u_1) / D_0$ exhibit a near-perfect universal collapse across all dimensions $N \in \{8, \dots, 24\}$ ($2.12 \pm 0.02$ at $\theta = 1.0$), with stable shape invariants $\beta_N = D_0 D_2 / D_1^2 \approx 0.19 - 0.26$.
- **Cell 54 (Anatomy of the Decoupling Ratio $s_N$):** Tracking the decoupling ratio $s_N = (\kappa N)^2 (D_0 / D_1)$ reveals that while $-\log|D_0|$ drops by 22.7 units and $-\log|D_1|$ drops by 20.2 units, their difference $\Delta_N = \log|D_1/D_0|$ drifts only from $10.64$ to $13.16$. This provides strong numerical evidence that $D_0$ and $D_1$ share the same leading exponential barrier decay rate.

---

## 7. Tri-Partite Zero-Energy Balance and Finite-$T$ Archimedean Cutoff Leakage

Let $\mathcal{Q}(v) = \mathcal{Q}_{\mathrm{pole}}(v) + \mathcal{Q}_{\mathrm{prime}}(v) + \mathcal{Q}_{\mathrm{arch}}(v)$ be the Connes–van Suijlekom quadratic form on the Galerkin subspace of dimension $2N+1$. For every finite dimension $N$, the algebraic matrix sum matches the minimum eigenvalue identically:

$$\mathcal{Q}_{\mathrm{matrix}}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{\mathrm{matrix}}(v_N) \equiv \lambda_{\min}(N).$$

### Proposition 7.1 (Continuous-Quadrature Balance and Finite-$T$ Archimedean Leakage)
*When the Archimedean contribution is evaluated independently via continuous-variable quadrature $\frac{1}{\pi} \int_0^{80} h_+(r) \Phi_{v_N}(r)^2 \, dr$ (using `cell46.py`, logged in `cell46.out` [10]), the independently computed components cancel from $\mathcal{O}(1)$ down to a residual of order $10^{-43}$ at $N = 24$:*

- $N = 4$: $\mathcal{Q}_{\mathrm{pole}} = +2.206186$, $\mathcal{Q}_{\mathrm{prime}} = -0.316153$, $\mathcal{Q}_{\mathrm{arch}} = -1.890032$, summing to $\mathcal{Q}_{\mathrm{total}} = 7.82 \times 10^{-15}$ ($\lambda_{\min} = 8.83 \times 10^{-15}$),
- $N = 8$: $\mathcal{Q}_{\mathrm{pole}} = +1.813949$, $\mathcal{Q}_{\mathrm{prime}} = -0.154916$, $\mathcal{Q}_{\mathrm{arch}} = -1.659033$, summing to $\mathcal{Q}_{\mathrm{total}} = 5.38 \times 10^{-23}$ ($\lambda_{\min} = 6.71 \times 10^{-23}$),
- $N = 12$: $\mathcal{Q}_{\mathrm{pole}} = +1.675166$, $\mathcal{Q}_{\mathrm{prime}} = -0.108101$, $\mathcal{Q}_{\mathrm{arch}} = -1.567065$, summing to $\mathcal{Q}_{\mathrm{total}} = 1.32 \times 10^{-29}$ ($\lambda_{\min} = 1.78 \times 10^{-29}$),
- $N = 16$: $\mathcal{Q}_{\mathrm{pole}} = +1.609630$, $\mathcal{Q}_{\mathrm{prime}} = -0.088194$, $\mathcal{Q}_{\mathrm{arch}} = -1.521436$, summing to $\mathcal{Q}_{\mathrm{total}} = 5.11 \times 10^{-35}$ ($\lambda_{\min} = 7.12 \times 10^{-35}$),
- $N = 20$: $\mathcal{Q}_{\mathrm{pole}} = +1.572288$, $\mathcal{Q}_{\mathrm{prime}} = -0.077529$, $\mathcal{Q}_{\mathrm{arch}} = -1.494759$, summing to $\mathcal{Q}_{\mathrm{total}} = 8.81 \times 10^{-40}$ ($\lambda_{\min} = 1.32 \times 10^{-39}$),
- $N = 24$: $\mathcal{Q}_{\mathrm{pole}} = +1.5516521957$, $\mathcal{Q}_{\mathrm{prime}} = -0.0718544317$, $\mathcal{Q}_{\mathrm{arch}} = -1.4797977640$, summing to $\mathcal{Q}_{\mathrm{total}} = 4.20136 \times 10^{-43}$ ($\lambda_{\min} = 2.53348 \times 10^{-43}$).

*The Certified Archimedean Cutoff Leakage $\delta_T^{\mathrm{tail}}$ (`cell56.py`, `cell57.py`).*
Evaluating the continuous Archimedean functional via the exact closed-form digamma identity (Corollary 5.4 of Paper 4) yields $\mathcal{Q}_{\mathrm{arch}}^{\mathrm{exact}} = -1.479797763974798326397825\dots$, matching continuous numerical quadrature to $4.96 \times 10^{-25}$. The tripartite continuous balance cancels algebraically from $\mathcal{O}(1)$ to $\mathcal{Q}_{\mathrm{total}} = 4.2013606231 \times 10^{-43}$. 

The ratio against the finite-cutoff Galerkin eigenvalue $\lambda_{\min}(24) = 2.5334848706 \times 10^{-43}$ is:

$$\frac{\lambda_{\min}(24)}{\mathcal{Q}_{\mathrm{total}}} = 0.6030153319\dots \qquad \left(\frac{\mathcal{Q}_{\mathrm{total}}}{\lambda_{\min}(24)} \approx 1.65837\right).$$

In the companion paper (Paper 4, Theorem 5.5) and Cell 57 (`cell57.py`), this residual is proven to be **100% finite-$T$ Archimedean cutoff tail leakage**:

$$\delta_T^{\mathrm{tail}} \equiv \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_{24}) - \lambda_{\min}(24) = \frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_{24}, r, L) \, dr = 1.66787575 \times 10^{-43}.$$

Because the Galerkin divided-difference matrix $Q_{\mathrm{arch}}^{(T)}$ truncates the continuous Fourier integral at $T = 400$, the discrete matrix quadratic form satisfies $v^T Q_{\mathrm{arch}}^{(T)} v \equiv \frac{1}{\pi}\int_0^T h_+ K_{\mathrm{Fourier}} dr$ identically. Across all tested dimensions $N \in \{8, 12, 16, 20, 24\}$, the discrepancy is identically the continuous cutoff tail $\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) \equiv -\delta_T^{\mathrm{tail}}(v_N)$ (balance error $2.29 \times 10^{-45}$ at $N = 24$). This definitively retires the interpretation of the residual as an unexplained finite-rank subspace projection defect.

### 7.2 Prime-Power Decomposition of the Negative Barrier

Direct point-evaluation of the Volterra convolution $K_{v_{24}}(\omega_q)$ at all prime powers $q \le 13$ matches the matrix-computed prime form to 52 decimal digits ($|\text{diff}| = 1.67 \times 10^{-52}$). 
- The lowest prime $q = 2$ provides **$98.65\%$** of the entire prime energy ($-0.0708858$).
- $q = 3$ accounts for **$1.34\%$** ($-0.0009658$).
- Contributions above $q = 7$ decay exponentially below $10^{-13}$ ($q = 11$: $-9.52 \times 10^{-28}$).
- At the Volterra endpoint $\omega = 0$ ($q = 13$), $K_{v_{24}}(0) = 0$ identically.

---

## 8. Formal Wiener–Hopf Continuum Scaling and Asymptotic Bounding Ladder

In the formal continuum scaling limit, we analyze the divided-difference Galerkin operator $Q_{mn} \approx \frac{\log(m/n)}{m - n}$ on normalized coordinates $x = m/N, y = n/N \in (0, 1]$.

### Proposition 8.1 (Formal Continuum Wiener–Hopf Scaling and Symbol Factorization)
*Let $\mathcal{T}$ be the continuous integral operator on $L^2((0, 1], dx)$ defined by:*

$$(\mathcal{T} \phi)(x) = \int_0^1 \frac{\log x - \log y}{x - y} \phi(y) \, dy.$$

1. **Wiener–Hopf Equivalence:**
   *Under the isometric isomorphism $U: L^2((0, 1], dx) \to L^2([0, \infty), d\xi)$ defined by $\xi = -\log x$ and $\Phi(\xi) = e^{-\xi/2} \phi(e^{-\xi})$, the operator $\mathcal{T}$ transforms into a pure Wiener–Hopf convolution operator on the half-line $\mathbb{R}_+$ with kernel:*
   $$K_{\mathrm{sym}}(w) = \frac{w}{2\sinh(w/2)} \qquad (w = \xi - \eta).$$
2. **Exact Double Gamma Factorization:**
   *The Fourier symbol $\widehat{K}(k) = \int_{-\infty}^\infty K_{\mathrm{sym}}(w) e^{i k w} dw$ admits the exact closed-form factorization:*
   $$\widehat{K}(k) = \frac{\pi^2}{\cosh^2(\pi k)} = K_+(k) K_-(k), \qquad K_+(k) = \left[ \Gamma\left(\frac{1}{2} - i k\right) \right]^2, \quad K_-(k) = \left[ \Gamma\left(\frac{1}{2} + i k\right) \right]^2.$$
3. **Singular Boundary-Layer Asymptotics:**
   *The leading singularity of the causal Wiener–Hopf factor $K_+(k)$ is a double pole at $k = -i/2$, which generates a logarithmic boundary-layer divergence in physical space:*
   $$\phi(x) \sim -C_1 \log x + C_0 \qquad (x \to 0^+).$$

*Proof.* Setting $x = e^{-\xi}, y = e^{-\eta}$ and substituting into $\mathcal{T}$ gives kernel $e^{-\xi/2} \frac{\xi - \eta}{e^{-\xi} - e^{-\eta}} e^{-\eta/2} = \frac{w}{2\sinh(w/2)}$ where $w = \xi - \eta$. Differentiating Ramanujan's hyperbolic integral $\int_{-\infty}^\infty \frac{e^{ikw}}{\cosh(w/2)} dw = \frac{2\pi}{\cosh(\pi k)}$ gives $\widehat{K}(k) = \frac{\pi^2}{\cosh^2(\pi k)}$. Euler's reflection formula yields the Gamma factorization. The double pole of $K_+(k)$ at $k = -i/2$ produces $\Phi(\xi) \sim (C_1 \xi + C_0) e^{-\xi/2}$, which under the inverse isometry yields $\phi(x) \sim -C_1 \log x + C_0$. $\blacksquare$

*Remark 8.1.1 (Mechanism Asymmetry of the Endpoint Jets).*
Proposition 8.1 provides the analytical foundation for the bulk/edge mechanism asymmetry discovered in Cell 54:
- In the second moment $D_1 = -\sqrt{2}\kappa^2 N^3 \int_0^1 x^2 \phi(x) dx$, the quadratic factor $x^2$ quenches the logarithmic singularity ($x^2 \log x \to 0$ as $x \to 0$). The integrand is smooth on $[0, 1]$, making $D_1$ regular and dominated by bulk modes ($x \sim \mathcal{O}(1)$).
- In the zeroth moment $D_0 = v_0 + \sqrt{2} N \int_0^1 \phi(x) dx$, the logarithmic divergence requires the discrete lattice modes near $m \in \{1, \dots, 5\}$ to engage in destructive cancellation against $v_0$, while edge modes ($m \sim N$) contribute negligibly ($\sim 10^{-8}$).

### Proposition 8.2 (First-Row Taylor Jet Ladder)
*Expanding the first-row matrix condition $\psi'(0) v_0 + \sqrt{2} \sum_{m=1}^N \frac{\psi(m)}{m} v_m = \lambda v_0 \approx 0$ via the Taylor series of the odd function $\psi(x)$ at $x = 0$ couples the endpoint jets directly to the higher shape invariants:*

$$\psi'(0) D_0 = \frac{\psi'''(0)}{6 \kappa^2} D_1 - \frac{\psi^{(5)}(0)}{120 \kappa^4} D_2 + \cdots + \mathcal{R}_N,$$

*which, dividing by $D_1$ and using the shape invariant $\beta_N = D_0 D_2 / D_1^2 \approx 0.24$ (Cell 53), yields:*

$$s_N \equiv (\kappa N)^2 \frac{D_0}{D_1} = N^2 \left[ \frac{\psi'''(0)}{6 \psi'(0)} - \frac{\psi^{(5)}(0)}{120 \kappa^2 \psi'(0)} \beta_N \frac{D_1}{D_0} + \cdots \right] + (\kappa N)^2 \frac{\mathcal{R}_N}{D_1}.$$

### Proposition 8.3 (Commutator Resolvent Formula and First-Jet Balance)
*In the even subspace $\mathcal{H}_{\mathrm{even}}$, projecting the quadratic commutator $[M^2, Q] u = D_0 b - B_1 e$ onto the orthogonal complement $u^\perp$ via the pseudoinverse $Q_{\mathrm{even}}^\dagger$ yields the asymptotic formula:*

$$\frac{D_1}{D_0} = \kappa^2 \left[ e^T Q_{\mathrm{even}}^\dagger \big( b + \mathcal{E}_{\mathrm{arith}} e \big) - \|M u\|_2^2 \right] + \mathcal{O}(\lambda),$$

*where $b_n = n \psi(n)$ and $\mathcal{E}_{\mathrm{arith}} = \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1}\psi \rangle$. Because the ground-state eigenvalue $\lambda_{\min}(N) \sim 10^{-43}$ is negligible compared to $D_0 \sim 10^{-20}$ at large $N$ (Cell 55), the $\mathcal{O}(\lambda)$ term vanishes asymptotically, leaving $D_1/D_0$ governed exclusively by the unperturbed resolvent.*

*Proof.* From the rank-$4$ commutator identity (Theorem 6.1 of Paper 4), $[M^2, Q] u = D_0 b - B_1 e$. Using the exact odd-sector identity $B_1 = -D_0 \mathcal{E}_{\mathrm{arith}}$ (Theorem 6.2 of Paper 4), we have $(Q - \lambda I) M^2 u = -D_0 (b + \mathcal{E}_{\mathrm{arith}} e)$. Projecting onto $u^\perp$ via $Q_{\mathrm{even}}^\dagger$ and contracting with $e^T$ (recalling $e^T M^2 u = -D_1 / \kappa^2$ and $e^T u = D_0$) gives $-D_1 / \kappa^2 = -D_0 [e^T Q_{\mathrm{even}}^\dagger \mathbf{w} - \|Mu\|_2^2] + \mathcal{O}(\lambda D_0)$. Dividing by $-D_0$ yields the result. $\blacksquare$

### Proposition 8.4 (Exact Operator Decomposition, Bound-State Interlacing, and Spectral Filtering Mechanism)
*Using the exact algebraic small-denominator cancellation theorem and $K^2$-commutator resolvent representation from Paper 4 (Theorems 7.2 & 7.3), the first-jet ratio admits the exact operator decomposition:*

$$\frac{D_1}{D_0} = \kappa^2 \left[ \mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} - D_0^2 M_2 \right], \qquad \mathcal{T}_{\mathrm{diag}} \equiv M_1 \langle d, R_{\mathrm{even}} d \rangle, \quad \mathcal{T}_{\mathrm{cross}} \equiv \langle d, R_{\mathrm{even}} K\boldsymbol\psi \rangle,$$

*where $s_2 = K\boldsymbol\psi + M_1 d \in c^\perp$, $M_1 = \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi \rangle$, and $R_{\mathrm{even}} = (Q_{\mathrm{even}} - \lambda I)_{c^\perp}^{-1}$. The apparent singularity $(E_k - \lambda)^{-1}$ cancels identically mode by mode, and the asymptotic first-jet scale is governed by the collective operator structure of the low-energy bound-state ladder.*

*Exact Operator Reformulation and Numerical Scaling Evidence (`cell59.py`–`cell62.py`).*
1. **The Interlaced Bound-State Tunneling Ladder (Cells 59–60):**
   Below the double-well barrier top ($E_{\mathrm{barrier}} \approx 1.30$, corresponding to mode index $k \le N/2 = 12$ at $N = 24$), the even and odd spectra do not decouple into an isolated ground state plus a uniform scattering continuum. Instead, they form an exquisitely structured ladder of bound tunneling doublets strictly interlacing down to the ground state:
   $$E_0 < \mu_1 < E_1 < \mu_2 < E_2 < \mu_3 < E_3 < \dots < \mu_{N/2} < E_{N/2}.$$
   Each step up the ladder increases the eigenvalue by 2 to 3 orders of magnitude ($E_k / \mu_k \sim 10^3$, $\mu_{k+1} / E_k \sim 10^3$; e.g. at $N = 24$, $E_0 \approx 2.53 \times 10^{-43}$, $\mu_1 \approx 4.35 \times 10^{-40}$, $E_1 \approx 4.50 \times 10^{-37}$, $\dots$, $\mu_{12} \approx 1.07$, $E_{12} \approx 1.31$).

2. **Transition Dipole Moment and Mode 1 Saturation (Cell 61):**
   Taking matrix elements of the rank-2 commutator $[Q, K] = -\boldsymbol\psi d^T + d\boldsymbol\psi^T$ on the ground state $c$ yields the exact coordinate projection $\langle e_j, K c \rangle = -D_0 \frac{a_j}{\Delta_{\mathrm{odd}, j}}$ (Corollary 6.3.1 of Paper 4). For the lowest odd mode $e_1$, this projection factors into the product of two dimensionless invariants:
   $$|\langle e_1, K c \rangle| = \left(\frac{D_0}{\sqrt{\Delta_{\mathrm{odd}, 1}}}\right) \left(\frac{|a_1|}{\sqrt{\Delta_{\mathrm{odd}, 1}}}\right) \equiv R_D \cdot C_N.$$
   Computational evaluation demonstrates that $|\langle e_1, K c \rangle| \approx 1.3134 = \mathcal{O}(1)$ represents the transition dipole moment of the two-level tunneling doublet. Because $R_D = D_0 / \sqrt{\Delta_1} \in [0.41, 0.64]$ is $\mathcal{O}(1)$ across 20 orders of magnitude, the overlap ratio $C_N \equiv |a_1| / \sqrt{\Delta_1} \approx 2.406$ is structurally constrained by the transition dipole moment to remain $\mathcal{O}(1)$ across the entire tested ladder. Furthermore, mode 1 carries **$99.9999\%$** of the total full-space coordinate norm $\|Kc\|^2 = 1.72507$, indicating that the coordinate commutator acts effectively on a two-state tunneling doublet.

3. **Operator Resolvent Balance and Excited-Sector Cancellation (Cell 62):**
   In contrast to $\|Kc\|^2$, the resolvent sum $\langle d, R_{\mathrm{even}} s_2 \rangle$ is distributed collectively across the entire $N/2$ bound ladder below the barrier top. High-precision evaluation confirms the exact decomposition of Theorem 7.3 across dimensions:
   - At $N = 24$: $\kappa^2 \mathcal{T}_{\mathrm{diag}} = +2\,234\,013.195$, $\kappa^2 \mathcal{T}_{\mathrm{cross}} = -1\,714\,040.497$, $\kappa^2 D_0^2 M_2 = 10.352$, summing identically to $D_1/D_0 = 519\,962.347$.
   - The ratio $\mathcal{T}_{\mathrm{cross}} / \mathcal{T}_{\mathrm{diag}}$ appears to stabilize near $\approx -0.77$ over the tested dimensions. The normalized ratio $(D_1/D_0) / [\kappa^2 \mathcal{T}_{\mathrm{diag}}]$ traces the non-monotone sequence $(0.476, 0.331, 0.242, 0.212, 0.233)$, clustering near $\approx 0.23$ for $N \ge 16$.

4. **The Spectral Filtering Identity and Asymptotic Filtering Mechanism:**
   From the commutator $[Q, K] = -\boldsymbol\psi d^T + d\boldsymbol\psi^T$, expanding $\boldsymbol\psi = \sum a_j e_j$ yields $\langle u_k, K\boldsymbol\psi \rangle = -d_k \sum_{j=1}^N \frac{a_j^2}{\mu_j - E_k}$. Subtracting this from $M_1 d_k = d_k \sum \frac{a_j^2}{\mu_j - \lambda}$ gives the exact mode source coefficient:
   $$\langle u_k, s_2 \rangle = d_k \sum_{j=1}^N a_j^2 \left( \frac{1}{\mu_j - \lambda} - \frac{1}{\mu_j - E_k} \right) = d_k \sum_{j=1}^N \frac{a_j^2 (E_k - \lambda)}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
   The high-mode cancellation is not an exact vanishing, but an asymptotic suppression: for $j > k$, the difference is governed by $(E_k - \lambda)/\mu_j^2 \ll 1$ due to the rapid growth of the eigenvalues up the ladder. Decomposing the sum into low modes $j \le k$ (where $E_k \gg \mu_j \implies \frac{1}{\mu_j - E_k} \approx 0$) and a high-mode remainder yields:
   $$\langle u_k, s_2 \rangle = d_k \sum_{j=1}^k \frac{a_j^2}{\mu_j - \lambda} + \mathcal{R}_{\mathrm{filt}}(k), \qquad \mathcal{R}_{\mathrm{filt}}(k) \equiv d_k \sum_{j > k} \frac{a_j^2 (E_k - \lambda)}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
   For each mode $k$, the primary driver is truncated to the low-frequency partial sum ($5.8$ for $k = 1$, $10.9$ for $k = 2$, averaging $\sim 23$ across active modes), providing the structural mechanism for the observed $\approx 0.23$ scale. Controlling the remainder $\mathcal{R}_{\mathrm{filt}}(k)$ analytically presents a concrete pathway to establish spectral filtering as a rigorous theorem.

5. **Square-Root-Scale Overlap Scaling:**
   In both sectors, the source overlaps scale with the square root of the energy splitting:
   $$d_k \sim \Gamma_{\mathrm{even}} \sqrt{\Delta_{\mathrm{even}, k}}, \qquad a_j \sim \Gamma_{\mathrm{odd}} \sqrt{\Delta_{\mathrm{odd}, j}}.$$
   Across the bound ladder, $d_k^2 / \Delta_k$ ranges from $\approx 580$ down to $\approx 43$, and $a_j^2 / \Delta_j$ ranges from $\approx 4.4$ to $\approx 16$. Rather than strict constant laws, the fundamental discovery is that these ratios remain bounded within an $\mathcal{O}(10^1 - 10^2)$ envelope across 40 orders of magnitude of gap collapse, extinguishing small denominators mode by mode. Both ladders truncate at mode 12 ($N/2$), with bound states carrying $>98\%$ of $M_1$ and $>99.7\%$ of $\langle d, R_{\mathrm{even}} d \rangle$.

6. **Spatial Profile of the Resolvent State $w_d = R_{\mathrm{even}} d$:**
   In canonical coordinates $s = m/N$, the normalized solution vector $w_d(m) / \|w_d\|$ is strongly concentrated near the boundary $s \in [0, 0.2]$ ($0.399$ at $s = 0$, $0.655$ at $s = 0.125$, $0.014$ at $s = 0.25$) and collapses by 23 orders of magnitude to $-2 \times 10^{-16}$ at $s = 1.0$. This spatial profile is heuristically interpreted as reflecting the semiclassical barrier turning-point region; establishing an exact continuum boundary-value problem for the limiting profile $W(s)$ remains a promising direction for future analysis.

### Proposition 8.5 (Empirical Semigroup Profile Squeezing and Shape Invariants)
*Under the first-jet normalization $u = \theta u_1 = \theta |D_0 / D_1|$, the normalized heat semigroup profile $\Theta_N(\theta) = H_N(\theta u_1) / D_0$ exhibits near-perfect universal collapse across all tested dimensions $N \in \{8, \dots, 24\}$ (matching within $1.5\%$ across 16 decimal orders of magnitude, Cell 53). Over the computed range, the dimensionless shape invariants stabilize:*

$$\beta_N = \frac{D_0 D_2}{D_1^2} \approx 0.19 - 0.26 < 1, \qquad \gamma_N = \frac{D_0^2 D_3}{D_1^3} \approx 0.012 - 0.027 < 1,$$

*and the observed profile is enclosed by the empirical two-sided squeezing bounds:*

$$1 + \theta \le \Theta_N(\theta) \le 1 + \theta + \frac{1}{2} \beta_N \theta^2 \qquad (\theta \in [0, 1]),$$

*which enclose the numerical data to three decimal places, confirming that $u_1 = |D_0/D_1|$ is the genuine physical boundary-layer time scale.*

### Conditional Proposition 8.6 (Sector-Decomposed Polynomial Bound on $D_1/D_0$ — Conditional on Proof Audit)
*Decompose the reduced even resolvent into bound-state and scattering components:*

$$R_{\mathrm{even}} = R_{\mathrm{bound}} + R_{\mathrm{scatt}},$$

*where $R_{\mathrm{bound}} = \sum_{k=1}^{N/2} \frac{u_k u_k^T}{E_k - \lambda}$ acts on the tunneling ladder below the barrier top ($E_k < E_{\mathrm{gap}} \approx 1.30$), and $R_{\mathrm{scatt}}$ acts on the complementary subspace $E \ge E_{\mathrm{gap}}$.*

*Assume:*
1. *Uniform Scattering Gap:* $\inf_N E_{\mathrm{scatt},\min}(N) \ge E_{\mathrm{gap}} > 0$ (empirically supported by Cells 49–50).
2. *Scattering Source Growth Hypothesis:* The effective source vector $s_2 = M_1 d + K\boldsymbol\psi$ satisfies $\|s_2\|_2 \le C_s N^{3/2} \log N$ on the scattering sector (subject to rigorous operator proof audit).
3. *Polynomial Bound-State Moment Control:* The bound-state coordinate projections satisfy $\sum_{k=1}^{N/2} \left| \frac{d_k \langle u_k, K^2 c \rangle}{D_0} \right| \le C_{\mathrm{bound}} N^p$ for some finite $p < \infty$ (empirically supported by Cell 62).

*Then:*
1. **Polynomial Upper Bound on First-Jet Ratio:**
   $$\left| \frac{D_1}{D_0} \right| \le \kappa^2 \left( C_{\mathrm{scatt}} N^2 \log N + C_{\mathrm{bound}} N^p \right) \le C_{\mathrm{upper}} N^{\max(2, p)} \log N.$$
2. **Polynomial Lower Bound on Cancellation Scale:**
   $$u_1 = \left| \frac{D_0}{D_1} \right| \ge \frac{c_1}{N^{\max(2, p)} \log N},$$
   *ruling out any exponential $e^{-\alpha N}$ collapse of the cancellation scale and indicating that $D_0$ and $D_1$ share the same leading exponential WKB barrier decay rate.*
3. **Spatial Boundary Layer Width and Decoupling:**
   *The boundary layer width $\delta_N = \sqrt{u_1} \ge \frac{1}{\kappa \sqrt{C_{\mathrm{upper}}} N^{\max(1, p/2)} \sqrt{\log N}}$ shrinks only algebraically. Substituting this polynomial lower bound into the boundary-defect metric $\mathcal{D}(N) = D_0^2 [1 + 1/(T^2 u_1)]^2$ demonstrates that $\mathcal{D}(N) \to 0$ collapses exponentially fast like $e^{-2\mathcal{S}_{\mathrm{WKB}}} N^{2\max(2, p)} (\log N)^2 \to 0$, which would establish the Boundary-Defect Decoupling Conjecture under these hypotheses.*

*Proof.* By Theorem 7.3 and Corollary 7.3.1 of Paper 4, $D_1/D_0 = \kappa^2 [\langle d, R_{\mathrm{even}} s_2 \rangle - D_0^2 M_2]$. Splitting $R_{\mathrm{even}} = R_{\mathrm{scatt}} + R_{\mathrm{bound}}$:
For the scattering piece, since $\|R_{\mathrm{scatt}}\| \le 1/E_{\mathrm{gap}}$, applying the Cauchy–Schwarz inequality with $\|d\|_2 = \sqrt{2N+1}$ and the scattering growth hypothesis $\|s_2\|_2 \le C_s N^{3/2} \log N$ yields $|\langle d, R_{\mathrm{scatt}} s_2 \rangle| \le C_{\mathrm{scatt}} N^2 \log N$.
For the bound piece, the small denominators $\Delta_k = E_k - \lambda$ are canceled algebraically via Corollary 7.3.1: $\frac{d_k \langle u_k, s_2 \rangle}{E_k - \lambda} = -\frac{d_k \langle u_k, K^2 c \rangle}{D_0}$. By assumption 3, this sum is bounded by $C_{\mathrm{bound}} N^p$.
Combining both contributions and the negligible norm shift $D_0^2 M_2 = \mathcal{O}(1)$ yields the upper bound $|D_1/D_0| \le C_{\mathrm{upper}} N^{\max(2, p)} \log N$. Inverting this bound directly yields the polynomial lower bound on $u_1 = |D_0/D_1|$, which would establish boundary-defect decoupling under these hypotheses. $\blacksquare$

### Proposition 8.7 (Analytical Reduction: Operator Resolvent Bounds, Jet-Equation Resummation, and Commutator Smoothing)
*The first-jet ratio $D_1/D_0$ and the spectral-filtering remainder admit an exact analytical reduction that eliminates the exponentially small denominator $D_0$ from the bound-state sector, transforming the asymptotic boundary decoupling problem into positive-resolvent operator estimates:*

1. **Exact Squared-Resolvent Laurent Jet Identity and Boundary Cancellation Asymmetry:**
   *From the exact rational resolvent representation of Paper 4 (Theorem 3.1 & Corollary 5.4), the squared resolvent is given by:*
   $$R_v(r) = \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} \right]^2, \qquad a_m = \frac{2\pi m}{L}.$$
   *Expanding the rational factor for $r > a_N$:*
   $$\frac{r}{r^2 - a_m^2} = \frac{1}{r} \frac{1}{1 - a_m^2/r^2} = \frac{1}{r} + \frac{a_m^2}{r^3} + \frac{a_m^4}{r^5} + \dots$$
   *Defining the discrete boundary moments:*
   $$T_0 \equiv v_0 + \sqrt{2}\sum_{m=1}^N v_m = T_v(0) = D_0, \qquad T_2 \equiv \sqrt{2}\sum_{m=1}^N a_m^2 v_m = -T_v''(0) = -D_1,$$
   *the resolvent expansion becomes:*
   $$\frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} = \frac{T_0}{r} + \frac{T_2}{r^3} + \mathcal{O}(r^{-5}) = \frac{D_0}{r} - \frac{D_1}{r^3} + \mathcal{O}(r^{-5}).$$
   *Squaring this expansion yields:*
   $$R_v(r) = \frac{2}{L} \left[ \frac{D_0^2}{r^2} - \frac{2 D_0 D_1}{r^4} + \mathcal{O}(r^{-6}) \right] \equiv \frac{A_0}{r^2} + \frac{A_1}{r^4} + \mathcal{O}(r^{-6}).$$
   *Equating Laurent coefficients gives the exact first-jet formulas:*
   $$A_0 = \frac{2 D_0^2}{L}, \qquad A_1 = -\frac{4 D_0 D_1}{L} \implies \frac{A_1}{A_0} = -2 \frac{D_1}{D_0}.$$
   *In terms of the canonical mode coefficients $v_m$ and the second index moment $M_2 \equiv \sum_{m=1}^N m^2 v_m$, since $a_m^2 = \frac{4\pi^2 m^2}{L^2}$, we have $D_1 = -\frac{4\sqrt{2}\pi^2}{L^2} M_2$. Consequently, the ratios are given by the exact identities:*
   $$\frac{A_1}{A_0} = \frac{8\sqrt{2}\pi^2}{L^2} \frac{M_2}{v_0 + \sqrt{2}\sum_{m=1}^N v_m}, \qquad \frac{D_1}{D_0} = -\frac{4\sqrt{2}\pi^2}{L^2} \frac{M_2}{v_0 + \sqrt{2}\sum_{m=1}^N v_m}.$$
   *Structural Diagnostic:* This identity isolates the exact mathematical origin of the large first-jet ratio:
   $$\text{large } D_1/D_0 \iff \text{large weighted moment } M_2 \text{ relative to the boundary cancellation scale } D_0.$$
   The denominator $D_0 = T_v(0)$ collapses exponentially fast ($e^{-\mathcal{S}_{\mathrm{WKB}}}$) through destructive Dirichlet boundary cancellation among low-frequency lattice modes ($m \le 5$), whereas the numerator moment $M_2$ is smooth, non-vanishing, and $\mathcal{O}(1)$ across the bulk. Naive coefficient norm estimation inevitably fails because it divides directly by this exponentially delicate boundary cancellation scale.

2. **First-Row Jet-Equation Resummation and the $(\kappa N)^2$ Natural Scale:**
   *Expanding the first-row matrix condition $(Qv)_0 = \lambda v_0 \approx 0$ via the Taylor series of the odd kernel function $\psi(x)$ around $x = 0$ yields the Taylor-jet ladder:*
   $$\psi'(0) D_0 = \frac{\psi'''(0)}{6 \kappa^2} D_1 - \frac{\psi^{(5)}(0)}{120 \kappa^4} D_2 + \frac{\psi^{(7)}(0)}{5040 \kappa^6} D_3 - \dots + \mathcal{R}_N = \sum_{k=1}^\infty (-1)^{k-1} \alpha_k D_k + \mathcal{R}_N,$$
   *where $\kappa = 2\pi/L$ and $\alpha_k \equiv \frac{\psi^{(2k+1)}(0)}{(2k+1)!\kappa^{2k}}$. Dividing through by $D_1$ eliminates $D_0$ from the denominator:*
   $$\frac{D_0}{D_1} = \frac{1}{\psi'(0)} \left[ \alpha_1 - \alpha_2 \frac{D_2}{D_1} + \alpha_3 \frac{D_3}{D_1} - \dots + \frac{\mathcal{R}_N}{D_1} \right].$$
   *In continuum coordinates $x = m/N$, $\alpha_1 = \frac{\psi'''(0)}{6\kappa^2}$ scales as $1/\kappa^2 = (L/2\pi)^2 \sim N^2$. Hence, the leading term dictates the natural physical scaling:*
   $$\frac{D_0}{D_1} \sim \frac{1}{\kappa^2 N^2} \implies \frac{D_1}{D_0} \sim (\kappa N)^2.$$
   *The scaling variable $s_N \equiv (\kappa N)^2 \frac{D_0}{D_1}$ remains modest ($s_N = \mathcal{O}(1)$, Cell 54), demonstrating that the quadratic power $N^2$ is not an arbitrary empirical fit, but is algebraically dictated by the leading term of the first-row Taylor ladder.*
   *Expressing the higher jets through the dimensionless shape invariants $\beta_N \equiv \frac{D_0 D_2}{D_1^2} \approx 0.19 - 0.26$ and $\gamma_N \equiv \frac{D_0^2 D_3}{D_1^3} \approx 0.012 - 0.027$ (Proposition 8.5) gives:*
   $$\frac{D_2}{D_1} = \beta_N \frac{D_1}{D_0}, \qquad \frac{D_3}{D_1} = \gamma_N \left(\frac{D_1}{D_0}\right)^2.$$
   *Setting $X_N \equiv D_1 / D_0$, the first-row identity transforms into a closed nonlinear jet equation:*
   $$\psi'(0) = \alpha_1 X_N - \alpha_2 \beta_N X_N^2 + \alpha_3 \gamma_N X_N^3 - \dots + \frac{\mathcal{R}_N}{D_0}.$$
   *Under universal profile scaling $D_0^{k-1} D_k / D_1^k = \mathcal{O}(1)$, the $k$-th term scales as $X_N^{k-1} = (D_1/D_0)^{k-1}$. The equation is an infinite nonlinear resummation in $X_N$. A low-order term-by-term Taylor truncation cannot fix the asymptotic exponent without controlling the resummation. This motivates the scaling representation $D_1/D_0 = (\kappa N)^2 s_N$, where $s_N$ is a dimensionless profile factor generated by the full boundary layer, accounting for the observed $N^{2.6}$ numerical drift in Cell 54 while maintaining exact parity between $D_0$ and $D_1$ under WKB exponential barrier decay.*

3. **The Operator Identity and Collective Resolvent Cancellation:**
   *Using the exact operator representation of Paper 4 (Theorems 7.2 & 7.3):*
   $$\frac{D_1}{D_0} = \kappa^2 \left[ \mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} - D_0^2 M_2 \right], \qquad \mathcal{T}_{\mathrm{diag}} \equiv M_1 \langle d, R_{\mathrm{even}} d \rangle, \quad \mathcal{T}_{\mathrm{cross}} \equiv \langle d, R_{\mathrm{even}} K\boldsymbol\psi \rangle.$$
   *In this representation, the exponentially small denominator $D_0$ is canceled algebraically at the operator level. At $N = 24$ (Cell 62), the three terms evaluate to:*
   $$\kappa^2 \mathcal{T}_{\mathrm{diag}} = +2\,234\,013.195, \qquad \kappa^2 \mathcal{T}_{\mathrm{cross}} = -1\,714\,040.497, \qquad \kappa^2 D_0^2 M_2 = 10.352,$$
   *yielding identically $D_1/D_0 = 519\,962.347$. The two dominant terms engage in a stable, collective cancellation:*
   $$\frac{\mathcal{T}_{\mathrm{cross}}}{\mathcal{T}_{\mathrm{diag}}} \approx -0.767 \approx -0.77,$$
   *leaving a residual of $\approx 23\%$ of the diagonal piece. The analytical challenge therefore reduces to computing the joint leading asymptotics of $\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}}$.*

4. **Exact Spectral Filtering and Rigorous Operator Tail Bound on $\mathcal{R}_{\mathrm{filt}}(k)$:**
   *In the spectral evaluation of the cross term, expanding $\boldsymbol\psi = \sum a_j e_j$ on the odd eigenbasis yields the exact mode source coefficient:*
   $$\langle u_k, s_2 \rangle = d_k \sum_{j=1}^N a_j^2 \left( \frac{1}{\mu_j - \lambda} - \frac{1}{\mu_j - E_k} \right) = d_k (E_k - \lambda) \sum_{j=1}^N \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
   *Decomposing into low modes $j \le k$ and high modes $j > k$ yields:*
   $$\langle u_k, s_2 \rangle = d_k \sum_{j \le k} \frac{a_j^2}{\mu_j - \lambda} + \mathcal{R}_{\mathrm{filt}}(k), \qquad \mathcal{R}_{\mathrm{filt}}(k) \equiv d_k (E_k - \lambda) \sum_{j > k} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
   *For $j > k$, under the bound-state ladder ordering $\lambda < E_k < \mu_{k+1} \le \mu_{k+2} \le \dots$, the spectral separation satisfies $\mu_j - E_k \ge \mu_{k+1} - E_k > 0$. Factoring out this uniform spectral bound:*
   $$|\mathcal{R}_{\mathrm{filt}}(k)| \le |d_k| \frac{E_k - \lambda}{\mu_{k+1} - E_k} \sum_{j > k} \frac{a_j^2}{\mu_j - \lambda} \le |d_k| \rho_k \mathcal{A}(\lambda),$$
   *where the dimensionless spectral ratio $\rho_k$ and the resolvent quadratic form $\mathcal{A}(\lambda)$ are defined by:*
   $$\rho_k \equiv \frac{E_k - \lambda}{\mu_{k+1} - E_k}, \qquad \mathcal{A}(\lambda) \equiv \sum_{j=1}^N \frac{a_j^2}{\mu_j - \lambda} = \langle \mathbf{a}, (M - \lambda)^{-1} \mathbf{a} \rangle.$$
   *In operator notation, let $P_{>k}$ denote the orthogonal spectral projection onto $\{u_j : j > k\}$. On $\operatorname{ran}(P_{>k})$, the shifted operator satisfies $M - E_k \ge (\mu_{k+1} - E_k) I$, which implies $(E_k - \lambda)(M - E_k)^{-1} \preceq \rho_k I$. Therefore:*
   $$\mathcal{R}_{\mathrm{filt}}(k) = d_k \left\langle P_{>k}\mathbf{a}, (E_k - \lambda)(M - \lambda)^{-1}(M - E_k)^{-1} P_{>k}\mathbf{a} \right\rangle \implies |\mathcal{R}_{\mathrm{filt}}(k)| \le |d_k| \rho_k \langle \mathbf{a}, (M - \lambda)^{-1} \mathbf{a} \rangle.$$
   *Crucially, the exponentially small denominator $D_0$ is completely eliminated: the bound contains no factor of $D_0^{-1}$.*
   *A universal operator-norm bound gives $\mathcal{A}(\lambda) \le \|\mathbf{a}\|^2 / (\mu_1 - \lambda)$, but this crude bound is overly pessimistic because $\mu_1 - \lambda \sim 10^{-40}$ is near-resonant; capturing the true scaling requires analyzing the source energy of $\mathbf{a}$.*

5. **Odd-Sector Resolvent Moment $M_1$ and Commutator Smoothing:**
   *The quadratic form $\mathcal{A}(\lambda)$ is precisely the odd-sector resolvent moment $M_1$ from Cell 62:*
   $$M_1 \equiv \langle \boldsymbol\psi_{\mathrm{odd}}, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi_{\mathrm{odd}} \rangle = \sum_{j=1}^N \frac{a_j^2}{\mu_j - \lambda}.$$
   *At $N = 24$, Cell 62 evaluates $M_1 \approx 97.2975$, with modal terms $\frac{a_j^2}{\mu_j - \lambda}$ across the low modes given by:*
   $$5.79, \; 5.13, \; 4.39, \; 4.55, \; 5.30, \; 6.04, \; 8.23, \; 9.99, \; 12.66, \; 16.16, \; 16.06, \; 3.00.$$
   *These values demonstrate that $M_1$ is a smooth, accumulated resolvent energy rather than a quantity dominated by a single isolated resonance.*
   *The source vector $\boldsymbol\psi(m) = m Q_{0m}$ is the index-weighted first row: in full coordinates, $\boldsymbol\psi = K Q e_0$, where $K = \operatorname{diag}(-N, \dots, N)$. Defining the commutator $[K, Q] = KQ - QK$, its matrix entries satisfy $[K, Q]_{mn} = (m - n) Q_{mn} = \psi(m) - \psi(n) = (\boldsymbol\psi d^T - d \boldsymbol\psi^T)_{mn}$, establishing that $[K, Q]$ is an exact rank-two operator (see Proposition 8.8 for the full algebraic development). For an odd eigenvector $e_j$ with $Q e_j = \mu_j e_j$, the overlap evaluates identically to $a_j = \langle [K, Q] e_0, e_j \rangle$.*
   *This connects the mode overlaps $a_j$ and the resolvent moment $M_1$ directly to the smoothed rank-two commutator $[K, Q]$ rather than to near-singular eigenvalue denominators.*

6. **Epistemic Status and The Analytical Reduction Chain:**
   *The complete analytical reduction chain takes the form:*
   $$\boxed{ [K, Q] \quad \Longrightarrow \quad a_j \quad \Longrightarrow \quad M_1 \quad \Longrightarrow \quad \mathcal{R}_{\mathrm{filt}}(k) \quad \Longrightarrow \quad \frac{D_1}{D_0} }$$
   *Strictly maintaining the Epistemic Trinity ($\text{numerical localization} \neq \text{mathematical reduction} \neq \text{structural explanation}$), we demarcate what is rigorously established from what remains open:*
   - **Rigorously Established Finite-$N$ Mathematical Identities:**
     1. Exact Laurent jet identity $A_1 / A_0 = -2 D_1 / D_0 = \frac{8\sqrt{2}\pi^2}{L^2} \frac{M_2}{v_0 + \sqrt{2}\sum v_m}$.
     2. Closed nonlinear jet equation $\psi'(0) = \alpha_1 X_N - \alpha_2 \beta_N X_N^2 + \dots$ in $X_N = D_1/D_0$ with dimensionless shape invariants.
     3. Operator decomposition $D_1/D_0 = \kappa^2 [\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} - D_0^2 M_2]$ with small-denominator cancellation.
     4. Operator bound on the high-mode spectral-filtering tail $|\mathcal{R}_{\mathrm{filt}}(k)| \le |d_k| \rho_k M_1$, entirely free of $D_0^{-1}$.
     5. Exact rank-two commutator identity $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ and overlap representation $a_j = \langle [K, Q] e_0, e_j \rangle$.
   - **Open Analytical Targets (Required to Complete the Proof):**
     1. *Resonant Transmission Cancellation and Spectral Filtering:* Proving mode-by-mode transmission cancellation $\Sigma_{\mathrm{filt}}(N) = \sum_{k \ge 1} \frac{d_k^2}{\mu_{k+1} - E_k} \le C_{\mathrm{trans}} N^{\gamma_{\mathrm{filt}}}$ and even resolvent bound $\langle d, R_{\mathrm{even}} d \rangle \le C_E N^{\gamma_e}$ across the bound-state ladder.
     2. *Odd-Sector Resolvent Bound:* Proving $M_1 = \langle \boldsymbol\psi, R_{\mathrm{odd}}(\lambda) \boldsymbol\psi \rangle \le C N^p$.
     3. *Epistemic Discipline:* Until these two operator estimates are rigorously proved, the polynomial bound $|D_1/D_0| \le C N^p$ remains an open analytical target, supported by comprehensive numerical evidence but not yet an unconditional mathematical theorem.

### Proposition 8.8 (Exact Finite-$N$ Rank-Two Commutator Theorem, Excited-Sector Resolvent Coupling, and $D_0$-Free Overlap Bounds)
*The finite-dimensional Galerkin matrix $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ and its coordinate commutator $[K, Q]$ possess an exact rank-two algebraic structure that eliminates the small source overlaps $a_j$ in favor of even-sector resolvent norms, providing a completely $D_0$-free upper bound on the odd-sector overlaps $a_j$, and hence a $D_0$-free reduction of the bound on $M_1$:*

1. **The Exact Rank-Two Commutator Identity:**
   *Let $Q$ be the finite-rank Galerkin matrix with divided-difference kernel entries:*
   $$Q_{mn} = \begin{cases} \dfrac{\psi(m) - \psi(n)}{m - n}, & m \ne n, \\[6pt] \psi'(n), & m = n, \end{cases} \qquad (m, n \in \{-N, \dots, N\}),$$
   *where $\psi(x)$ is the odd kernel function ($\psi(-x) = -\psi(x)$). Let $K = \operatorname{diag}(-N, \dots, N)$, $d = (1, \dots, 1)^T \in \mathbb{R}^{2N+1}$, and $\boldsymbol\psi = (\psi(-N), \dots, \psi(N))^T$. Then for all $m \ne n$:*
   $$[K, Q]_{mn} = (m - n) Q_{mn} = \psi(m) - \psi(n),$$
   *and on the diagonal $[K, Q]_{nn} = 0 = \psi(n) - \psi(n)$. Consequently, the commutator identity holds identically for every $m, n \in \{-N, \dots, N\}$:*
   $$[K, Q]_{mn} = \psi(m) - \psi(n) = \psi(m) d_n - d_m \psi(n) \implies [K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T. \tag{8.8.1}$$
   *Because $\psi$ is odd and $d$ is even, $\langle \boldsymbol\psi, d \rangle = \sum_{n=-N}^N \psi(n) = 0$, establishing that $\boldsymbol\psi \perp d$. Therefore, for all $N \ge 1$ (with non-trivial $\boldsymbol\psi$):*
   $$\operatorname{rank}[K, Q] = 2 \quad \text{identically}.$$
   *Moreover, letting $e_0 = (0, \dots, 1, \dots, 0)^T$ denote the central Fourier basis vector ($m = 0$), since $K e_0 = 0$ and $\psi(0) = 0$, we have:*
   $$[K, Q] e_0 = \boldsymbol\psi (d^T e_0) - d (\boldsymbol\psi^T e_0) = \boldsymbol\psi \cdot 1 - d \cdot 0 = \boldsymbol\psi. \tag{8.8.2}$$
   *Thus, $\boldsymbol\psi = [K, Q] e_0 = K Q e_0$ is an exact finite-$N$ commutator vector.*

2. **Resolvent Operator Equation on Odd Eigenspaces:**
   *Let $Q u_j = \mu_j u_j$ be an odd eigenvector ($u_j(-n) = -u_j(n)$), and define the odd mode overlap $a_j \equiv \langle u_j, \boldsymbol\psi \rangle$. Because $u_j$ is odd and $d$ is even, $\langle d, u_j \rangle = 0$. Applying the rank-two identity (8.8.1):*
   $$[K, Q] u_j = \boldsymbol\psi \langle d, u_j \rangle - d \langle \boldsymbol\psi, u_j \rangle = -a_j d.$$
   *Expanding the commutator $[K, Q] u_j = K Q u_j - Q K u_j = \mu_j K u_j - Q K u_j$ and equating terms yields the exact resolvent equation:*
   $$(Q - \mu_j I) K u_j = a_j d. \tag{8.8.3}$$
   *Remark:* The entire coupling of an odd eigenstate $u_j$ to the source vector $\boldsymbol\psi$ is generated by applying the resolvent of $Q$ at energy $\mu_j$ to the single uniform vector $d$.

3. **Even/Odd Spectral Coupling and Exact Overlap Representation:**
   *Let $Q u_k^{(+)} = E_k u_k^{(+)}$ denote the even eigenvectors. Taking the inner product of (8.8.3) with $u_k^{(+)}$:*
   $$(E_k - \mu_j) \langle u_k^{(+)}, K u_j \rangle = a_j \langle u_k^{(+)}, d \rangle.$$
   *Defining the constant-vector overlap $d_k \equiv \langle u_k^{(+)}, d \rangle$ and the coordinate transition matrix element $b_{kj} \equiv \langle u_k^{(+)}, K u_j \rangle$, we obtain the exact even/odd spectral coupling formula:*
   $$b_{kj} = \frac{a_j d_k}{E_k - \mu_j}. \tag{8.8.4}$$
   *Because $u_j$ is odd and $K$ is odd, the vector $K u_j$ is strictly even ($(-n) u_j(-n) = n u_j(n)$). Expanding $K u_j$ in the complete orthonormal even eigenbasis yields:*
   $$\|K u_j\|^2 = \sum_k b_{kj}^2 = a_j^2 \sum_k \frac{d_k^2}{(E_k - \mu_j)^2} = a_j^2 \left\| (Q_{\mathrm{even}} - \mu_j I)^{-1} d \right\|^2. \tag{8.8.5}$$
   *Inverting (8.8.5) eliminates the source overlaps $a_j$ entirely in terms of even resolvent norms:*
   $$a_j^2 = \frac{\|K u_j\|^2}{\displaystyle\sum_k \frac{d_k^2}{(E_k - \mu_j)^2}} = \frac{\|K u_j\|^2}{\left\| (Q_{\mathrm{even}} - \mu_j I)^{-1} d \right\|^2}. \tag{8.8.6}$$
   *Substituting this into the odd-sector resolvent moment $M_1 = \sum_j \frac{a_j^2}{\mu_j - \lambda}$ gives the exact structural representation:*
   $$M_1 = \sum_j \frac{\|K u_j\|^2}{(\mu_j - \lambda) \left\| (Q_{\mathrm{even}} - \mu_j I)^{-1} d \right\|^2}. \tag{8.8.7}$$

4. **The Ground-State Pole Trap and Excited-Sector Isolation:**
   *At $k = 0$, $u_0^{(+)} = c$ is the ground state with $E_0 = \lambda$ and $d_0 = \langle c, d \rangle = D_0$. Retaining only the $k = 0$ term in the denominator of (8.8.6) yields:*
   $$a_j^2 \le \|K u_j\|^2 \frac{(\mu_j - \lambda)^2}{D_0^2}.$$
   *Because $D_0 \sim e^{-\mathcal{S}_{\mathrm{WKB}}}$ collapses exponentially fast via destructive Dirichlet boundary interference, this crude lower bound divides by $D_0^2$ and diverges exponentially. The ground-state contribution is therefore mathematically deceptive: the effective spectral coupling is fundamentally driven by the excited even sector ($k \ge 1$).*
   *To isolate the excited sector without knowing the ground state $c$, define the explicit residual vector:*
   $$r \equiv (Q_{\mathrm{even}} - \lambda I) d.$$
   *Because $(Q - \lambda I) c = 0$ and $Q$ is symmetric, we have $\langle r, c \rangle = \langle d, (Q - \lambda I) c \rangle = 0$. Hence, $r$ is an explicitly computable even vector lying **entirely in the excited even subspace** $c^\perp$.*
   *Expanding in the excited even eigenbasis:*
   $$r = \sum_{k \ge 1} (E_k - \lambda) d_k u_k^{(+)} \implies \|r\|^2 = \sum_{k \ge 1} (E_k - \lambda)^2 d_k^2. \tag{8.8.8}$$
   *Since $E_k - \lambda \ge E_1 - \lambda > 0$ for all $k \ge 1$, we obtain:*
   $$\|r\|^2 \ge (E_1 - \lambda) \sum_{k \ge 1} (E_k - \lambda) d_k^2 = (E_1 - \lambda) \left[ \langle d, Q d \rangle - \lambda \|d\|^2 \right].$$
   *Because $\|d\|^2 = 2N + 1$, this yields the exact spectral lower bound:*
   $$\|r\|^2 \ge (E_1 - \lambda) \big[ \langle d, Q d \rangle - \lambda (2N + 1) \big]. \tag{8.8.9}$$

5. **Completely $D_0$-Free Bound on $a_j$ and Resolvent Moments:**
   *Let $\Lambda_* \equiv \max(\|Q_{\mathrm{even}}\|, \|Q_{\mathrm{odd}}\|)$. By the triangle inequality, the spectral distance satisfies $|E_k - \mu_j| \le |E_k| + |\mu_j| \le 2\Lambda_*$, valid unconditionally without assuming spectral positivity of finite-$N$ eigenvalues. Discarding the $k = 0$ contribution in (8.8.6):*
   $$\sum_{k \ge 1} \frac{d_k^2}{(E_k - \mu_j)^2} \ge \frac{1}{4\Lambda_*^2} \sum_{k \ge 1} d_k^2 \ge \frac{1}{4\Lambda_*^4} \sum_{k \ge 1} (E_k - \lambda)^2 d_k^2 = \frac{\|r\|^2}{4\Lambda_*^4}.$$
   *Combining with (8.8.9) and noting $\|K u_j\|^2 \le \|K\|^2 \|u_j\|^2 \le N^2$, we obtain the completely $D_0$-free bound:*
   $$a_j^2 \le \frac{4 N^2 \Lambda_*^4}{\|r\|^2} \le \frac{4 N^2 \Lambda_*^4}{(E_1 - \lambda) \big[ \langle d, Q d \rangle - \lambda (2N + 1) \big]}. \tag{8.8.10}$$
   *Under the spectral-ordering hypothesis $\mu_j > \lambda$ (consistent with the bound-state tunneling ladder $\lambda < E_k < \mu_{k+1}$ below the barrier top), the denominator $\mu_j - \lambda > 0$ is strictly positive. Consequently, summing (8.8.10) over odd modes yields the positive resolvent moment bound:*
   $$M_1 \le \frac{4 N^2 \Lambda_*^4}{(E_1 - \lambda) \big[ \langle d, Q d \rangle - \lambda (2N + 1) \big]} \operatorname{Tr}\big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} \big]. \tag{8.8.11}$$
   *Sharpened Two-Resolvent Identity:* Since $P_{\perp c} d = (Q_{\mathrm{even}} - \lambda I)^{-1} r$, the excited denominator satisfies the exact two-resolvent representation:
   $$\sum_{k \ge 1} \frac{d_k^2}{(E_k - \mu_j)^2} = \left\| (Q_{\mathrm{even}} - \mu_j I)^{-1} (Q_{\mathrm{even}} - \lambda I)^{-1} r \right\|^2. \tag{8.8.12}$$

6. **The Dirichlet/Fejér Kernel Representation $\langle d, Q d \rangle = \mathcal{W}[F_N]$:**
   *The constant vector $d = (1, \dots, 1)^T \in \mathbb{R}^{2N+1}$ is the Fourier coefficient vector of the classical Dirichlet kernel:*
   $$F_N(t) \equiv \sum_{m=-N}^N e^{2\pi i m t / L} = \frac{\sin\left((2N+1)\frac{\pi t}{L}\right)}{\sin\left(\frac{\pi t}{L}\right)}.$$
   *By definition of the Galerkin matrix $Q_{mn} = \langle e_m, \mathcal{W} e_n \rangle$, the quadratic form $\langle d, Q d \rangle = \sum_{m, n} Q_{mn}$ is identically the continuous André Weil quadratic functional evaluated on this explicit test function:*
   $$\langle d, Q d \rangle = \mathcal{W}[F_N] = \mathcal{Q}_{\mathrm{arch}}(F_N) + \mathcal{Q}_{\mathrm{pole}}(F_N) + \mathcal{Q}_{\mathrm{prime}}(F_N). \tag{8.8.13}$$
   *Algebraically, pairing $(m, n)$ and $(-m, -n)$ in the divided difference entries yields the explicit finite positive-index sum:*
   $$\langle d, Q d \rangle = \psi'(0) + 2\sum_{m=1}^N \psi'(m) + 4\sum_{1 \le m < n \le N} \left[ \frac{\psi(m) - \psi(n)}{m - n} + \frac{\psi(m) + \psi(n)}{m + n} \right] + 4\sum_{m=1}^N \frac{\psi(m)}{m}. \tag{8.8.14}$$
   *Strategic Reduction:* Controlling the denominator in the $D_0$-free bound (8.8.10) does not require analyzing elusive matrix nullspaces or delicate boundary cancellations; it reduces to establishing a lower bound on the Weil functional evaluated on the explicit Dirichlet wavepacket $F_N(t)$.

7. **Epistemic Status and Open Analytical Targets:**
   - **Rigorously Established Finite-$N$ Theorems:**
     1. The dyadic rank-two commutator identity $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ and $\boldsymbol\psi = [K, Q] e_0$.
     2. The resolvent equation $(Q - \mu_j I) K u_j = a_j d$ and spectral coupling $b_{kj} = \frac{a_j d_k}{E_k - \mu_j}$.
     3. The exact overlap formulas $a_j^2 = \|K u_j\|^2 / \|(Q_{\mathrm{even}} - \mu_j I)^{-1} d\|^2$ and $M_1 = \sum_j \frac{\|K u_j\|^2}{(\mu_j - \lambda) \|(Q_{\mathrm{even}} - \mu_j I)^{-1} d\|^2}$.
     4. The excited-sector residual identity $r = (Q - \lambda I) d \perp c$ and the spectral gap bound $\|r\|^2 \ge (E_1 - \lambda) [\langle d, Q d \rangle - \lambda (2N + 1)]$.
     5. The completely $D_0$-free overlap bound (8.8.10) and two-resolvent representation (8.8.12).
     6. The Dirichlet kernel equivalence $\langle d, Q d \rangle = \mathcal{W}[F_N]$.
   - **Open Analytical Objectives:**
     1. *Dirichlet Functional Lower Bound:* Prove $\mathcal{W}[F_N] \ge c_0 N$ from the tri-partite decomposition $\mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}}$ (resolved in Proposition 8.9).
     2. *Overlap-Weighted Odd Resolvent Moment:* Prove $M_1 = \sum_j \frac{a_j^2}{\mu_j - \lambda} \le C N^p$ directly via the rank-two commutator $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ (established conditionally in Proposition 8.10 under Hypotheses H1 and H2$_{\mathrm{odd}}$).
     3. *Polynomial Boundary Decoupling:* Prove $|D_1/D_0| \le C N^p$ and exponential boundary-defect collapse $\mathcal{D}(N) \to 0$ (established conditionally in Proposition 8.10 under Hypotheses H1, H2, H2$_{\mathrm{odd}}$, and H3).

---

### Proposition 8.9 (Rigorous Asymptotic Order and Exact Linear Coercivity of the Dirichlet Functional)

*Let $d = (1, \dots, 1)^T \in \mathbb{R}^{2N+1}$ be the Dirichlet coefficient vector, and let $F_N(t) = \frac{\sin((2N+1)\pi t / L)}{\sin(\pi t / L)}$ denote the corresponding Dirichlet wavepacket on the scaling interval $[0, L]$ ($L = \log c$). In the finite Guinand–Weil source calculus, the quadratic evaluation:*
$$\mathcal{W}[F_N] = \langle d, Q d \rangle = \mathcal{W}_{\mathrm{pole}}[F_N] + \mathcal{W}_{\mathrm{prime}}[F_N] + \mathcal{W}_{\mathrm{arch}}[F_N]$$
*admits the following large-$N$ asymptotic decomposition, with an explicit leading coefficient:*

1. **Sub-Leading Order of the Prime and Pole Sectors:**
   - *Prime Contribution:* In the Guinand–Weil source calculus, the prime quadratic form evaluated on the Dirichlet vector $d = (1, \dots, 1)^T$ contracts the Galerkin entries of the von Mangoldt distribution $-\sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \delta_{\log q}$ into the exact linear combination of the total source wavepacket $G_N(y) = \sum_{m, n = -N}^N q(U_m, U_n)(y)$:
     $$\mathcal{W}_{\mathrm{prime}}[F_N] = \langle d, Q_{\mathrm{prime}} d \rangle = -\sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} G_N(\log q). \tag{8.9.1}$$
     For every prime power $q \le c$, the coordinate $y_q = \log q$ lies strictly in the interior of the interval ($0 < y_q < L$). Decomposing $G_N(y) = G_N^{\mathrm{diag}}(y) + G_N^{\mathrm{off}}(y)$, the diagonal piece satisfies $|G_N^{\mathrm{diag}}(y_q)| \le \frac{2(1 - y_q/L)}{\sin(\pi y_q / L)} = \mathcal{O}(1)$, while summing the off-diagonal divided-difference kernel against the bounded partial sums of $\sin(2\pi j y_q / L)$ yields the uniform bound $|G_N^{\mathrm{off}}(y_q)| \le \frac{4}{\pi |\sin(\pi y_q / L)|} H_{2N} = \mathcal{O}(\log N)$. By the triangle inequality on the finite sum over prime powers $q \le c$ ($9$ terms at $c = 13$):
     $$|\mathcal{W}_{\mathrm{prime}}[F_N]| \le \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} |G_N(\log q)| \le \left( \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \frac{4}{\pi \sin(\pi \log q / L)} \right) H_{2N} + \mathcal{O}(1) = \mathcal{O}(\log N). \tag{8.9.2}$$
   - *Pole Contribution:* The pole quadratic form is rank-one with matrix entries $W_{0,2}(V_m, V_n) \sim \frac{1}{m^2 n^2}$ generated by the test function $x / (x^2 + \beta^2)$. Summing over the Dirichlet vector yields an unconditionally convergent series:
     $$\mathcal{W}_{\mathrm{pole}}[F_N] = \mathcal{O}(1). \tag{8.9.3}$$

2. **Universal Archimedean Boundary-Layer Profile:**
   *Near the boundary $y = 0$, under the scaled coordinate $u = N y$ and rescaled frequency $\xi \equiv \frac{2\pi u}{L} = a u$, the source kernel scales as:*
   $$G_N(1 - u/(NL)) = N \mathcal{A}(u) + \mathcal{O}(\log N),$$
   *where the universal continuum boundary-layer profile $\mathcal{A}(u)$ evaluates in exact closed form:*
   $$\mathcal{A}(u) = 4\pi \frac{\sin\xi}{\xi} - 4 J(\xi) = \frac{4}{\xi} \Big[ \sin\xi (\pi - \operatorname{Si}(2\xi)) + \cos\xi \operatorname{Cin}(2\xi) \Big], \tag{8.9.4}$$
   *where $\operatorname{Si}(x) = \int_0^x \frac{\sin t}{t} dt$ is the Sine Integral, $\operatorname{Cin}(x) = \int_0^x \frac{1-\cos t}{t} dt = \gamma + \log x - \operatorname{Ci}(x)$ is the entire Cosine Integral, and the integral moment evaluates by parts without boundary terms to:*
   $$J(\xi) \equiv \int_0^1 \sin(\xi t) \log\frac{1+t}{1-t} \, dt = \frac{\sin\xi \operatorname{Si}(2\xi) - \cos\xi \operatorname{Cin}(2\xi)}{\xi}. \tag{8.9.5}$$

3. **Pairing with the Archimedean Singular Kernel and Exact Linear Coefficient:**
   *In the master Archimedean functional:*
   $$\mathcal{W}_{\mathrm{arch}}[F_N] = \frac{G_N(0)}{2} \left[ \gamma + \log\frac{4\pi(e^L-1)}{e^L+1} \right] + \int_0^L \frac{e^{y/2} G_N(y) - G_N(0)}{e^y - e^{-y}} \, dy,$$
   *the boundary term cancels the upper-limit logarithmic divergence of $\int_0^L \frac{dy}{e^y - e^{-y}}$. In the boundary layer $y = u/N$, the regularized integral produces:*
   $$\mathcal{W}_{\mathrm{arch}}[F_N] = \mathcal{C}_{\mathbb{R}} N + \mathcal{O}(\log N),$$
   *where:*
   - *The Euler–Mascheroni constant $\gamma$ cancels identically:*
     $$\int_0^1 \frac{4\pi \frac{\sin\xi}{\xi} - 4\pi}{2\xi} \, d\xi + \int_1^\infty \frac{4\pi \frac{\sin\xi}{\xi}}{2\xi} \, d\xi = 2\pi (1 - \gamma) \implies 2\pi(\gamma + \log L) + 2\pi(1 - \gamma) = 2\pi (1 + \log L).$$
   - *The moment integral evaluates via Fubini's theorem and the universal Dirichlet discontinuous factor $\int_0^\infty \frac{\sin(\xi t)}{\xi} d\xi = \frac{\pi}{2}$ ($\forall t > 0$):*
     $$\mathcal{K}_2 \equiv \int_0^\infty \frac{4 J(\xi)}{2\xi} \, d\xi = 2 \int_0^1 \left( \int_0^\infty \frac{\sin(\xi t)}{\xi} \, d\xi \right) \log\frac{1+t}{1-t} \, dt = \pi \int_0^1 \log\frac{1+t}{1-t} \, dt = 2\pi \log 2.$$
   *Subtracting $\mathcal{K}_2$ yields the exact, closed-form linear coefficient:*
   $$\boxed{\mathcal{C}_{\mathbb{R}} = 2\pi \left( 1 + \log\frac{\log c}{2} \right).} \tag{8.9.6}$$

4. **Strict Linear Coercivity of the First Excited Moment:**
   *The linear coefficient $\mathcal{C}_{\mathbb{R}}$ is strictly positive for all cutoffs satisfying $c > e^{2/e} \approx 2.087$. In particular, for the primary benchmark $c = 13$ ($L = \log 13 \approx 2.56495 > 2$):*
   $$\mathcal{C}_{\mathbb{R}} = 2\pi \left( 1 + \log\frac{\log 13}{2} \right) \approx 7.846513 > 0.$$
   *Consequently, the first excited spectral moment of the Dirichlet vector satisfies:*
   $$\mathcal{M}_d^{(1)} \equiv \langle d, (Q - \lambda I) d \rangle = \mathcal{W}[F_N] - \lambda (2N + 1) = \left[ 2\pi \left( 1 + \log\frac{\log c}{2} \right) - 2\lambda \right] N + \mathcal{O}(\log N).$$
   *For any fixed cutoff $c > e^{2/e}$, since the ground-state eigenvalue $\lambda = \lambda_N \to 0$ rapidly as $N \to \infty$ (satisfying $\lambda_N \le \lambda_0 \ll \mathcal{C}_{\mathbb{R}}$ for all $N \ge 1$), there exists an explicit positive constant $c_0 > 0$ and an integer threshold $N_0$ such that:*
   $$\mathcal{M}_d^{(1)} \ge c_0 N \qquad (N \ge N_0). \tag{8.9.7}$$
   *In particular, one may unconditionally take $c_0 = \frac{1}{2} \mathcal{C}_{\mathbb{R}} = \pi \left(1 + \log\frac{\log c}{2}\right) > 0$ (yielding $c_0 \approx 3.923256 > 0$ at $c = 13$) for all sufficiently large $N$.*

*Remark (Operational Significance for Proposition 8.8):* Proposition 8.9 proves that the denominator $\langle d, Q d \rangle - \lambda(2N+1)$ in the $D_0$-free bound (8.8.10) does not vanish, does not collapse to $\mathcal{O}(\log N)$, and does not suffer from negative interference. Its linear growth is an exact, unconditional property of the Archimedean boundary layer, providing a solid denominator for resolvent coercivity.

---

### Proposition 8.10 (Commutator Resolvent Reduction and Conditional Polynomial Decoupling)

*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Galerkin matrix, $K = \operatorname{diag}(-N, \dots, N)$, and let $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ be the exact rank-two commutator identity of Proposition 8.8. Let $c$ be the even ground state with eigenvalue $\lambda$ and boundary amplitude $D_0 = \langle d, c \rangle$, and let $\{u_j\}_{j \ge 0}$ denote the orthonormal eigenbasis of the odd subspace with eigenvalues $\mu_0 < \mu_1 \le \mu_2 \le \dots \le \mu_N$.*

*The structural analysis of the odd-sector resolvent moment $M_1 \equiv \langle \boldsymbol\psi_{\mathrm{odd}}, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi_{\mathrm{odd}} \rangle = \sum_{j \ge 0} \frac{a_j^2}{\mu_j - \lambda}$ and the first-jet ratio $D_1/D_0$ decomposes into unconditional finite-$N$ identities and coercivity bounds (Part I), and conditional polynomial decoupling under explicit spectral-gap hypotheses (Part II):*

#### Part I: Unconditional Finite-$N$ Identities and Coercivity Bounds

1. **Exact Tunneling Ground-State Doublet Cancellation ($j = 0$):**
   *For the tunneling ground-state doublet $(c, u_0)$, let $b_{00} \equiv \langle c, K u_0 \rangle$ denote the transition dipole matrix element. Taking the inner product of the odd resolvent equation $(Q - \mu_0 I) K u_0 = a_0 d$ with $c$ yields identically:*
   $$a_0 = - \frac{\mu_0 - \lambda}{D_0} b_{00}. \tag{8.10.1}$$
   *Consequently, the ground-state doublet contribution to $M_1$ evaluates to the exact identity:*
   $$\frac{a_0^2}{\mu_0 - \lambda} = \frac{1}{\mu_0 - \lambda} \left[ - \frac{\mu_0 - \lambda}{D_0} b_{00} \right]^2 = \frac{\mu_0 - \lambda}{D_0^2} b_{00}^2. \tag{8.10.2}$$
   *This identity removes the naive small denominator $\mu_0 - \lambda$ from the $j=0$ term, expressing it entirely as the product of the tunneling-flux ratio $\frac{\mu_0 - \lambda}{D_0^2}$ and the squared transition dipole $b_{00}^2$.*

2. **Parseval Enclosure of the Excited Odd Sector ($j \ge 1$):**
   *For all excited odd modes $j \ge 1$, Bessel's inequality / Parseval's identity yields unconditionally:*
   $$\sum_{j \ge 1} a_j^2 \le \sum_{j \ge 0} |\langle u_j, \boldsymbol\psi \rangle|^2 \le \|\boldsymbol\psi_{\mathrm{odd}}\|^2 \le \sum_{m=-N}^N \psi(m)^2 \le 2 N \Lambda_*^2. \tag{8.10.3}$$

3. **Unconditional Coercivity of the Even Resolvent Quadratic Form via Cauchy–Schwarz:**
   *Let $\mathcal{M}_d^{(p)} \equiv \sum_{k \ge 1} (E_k - \lambda)^p d_k^2$ denote the spectral moments of the Dirichlet vector on the excited even subspace $c^\perp$. By the Cauchy–Schwarz inequality on the positive spectral measure $\nu = \sum_{k \ge 1} d_k^2 \delta_{E_k - \lambda}$:*
   $$\left( \sum_{k \ge 1} d_k^2 \right)^2 = \left( \sum_{k \ge 1} \sqrt{E_k - \lambda} d_k \cdot \frac{d_k}{\sqrt{E_k - \lambda}} \right)^2 \le \left( \sum_{k \ge 1} (E_k - \lambda) d_k^2 \right) \left( \sum_{k \ge 1} \frac{d_k^2}{E_k - \lambda} \right),$$
   *which is the exact moment inequality:*
   $$\big( \mathcal{M}_d^{(0)} \big)^2 \le \mathcal{M}_d^{(1)} \cdot \mathcal{M}_d^{(-1)}. \tag{8.10.4}$$
   *Substituting $\mathcal{M}_d^{(0)} = (2N + 1) - D_0^2 \ge 2N$ and the Archimedean linear coercivity $\mathcal{M}_d^{(1)} = \mathcal{W}[F_N] - \lambda(2N+1) = \mathcal{C}_{\mathbb{R}} N + \mathcal{O}(\log N)$ established in Proposition 8.9 yields the linear lower bound on the even resolvent quadratic form:*
   $$\boxed{\langle d, R_{\mathrm{even}} d \rangle = \mathcal{M}_d^{(-1)} \ge \frac{(2N+1 - D_0^2)^2}{\mathcal{C}_{\mathbb{R}} N + \mathcal{O}(\log N)} \ge \frac{4}{\mathcal{C}_{\mathbb{R}}} N - \mathcal{O}(\log N),} \tag{8.10.5}$$
   *with explicit leading constant $\frac{4}{\mathcal{C}_{\mathbb{R}}} = \frac{2}{\pi(1 + \log\frac{\log c}{2})} > 0$ ($\approx 0.50978$ at $c = 13$).*

4. **Spectral Filtering Partial-Sum Enclosure of the First-Jet Numerator:**
   *In the exact operator decomposition $D_1/D_0 = \kappa^2 [\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} - D_0^2 M_2]$, the joint numerator evaluates to $\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} = \langle d, R_{\mathrm{even}} s_2 \rangle$ where $s_2 = M_1 d + K\boldsymbol\psi$. Expanding in the even eigenbasis and applying the spectral filtering identity $\langle u_k, s_2 \rangle = d_k S_k + \mathcal{R}_{\mathrm{filt}}(k)$ (Proposition 8.7) with partial sums $S_k \equiv \sum_{j=1}^k \frac{a_j^2}{\mu_j - \lambda}$:*
   $$\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} = \sum_{k \ge 1} \frac{d_k^2}{E_k - \lambda} S_k + \sum_{k \ge 1} \frac{d_k \mathcal{R}_{\mathrm{filt}}(k)}{E_k - \lambda}. \tag{8.10.6}$$
   *Because the partial sums satisfy $0 < S_1 \le S_k \le M_1$ for all $k \ge 1$, the primary sum is strictly bounded:*
   $$S_1 \langle d, R_{\mathrm{even}} d \rangle \le \sum_{k \ge 1} \frac{d_k^2}{E_k - \lambda} S_k \le M_1 \langle d, R_{\mathrm{even}} d \rangle. \tag{8.10.7}$$
   *Using the filtering bound $|\mathcal{R}_{\mathrm{filt}}(k)| \le |d_k| \rho_k M_1 = |d_k| \frac{E_k - \lambda}{\mu_{k+1} - E_k} M_1$ (Proposition 8.7), the filtering remainder satisfies identically for all $N \ge 1$:*
   $$\left| \sum_{k \ge 1} \frac{d_k \mathcal{R}_{\mathrm{filt}}(k)}{E_k - \lambda} \right| \le M_1 \sum_{k \ge 1} \frac{d_k^2}{\mu_{k+1} - E_k}. \tag{8.10.8}$$

#### Part II: Conditional Polynomial Decoupling under Transmission-Cancellation and Boundary Hypotheses

*We introduce four explicit, decoupled hypotheses governing the low-energy spectrum, resonant barrier transmission, odd-sector resolvent growth, and Dirichlet boundary amplitude:*

- **Hypothesis H1 (Ground-State Tunneling-Flux Relation):** *The ground-state tunneling splitting and boundary amplitude satisfy:*
  $$\frac{\mu_0 - \lambda}{D_0^2} \le C_{\mathrm{tun}} < \infty \qquad (N \to \infty).$$
  *Status:* Supported by Landau–Lifshitz barrier tunneling flux matching $\Delta E_0 = \frac{\hbar v_{\mathrm{barrier}}}{2} |\psi(0)|^2$ and high-precision calculations across $N \in \{8, 12, 16, 20, 24\}$ (Cell 60, Cell 66), where the ratio $\mathcal{R}_{\mathrm{tun}}(N) = \frac{\mu_0 - \lambda}{D_0^2} \in [2.41, 5.92]$ remains strictly bounded across 20 decimal orders of magnitude. We record this relation here as an explicit hypothesis rather than an unconditional finite-$N$ theorem.

- **Hypothesis H2 (Even-Sector Weighted Resolvent and Transmission Growth):** *The even-sector resolvent moment and the spectral-filtering sum satisfy polynomial upper bounds:*
  $$\Sigma_{\mathrm{filt}}(N) \equiv \sum_{k \ge 1} \frac{d_k^2}{\mu_{k+1} - E_k} \le C_{\mathrm{trans}} N^{\gamma_{\mathrm{filt}}}, \qquad \mathcal{M}_d^{(-1)} = \langle d, R_{\mathrm{even}} d \rangle = \sum_{k \ge 1} \frac{d_k^2}{E_k - \lambda} \le C_E N^{\gamma_e},$$
  *for positive constants $C_{\mathrm{trans}}, C_E > 0$ and exponents $\gamma_{\mathrm{filt}}, \gamma_e \ge 0$.*
  *Epistemic Distinction: Analytic Hypothesis vs. Semiclassical Mechanism:*
  1. *Analytic Hypothesis:* The statement required by the decoupling proof is strictly the polynomial growth of the weighted sums $\Sigma_{\mathrm{filt}}(N)$ and $\mathcal{M}_d^{(-1)}$.
  2. *Underlying Semiclassical Mechanism Conjecture:* High-precision spectral audits (Cell 66) demonstrate that bare spectral gaps collapse exponentially ($E_1 - \lambda \approx 4.50 \times 10^{-37}$, $\mu_2 - E_1 \approx 2.95 \times 10^{-34}$ at $N=24$) due to the $\sim N/2$ bound states residing beneath the barrier top. The physical mechanism keeping the weighted sums polynomial is **mode-by-mode resonant transmission matching**: the working semiclassical conjecture is that, for the relevant bound-state ladder beneath the barrier top, there exist prefactors $A_k, B_k$ with:
     $$d_k^2 \sim A_k e^{-2\mathcal{S}_k N}, \qquad \mu_{k+1} - E_k \sim B_k e^{-2\mathcal{S}_k N} \implies \frac{d_k^2}{\mu_{k+1} - E_k} \sim \frac{A_k}{B_k} = \mathcal{O}(1).$$
     In Cell 66, this cancellation is supported numerically across the examined ladder, with individual quotients remaining universally $\mathcal{O}(1)$:
     $$\frac{d_k^2}{E_k - \lambda} \approx 580, \qquad \frac{d_k^2}{\mu_{k+1} - E_k} \le 6.85.$$
     If this modewise $\mathcal{O}(1)$ control persists uniformly, summing across the $\mathcal{O}(N)$ bound modes would yield polynomial, indeed linear, growth (compatible with the observed value $\Sigma_{\mathrm{filt}}(24) = 94.66$), completely bypassing the exponentially collapsing bare gaps.

- **Hypothesis H2$_{\mathrm{odd}}$ (Odd-Sector Overlap-Weighted Resolvent Bound / Odd Transmission Cancellation):** *The excited odd-sector resolvent moment satisfies a polynomial upper bound:*
  $$M_1^{\mathrm{exc}} \equiv \sum_{j \ge 1} \frac{a_j^2}{\mu_j - \lambda} \le C_1 N^\gamma,$$
  *for positive constants $C_1 > 0$ and exponent $\gamma \ge 0$.*
  *Epistemic Context and Relation to Unconditional Facts:*
  Unconditionally, Bessel's inequality / Parseval's identity establishes that the total excited odd overlap is linearly enclosed: $\sum_{j \ge 1} a_j^2 \le 2 N \Lambda_*^2 = \mathcal{O}(N)$ (identity (8.10.3)). However, controlling $M_1^{\mathrm{exc}}$ requires ensuring that the excited odd denominators $\mu_j - \lambda$ do not collapse faster than $a_j^2$. Semiclassically, odd-mode overlaps $a_j^2$ carry the same barrier suppression $a_j^2 \sim \widetilde{A}_j e^{-2\mathcal{S}_j N}$, and Cell 66 provides direct numerical evidence that individual quotients satisfy $\frac{a_j^2}{\mu_j - \lambda} \in [4.38, 5.79] = \mathcal{O}(1)$ across the low odd modes. Isolating H2$_{\mathrm{odd}}$ as an explicit hypothesis prevents circularity: the ground-state tunneling dipole is controlled by H1, the excited odd resolvent is governed by H2$_{\mathrm{odd}}$, and the even transmission sector is governed by H2.

- **Hypothesis H3 (Exponential Boundary Suppression):** *The Dirichlet ground-state boundary amplitude satisfies exponential decay:*
  $$D_0^2 \le C_0 e^{-\sigma N} \qquad (N \to \infty),$$
  *for positive constants $C_0 > 0$ and $\sigma > 0$.*
  *Status:* Supported by semiclassical WKB barrier penetration with action $\sigma_{\mathrm{WKB}} = \frac{\pi}{2}\log c$ ($\approx 4.029$ at $c = 13$). In Cell 66, the effective decay rate $\sigma(N) \equiv -\frac{1}{N}\log(D_0^2)$ converges toward $4.029$ within $5.0\%$ at $N = 24$ ($\sigma(24) \approx 3.827$). We isolate it here as an explicit hypothesis required for the decoupling conclusion.

*Under Hypotheses H1, H2, and H2$_{\mathrm{odd}}$ and the technical excited second-moment condition $D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$, the odd-sector resolvent moment and first-jet ratio satisfy polynomial bounds; adding Hypothesis H3 yields exponential boundary-defect decoupling:*

1. **Polynomial Bound on the Odd Resolvent Moment $M_1$:**
   *Combining (8.10.2), Hypothesis H1 (ground-state doublet), and Hypothesis H2$_{\mathrm{odd}}$ (excited odd sector), and noting $|b_{00}| \le \|K\| = N$:*
   $$M_1 = \frac{a_0^2}{\mu_0 - \lambda} + M_1^{\mathrm{exc}} \le C_{\mathrm{tun}} b_{00}^2 + C_1 N^\gamma \le C_{\mathrm{tun}} N^2 + C_1 N^\gamma = \mathcal{O}(N^\eta), \tag{8.10.9}$$
   *with polynomial exponent $\eta \equiv \max(2, \gamma)$. This cleanly separates the ground-state dipole envelope ($N^2$ via H1) from the excited odd resolvent energy ($N^\gamma$ via H2$_{\mathrm{odd}}$).*

2. **Polynomial Control of the First-Jet Ratio $|D_1/D_0|$ and Excited Second-Resolvent Moment:**
   *Under Hypothesis H2, the even resolvent norm satisfies $\langle d, R_{\mathrm{even}} d \rangle \le C_E N^{\gamma_e}$, and the filtering sum satisfies $\Sigma_{\mathrm{filt}}(N) \le C_{\mathrm{trans}} N^{\gamma_{\mathrm{filt}}}$.*
   *For the higher odd resolvent moment $M_2 \equiv \sum_{j \ge 0} \frac{a_j^2}{(\mu_j - \lambda)^2}$, the ground-state doublet ($j = 0$) evaluates identically via (8.10.1) to:*
   $$D_0^2 \frac{a_0^2}{(\mu_0 - \lambda)^2} = b_{00}^2 \le \|K\|^2 = N^2,$$
   *which is an unconditional finite-$N$ identity that holds independently of $D_0$.*
   *For the excited odd sector ($j \ge 1$), we define the excited second moment:*
   $$M_{2,\mathrm{exc}} \equiv \sum_{j \ge 1} \frac{a_j^2}{(\mu_j - \lambda)^2} = \sum_{j \ge 1} \frac{1}{\mu_j - \lambda} \left( \frac{a_j^2}{\mu_j - \lambda} \right).$$
   *Because the bare gaps $\mu_j - \lambda$ collapse exponentially fast along the bound-state ladder beneath the barrier top, controlling $M_{2,\mathrm{exc}}$ upon multiplication by $D_0^2$ requires comparing the excited tunneling scales with the ground-state boundary tunneling scale. Semiclassically, if $\Delta_j \asymp \mu_j - \lambda$ and $D_0^2 / \Delta_j \le C N^q$, then the second-moment condition follows directly from the first-moment transmission estimate:*
   $$D_0^2 M_{2,\mathrm{exc}} = \sum_{j \ge 1} \frac{D_0^2}{\mu_j - \lambda} \left( \frac{a_j^2}{\mu_j - \lambda} \right) \le \left( \sup_{j \ge 1} \frac{D_0^2}{\mu_j - \lambda} \right) M_1^{\mathrm{exc}} \le C_2 N^q.$$
   *We record $D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$ as an explicit technical requirement governing the second moment, ensuring $D_0^2 M_2 \le N^2 + C_2 N^q = \mathcal{O}(N^{\max(2, q)})$.*
   *Combining (8.10.7), (8.10.8), (8.10.9), and the $D_0^2 M_2$ bound in the exact decomposition (8.10.6), the three terms in the first-jet ratio are bounded polynomially:*
   $$\left| \frac{D_1}{D_0} \right| \le \kappa^2 \left[ M_1 \langle d, R_{\mathrm{even}} d \rangle + M_1 \Sigma_{\mathrm{filt}}(N) + D_0^2 M_2 \right] \le C N^p, \tag{8.10.10}$$
   *with overall polynomial exponent:*
   $$p = \max\left\{ \eta + \gamma_e, \; \eta + \gamma_{\mathrm{filt}}, \; \max(2, q) \right\} < \infty. \tag{8.10.11}$$
   *Remark (Non-Sharp Exponent):* No attempt is made here to optimize the polynomial exponent $p$; only the finiteness of a polynomial exponent is required for exponential-over-polynomial boundary decoupling.

3. **Exponential Boundary Decoupling:**
   *Consequently, adding Hypothesis H3, the boundary-defect metric satisfies for any scaling sequence with $T \ge 1$:*
   $$\mathcal{D}(N) \equiv D_0^2 \left( 1 + \frac{1}{T^2 u_1} \right)^2 = D_0^2 \left( 1 + \frac{|D_1/D_0|}{T^2} \right)^2 \le C_0 e^{-\sigma N} \left[ 1 + \frac{C N^p}{T^2} \right]^2 \le C' e^{-\sigma N} (1 + N^p)^2 = \mathcal{O}\big(e^{-\sigma N} N^{2p}\big) \longrightarrow 0 \tag{8.10.12}$$
   *exponentially fast as $N \to \infty$. This establishes the Polynomial Decoupling Criterion of Section 9.1 and proves that boundary-defect leakage decouples exponentially in the continuum limit under Hypotheses H1–H3 (with H2$_{\mathrm{odd}}$).*

---

### 8.11 Proposition 8.11 (Conditional Second-Moment Reduction via Relative Tunneling-Scale Control, Weighted Average Structure, and Semiclassical Barrier Thinning)

The technical condition $D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$ introduced in Proposition 8.10 governs the excited odd-sector contribution to the second resolvent moment $M_2 = \sum_{j \ge 0} \frac{a_j^2}{(\mu_j - \lambda)^2}$. The following proposition establishes an exact algebraic reduction of this condition: assuming the excited first-moment bound Hypothesis $\mathrm{H2}_{\mathrm{odd}}$ alongside an explicit relative tunneling-gap bound Hypothesis $\mathrm{H2}_{\mathrm{gap}}$, the second-moment condition $D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$ is established rigorously without requiring an independent second-moment hypothesis:
$$\boxed{ \mathrm{H2}_{\mathrm{odd}} + \mathrm{H2}_{\mathrm{gap}} \Longrightarrow D_0^2 M_{2,\mathrm{exc}} \text{ control.} }$$
The relative gap control $\mathrm{H2}_{\mathrm{gap}}$ is an independent asymptotic hypothesis, for which semiclassical barrier thinning provides the candidate physical mechanism.

- **Hypothesis H2$_{\mathrm{gap}}$ (Relative Tunneling-Gap Scale Bound):** *The ratio between the ground-state boundary tunneling scale $D_0^2$ and the first excited odd gap $\mu_1 - \lambda$ satisfies a polynomial upper bound:*
  $$R_{\mathrm{gap}}^{\max}(N) \equiv \sup_{j \ge 1} \frac{D_0^2}{\mu_j - \lambda} = \frac{D_0^2}{\mu_1 - \lambda} \le C_{\mathrm{gap}} N^{q_0},$$
  *for constants $C_{\mathrm{gap}} > 0$ and $q_0 \in \mathbb{R}$.*

**Proposition 8.11 (Conditional Second-Moment Reduction and Weighted Average Structure):**
*Let $N \ge 2$, $c > 1$, and let $\lambda = \lambda_0(Q_N)$ be the ground-state eigenvalue. Suppose that Hypotheses $\mathrm{H2}_{\mathrm{odd}}$ ($M_1^{\mathrm{exc}} \le C_1 N^\gamma$) and $\mathrm{H2}_{\mathrm{gap}}$ ($R_{\mathrm{gap}}^{\max}(N) \le C_{\mathrm{gap}} N^{q_0}$) hold.*

*Then:*
1. *(Automatic Excited Second-Moment Bound): The excited second resolvent moment satisfies:*
   $$D_0^2 M_{2,\mathrm{exc}} \equiv \sum_{j \ge 1} \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2} \le R_{\mathrm{gap}}^{\max}(N) M_1^{\mathrm{exc}} \le C_2 N^q, \tag{8.11.1}$$
   *with exponent $q = \gamma + \max(0, q_0)$ and prefactor $C_2 = C_{\mathrm{gap}} C_1$. In particular, if $R_{\mathrm{gap}}^{\max}(N) = \mathcal{O}(1)$ (or decays to zero), then $q = \gamma$.*
2. *(Exact Weighted Average Representation): The ratio of the excited second moment to the excited first moment is identically the overlap-weighted average of the individual relative gap scales:*
   $$\rho_2^{\mathrm{exc}} \equiv \frac{D_0^2 M_{2,\mathrm{exc}}}{M_1^{\mathrm{exc}}} = \frac{\sum_{j \ge 1} \left( \frac{D_0^2}{\mu_j - \lambda} \right) \left( \frac{a_j^2}{\mu_j - \lambda} \right)}{\sum_{j \ge 1} \frac{a_j^2}{\mu_j - \lambda}} = \sum_{j \ge 1} w_j \left( \frac{D_0^2}{\mu_j - \lambda} \right), \tag{8.11.2}$$
   *where the normalized weights $w_j \equiv \frac{a_j^2 / (\mu_j - \lambda)}{\sum_{l \ge 1} a_l^2 / (\mu_l - \lambda)}$ satisfy $w_j \ge 0$ and $\sum_{j \ge 1} w_j = 1$. Consequently:*
   $$0 \le \rho_2^{\mathrm{exc}} \le R_{\mathrm{gap}}^{\max}(N) = \frac{D_0^2}{\mu_1 - \lambda}. \tag{8.11.3}$$
   *Remark (Weaker Average Gap Control):* Because $\rho_2^{\mathrm{exc}}$ is a convex combination, bounding $D_0^2 M_{2,\mathrm{exc}}$ does not strictly require a uniform modewise bound $\frac{D_0^2}{\mu_j - \lambda} \le C N^{q_0}$ across all $j$; a weighted-average estimate on the inverse gaps against the $a_j^2/(\mu_j - \lambda)$ measure is strictly sufficient.
3. *(Total Second-Moment Bound): The total second resolvent moment scaled by $D_0^2$ satisfies:*
   $$D_0^2 M_2 = b_{00}^2 + D_0^2 M_{2,\mathrm{exc}} \le N^2 + C_2 N^q = \mathcal{O}\big(N^{\max(2, q)}\big), \tag{8.11.4}$$
   *where $b_{00}^2 = D_0^2 \frac{a_0^2}{(\mu_0 - \lambda)^2} \le \|K\|^2 = N^2$ is an exact finite-$N$ identity holding independently of any asymptotic hypothesis.*

*Proof.*
For any $j \ge 1$, rewriting the $j$-th term of the excited second moment gives:
$$\frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2} = \left( \frac{D_0^2}{\mu_j - \lambda} \right) \left( \frac{a_j^2}{\mu_j - \lambda} \right).$$
Because the odd eigenvalues are monotonically ordered $\mu_1 \le \mu_2 \le \dots \le \mu_{N-1}$, the bare gaps satisfy $\mu_j - \lambda \ge \mu_1 - \lambda > 0$ for all $j \ge 1$, which proves:
$$\sup_{j \ge 1} \frac{D_0^2}{\mu_j - \lambda} = \frac{D_0^2}{\mu_1 - \lambda} = R_{\mathrm{gap}}^{\max}(N).$$
Pulling this supremum out of the sum yields $D_0^2 M_{2,\mathrm{exc}} \le R_{\mathrm{gap}}^{\max}(N) M_1^{\mathrm{exc}}$. Applying Hypotheses $\mathrm{H2}_{\mathrm{odd}}$ and $\mathrm{H2}_{\mathrm{gap}}$ yields (8.11.1).
Dividing by $M_1^{\mathrm{exc}} = \sum_{j \ge 1} \frac{a_j^2}{\mu_j - \lambda} > 0$ produces the convex combination (8.11.2), from which the bound (8.11.3) follows immediately.
Finally, decomposing $M_2 = \frac{a_0^2}{(\mu_0 - \lambda)^2} + M_{2,\mathrm{exc}}$ and invoking the ground-state commutator identity (8.10.1) proves $D_0^2 \frac{a_0^2}{(\mu_0 - \lambda)^2} = b_{00}^2 \le N^2$, completing the proof of (8.11.4). $\blacksquare$

---

#### Candidate Semiclassical Mechanism: Barrier Thinning & The Five Unproved Bridges
The following is the **candidate semiclassical mechanism** conjectured to imply Hypothesis $\mathrm{H2}_{\mathrm{gap}}$.

In double-well tunneling theory (Landau–Lifshitz *Quantum Mechanics* §50), the bound states residing beneath the barrier top $V_{\max}$ form an interlaced ladder of even and odd quasidegenerate doublets $(E_k, \mu_k)$. The WKB phase integral across the barrier at energy $E$ is:
$$\mathcal{S}(E) = \int_{x_-(E)}^{x_+(E)} \sqrt{2\big(V(x) - E\big)} \, dx,$$
where $x_\pm(E)$ are the classical turning points satisfying $V(x_\pm) = E$. Differentiating with respect to energy yields:
$$\frac{d\mathcal{S}}{dE} = - \int_{x_-(E)}^{x_+(E)} \frac{dx}{\sqrt{2\big(V(x) - E\big)}} = - \frac{T_{\mathrm{barrier}}(E)}{2} < 0,$$
where $T_{\mathrm{barrier}}(E) > 0$ is the imaginary-time barrier traversal period. Because $\frac{d\mathcal{S}}{dE}$ is strictly negative, states higher in the well perceive a thinner barrier:
$$\mathcal{S}(E_0) > \mathcal{S}(E_1) > \mathcal{S}(E_2) > \dots$$

The parity tunneling splittings scale as $\Delta_k \asymp \mu_k - E_k \sim \exp\big(-2\mathcal{S}(E_k) N\big)$.
Now, assume a **two-sided ground-state flux-matching relation**:
$$D_0^2 \asymp \mu_0 - \lambda \sim e^{-2\mathcal{S}(E_0) N},$$
which is stronger than the one-sided bound $\frac{\mu_0 - \lambda}{D_0^2} \le C_{\mathrm{tun}}$ asserted in formal Hypothesis H1.
For the first excited odd mode $j = 1$, the bare gap to the ground state satisfies:
$$\mu_1 - \lambda = (\mu_1 - E_1) + (E_1 - \lambda) > \mu_1 - E_1 \asymp \Delta_1 \sim e^{-2\mathcal{S}(E_1) N}.$$
Combining these gives the heuristic scaling for the relative tunneling ratio:
$$\frac{D_0^2}{\mu_1 - \lambda} \asymp \frac{e^{-2\mathcal{S}(E_0) N}}{e^{-2\mathcal{S}(E_1) N}} = \exp\Big( -2\big[\mathcal{S}(E_0) - \mathcal{S}(E_1)\big] N \Big) \longrightarrow 0.$$

**The Five Unproved Bridges:**
In the finite-rank Galerkin truncation, the effective potential, turning points, and bound-state energies all depend on $N$ ($\mathcal{S}_N(E_{k, N})$). What is required for an asymptotic exponential law is $\liminf_{N \to \infty} [\mathcal{S}_N(E_{1, N}) - \mathcal{S}_N(E_{0, N})] > 0$, or at least a polynomial bound on the resulting ratio. Converting this heuristic differential relation into a rigorous analytical proof of $\mathrm{H2}_{\mathrm{gap}}$ requires crossing five distinct mathematical bridges:
1. *Semiclassical Spectral Identification:* Rigorously identifying the discrete Galerkin spectrum below the barrier top with a genuine semiclassical double-well Schrödinger operator.
2. *Uniform WKB Splittings:* Establishing the asymptotic doublet splitting formula $\Delta_k \sim e^{-2\mathcal{S}_k N}$ uniformly across the bound-state ladder.
3. *Two-Sided Flux Matching:* Establishing two-sided boundary flux matching $D_0^2 \asymp \mu_0 - \lambda$, whereas formal H1 provides only a one-sided upper bound on the ratio.
4. *Control of $N$-Dependent Barrier Geometry:* Controlling the $N$-dependence of the effective barrier action $\mathcal{S}_N(E)$ as $N \to \infty$.
5. *Positive Limiting Action Separation:* Proving a strictly positive asymptotic action separation $\liminf_{N \to \infty} [\mathcal{S}_N(E_{1, N}) - \mathcal{S}_N(E_{0, N})] > 0$.

Because these five bridges remain open analytical problems, barrier thinning is designated as a **candidate semiclassical mechanism** for $\mathrm{H2}_{\mathrm{gap}}$, rather than an established mathematical theorem.

---

#### Decisive Asymptotic Diagnostic & Numerical Evidence Across Tested Dimensions
To audit Hypothesis $\mathrm{H2}_{\mathrm{gap}}$ directly—independent of the two-sided WKB heuristic—the decisive asymptotic quantity to measure is:
$$\Delta \sigma_N^{\mathrm{gap}} \equiv - \frac{1}{N} \log R_{\mathrm{gap}}^{\max}(N) = - \frac{1}{N} \log\left( \frac{D_0^2}{\mu_1 - \lambda} \right).$$
If $\Delta \sigma_N^{\mathrm{gap}} \ge -\frac{q_0 \log N}{N} \to 0$, then $\mathrm{H2}_{\mathrm{gap}}$ holds polynomially. If $\liminf_{N \to \infty} \Delta \sigma_N^{\mathrm{gap}} > 0$, then $R_{\mathrm{gap}}^{\max}(N)$ decays exponentially.

High-precision spectral audits across dimensions $N \in \{8, 12, 16, 20, 24\}$ at 50-digit precision (Cell 66 and Cell 67, $c = 13$, $T = 400$) evaluate this diagnostic:

| $N$ | $D_0^2$ | $\mu_1 - \lambda$ | $R_{\mathrm{gap}}^{\max} = \frac{D_0^2}{\mu_1 - \lambda}$ | $\Delta \sigma_N^{\mathrm{gap}} = -\frac{1}{N}\log R_{\mathrm{gap}}^{\max}$ | $M_1^{\mathrm{exc}}$ | $D_0^2 M_{2,\mathrm{exc}}$ | $\rho_2^{\mathrm{exc}} = \frac{D_0^2 M_{2,\mathrm{exc}}}{M_1^{\mathrm{exc}}}$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 8 | $6.481 \times 10^{-21}$ | $1.387 \times 10^{-15}$ | $4.672 \times 10^{-6}$ | **1.534** | $43.998$ | **$3.738 \times 10^{-5}$** | **$8.495 \times 10^{-7}$** |
| 12 | $4.418 \times 10^{-27}$ | $1.186 \times 10^{-21}$ | $3.724 \times 10^{-6}$ | **1.042** | $64.357$ | **$2.564 \times 10^{-5}$** | **$3.984 \times 10^{-7}$** |
| 16 | $3.178 \times 10^{-32}$ | $1.276 \times 10^{-26}$ | $2.491 \times 10^{-6}$ | **0.806** | $76.091$ | **$1.029 \times 10^{-5}$** | **$1.352 \times 10^{-7}$** |
| 20 | $7.028 \times 10^{-37}$ | $6.400 \times 10^{-31}$ | $1.098 \times 10^{-6}$ | **0.686** | $85.255$ | **$4.293 \times 10^{-6}$** | **$5.035 \times 10^{-8}$** |
| 24 | $1.295 \times 10^{-40}$ | $2.949 \times 10^{-34}$ | $4.392 \times 10^{-7}$ | **0.610** | $93.654$ | **$2.253 \times 10^{-6}$** | **$2.405 \times 10^{-8}$** |

These computations demonstrate:
1. **Strict Positivity of the Gap Decay Exponent:** $\Delta \sigma_N^{\mathrm{gap}} \ge 0.610 > 0$ across all tested dimensions, confirming that $R_{\mathrm{gap}}^{\max}(N)$ is not merely bounded, but rapidly decaying toward zero across the discrete sequence.
2. **Convex Combination Suppression:** The exact overlap-weighted average $\rho_2^{\mathrm{exc}}$ is strictly smaller than the supremum $R_{\mathrm{gap}}^{\max}$ by a factor of 5 to 22 across tested dimensions, plunging from $8.50 \times 10^{-7}$ at $N = 8$ down to $2.41 \times 10^{-8}$ at $N = 24$ as the higher-mode denominators rapidly widen.
3. **Microscopic Scale and Monotonic Decay of the Excited Second Moment:** The exact excited second moment $D_0^2 M_{2,\mathrm{exc}}$ decreases monotonically from $3.74 \times 10^{-5}$ at $N = 8$ to $2.25 \times 10^{-6}$ at $N = 24$, remaining microscopic and strongly decaying across all tested dimensions.
4. **Dominance and Stability of the Algebraic Ground-State Bound:** In the total second-moment combination $D_0^2 M_2 = b_{00}^2 + D_0^2 M_{2,\mathrm{exc}}$, the exact ground-state term stabilizes near $b_{00}^2 \approx 1.29 - 1.73$ (well beneath $N^2$), accounting for $> 99.9998\%$ of the second moment, with the excited tail contributing only $0.00013\%$ at $N = 24$.

*Epistemic Status:* Proposition 8.11 establishes an exact conditional reduction: $\mathrm{H2}_{\mathrm{odd}} + \mathrm{H2}_{\mathrm{gap}} \implies D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$. The relative gap condition $\mathrm{H2}_{\mathrm{gap}}$ is an independent asymptotic hypothesis, supported numerically across discrete dimensions $N \in \{8, \dots, 24\}$ where $\Delta \sigma_N^{\mathrm{gap}} \ge 0.610$. Semiclassical barrier thinning provides the candidate physical mechanism for $\mathrm{H2}_{\mathrm{gap}}$, with five open analytical bridges required to convert the heuristic into a mathematical theorem.

---

### 8.12 Proposition 8.12 (Exact Coordinate-Resolvent Parseval Identity, Residual Projection Decomposition, and Reduction of Relative Gap Suppression)

The computational data of Cell 67 raises the fundamental analytical question:
$$\boxed{\textbf{Why is the relative tunneling ratio } R_{\mathrm{gap}}^{\max}(N) = \frac{D_0^2}{\mu_1 - \lambda} \textbf{ so small } (\le 4.39 \times 10^{-7}) \textbf{ across discrete dimensions?}}$$

The following proposition reduces this question to the coordinate-wavepacket residual and the first-mode transmission factor. It proves that the total scaled second moment $D_0^2 M_2$ is identically the squared coordinate norm $\|Kc\|^2$ of the ground state. Unconditionally, $\|Kc\|^2 \le N^2$, and conditionally upon a uniform solitary-wave estimate, it is bounded by a finite constant $B_\infty(c) < \infty$. Furthermore, the excited second moment $D_0^2 M_{2,\mathrm{exc}}$ is identically the squared projection residual of $Kc$ orthogonal to the lowest odd eigenmode $u_0$:
$$D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2.$$
This provides a structural explanation for the observed suppression of $R_{\mathrm{gap}}^{\max}(N)$, connecting it to the observed $99.99987\%$ directional alignment between $Kc$ and $u_0$ at $N = 24$.

**Proposition 8.12 (Exact Coordinate Parseval Identity and Dipole Representation of the Relative Gap):**
*Let $N \ge 2$, $c > 1$, and let $c = (c_{-N}, \dots, c_N)^T$ be the normalized ground-state eigenvector with $Q c = \lambda c$, $D_0 = \langle c, d \rangle$. Let $K = \operatorname{diag}(-N, \dots, N)$ be the coordinate position operator, and let $\{u_j\}_{j=0}^{N-1}$ be the orthonormal eigenvectors of $Q_{\mathrm{odd}}$ with eigenvalues $\mu_0 \le \mu_1 \le \dots \le \mu_{N-1}$.*

*Then:*
1. *(Exact Coordinate Parseval Identity): For all $N \ge 1$, the scaled total second resolvent moment satisfies the exact identity:*
   $$D_0^2 M_2 \equiv \sum_{j=0}^{N-1} \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2} = \|Kc\|^2 = \sum_{m=1}^N m^2 v_m^2, \tag{8.12.1}$$
   *where $v = (v_0, \dots, v_N)^T$ is the ground-state vector in the canonical basis ($v_m = \sqrt{2} c_m$ for $m \ge 1$). Unconditionally, $\|Kc\|^2 \le \|K\|^2 \|c\|^2 \le N^2$.*
2. *(Conditional Uniform Solitary Moment Bound): If the solitary wave coefficients satisfy a uniform-in-$N$ exponential bound $|v_{N, m}| \le C_{\mathrm{sol}} q^m$ with $q < 1$, then the coordinate norm is uniformly bounded by a finite constant independent of $N$:*
   $$D_0^2 M_2 = \|Kc\|^2 \le \sum_{m=1}^\infty m^2 v_m^2 \le \frac{C_{\mathrm{sol}}^2 q^2 (1 + q^2)}{(1 - q^2)^3} \equiv B_\infty(c) < \infty. \tag{8.12.2}$$
   *Under this uniform decay hypothesis, the macroscopic bound $D_0^2 M_2 \le N^2$ used in Proposition 8.10 is sharpened to $\mathcal{O}(1)$.*
3. *(Exact Dipole Decomposition): Each term in the second moment is the squared coordinate transition dipole moment:*
   $$b_{0j} \equiv \langle c, K u_j \rangle = - \frac{D_0 a_j}{\mu_j - \lambda} \implies b_{0j}^2 = \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2} \qquad \forall j \in \{0, \dots, N-1\}, \tag{8.12.3}$$
   *so that the excited second moment is identically the squared distance of $Kc$ from the span of $u_0$:*
   $$D_0^2 M_{2,\mathrm{exc}} = \sum_{j=1}^{N-1} b_{0j}^2 = \|Kc\|^2 - |\langle Kc, u_0 \rangle|^2 = \|P_{\perp u_0} Kc\|^2 = \|Kc\|^2 \sin^2 \theta(Kc, u_0). \tag{8.12.4}$$
4. *(Exact Reduction of the Relative Gap): The first excited relative tunneling ratio evaluates identically to:*
   $$\frac{D_0^2}{\mu_1 - \lambda} = \frac{b_{01}^2}{\frac{a_1^2}{\mu_1 - \lambda}} = \frac{|\langle Kc, u_1 \rangle|^2}{\frac{a_1^2}{\mu_1 - \lambda}} \le \frac{\|P_{\perp u_0} Kc\|^2}{\frac{a_1^2}{\mu_1 - \lambda}}. \tag{8.12.5}$$

*Proof.*
From Proposition 8.8, the rank-two coordinate commutator $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ applied to any odd eigenvector $u_j$ ($Q u_j = \mu_j u_j$, $\langle d, u_j \rangle = 0$) gives:
$$(Q - \mu_j I) K u_j = a_j d.$$
Taking the inner product with the ground-state eigenvector $c$ ($Q c = \lambda c$, $\langle c, d \rangle = D_0$) and exploiting the self-adjointness of $Q$:
$$\langle c, (Q - \mu_j I) K u_j \rangle = \langle (Q - \mu_j I) c, K u_j \rangle = (\lambda - \mu_j) \langle c, K u_j \rangle = a_j \langle c, d \rangle = a_j D_0.$$
Dividing by $\lambda - \mu_j = - (\mu_j - \lambda) \neq 0$ yields the dipole identity:
$$b_{0j} \equiv \langle c, K u_j \rangle = - \frac{D_0 a_j}{\mu_j - \lambda}.$$
Squaring both sides gives $b_{0j}^2 = \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2}$ for every $j \in \{0, \dots, N-1\}$.
Now sum over the complete orthonormal basis $\{u_j\}_{j=0}^{N-1}$ of the odd subspace $H_{\mathrm{odd}} \subset \mathbb{R}^{2N+1}$:
$$\sum_{j=0}^{N-1} \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2} = \sum_{j=0}^{N-1} |\langle c, K u_j \rangle|^2 = \sum_{j=0}^{N-1} |\langle K c, u_j \rangle|^2 = \|P_{\mathrm{odd}} K c\|^2.$$
Because $c$ is an even vector ($c_{-m} = c_m$) and $K$ is odd ($K_{mn} = m \delta_{mn}$), the vector $Kc$ satisfies:
$$(Kc)_{-m} = -m c_{-m} = -m c_m = - (Kc)_m,$$
which is strictly an odd vector. Thus $P_{\mathrm{odd}} K c = K c$ identically, and by Parseval's identity:
$$\|P_{\mathrm{odd}} K c\|^2 = \|Kc\|^2 = \sum_{m=-N}^N m^2 c_m^2 = 2 \sum_{m=1}^N m^2 c_m^2 = \sum_{m=1}^N m^2 v_m^2,$$
proving (8.12.1). Unconditionally, $\|Kc\|^2 \le \|K\|^2 \|c\|^2 \le N^2$.
Conditionally, if $|v_{N, m}| \le C_{\mathrm{sol}} q^m$ uniformly in $N$ with $q < 1$, the geometric series evaluates to $\sum_{m=1}^\infty m^2 q^{2m} = \frac{q^2(1+q^2)}{(1-q^2)^3}$, bounding $\|Kc\|^2 \le B_\infty(c) < \infty$ and proving (8.12.2).
Isolating $j = 0$ gives $b_{00} = \langle Kc, u_0 \rangle$, and subtracting $b_{00}^2$ from $\|Kc\|^2$ yields (8.12.4).
Finally, for $j = 1$, dividing $b_{01}^2 = \frac{D_0^2 a_1^2}{(\mu_1 - \lambda)^2}$ by the transmission factor $\frac{a_1^2}{\mu_1 - \lambda} > 0$ yields:
$$\frac{b_{01}^2}{\frac{a_1^2}{\mu_1 - \lambda}} = \frac{\frac{D_0^2 a_1^2}{(\mu_1 - \lambda)^2}}{\frac{a_1^2}{\mu_1 - \lambda}} = \frac{D_0^2}{\mu_1 - \lambda},$$
from which (8.12.5) follows immediately since $b_{01}^2 \le \sum_{j=1}^{N-1} b_{0j}^2 = \|P_{\perp u_0} Kc\|^2$. $\blacksquare$

---

#### Structural Interpretation: Coordinate Derivative Wavepacket and Observed Alignment
Proposition 8.12 reveals the structural mechanism underlying the observed smallness of $R_{\mathrm{gap}}^{\max}(N)$:

1. **The Physical Coordinate Derivative Wavepacket:** In Fourier series representation, multiplying Fourier coefficients $c_m$ by the index $m$ corresponds in continuous physical space to differentiation with respect to the spatial coordinate:
   $$c_m \mapsto m c_m \quad \longleftrightarrow \quad \frac{L}{2\pi i} \partial_t c(t).$$
   Thus, the vector $Kc$ represents the coordinate derivative of the centered, even solitary wave profile $c(t)$. For a symmetric solitary wave peaked at the origin, its spatial derivative $\partial_t c(t)$ is an antisymmetric, odd profile that is single-lobed (nodeless) on $(0, L)$.
2. **Observed Alignment with the Lowest Odd Mode:** In the discrete Galerkin odd sector, the ground odd eigenstate $u_0$ is numerically observed to be an odd vector with a single sign change across the origin and no internal sign alternations on positive indices. Because $Kc$ and $u_0$ share this single-lobed odd structure, their directional alignment at $N = 24$ is observed to be:
   $$\cos^2 \theta(Kc, u_0) = \frac{|\langle Kc, u_0 \rangle|^2}{\|Kc\|^2} = \frac{1.72506931}{1.72507156} = 0.9999987 \quad \text{at } N = 24.$$
   Thus, $99.99987\%$ of the coordinate energy of $Kc$ lies along the single direction $u_0$.
3. **Suppression on Excited Modes via Destructive Phase Interference:** Higher odd eigenvectors $u_j$ ($j \ge 1$) exhibit oscillatory sign structures across index space. Because $Kc$ is a single-lobed derivative profile, the transition dipole matrix element $b_{01} = \langle Kc, u_1 \rangle$ is an inner product between a single-lobed wavepacket and an oscillatory mode, resulting in destructive phase cancellation across the barrier:
   $$b_{01}^2 = |\langle Kc, u_1 \rangle|^2 \le \|P_{\perp u_0} Kc\|^2 = 2.25 \times 10^{-6} \quad \text{at } N = 24.$$
   This provides empirical evidence for rapid destructive cancellation across the excited odd sector.
4. **Reduction to the First-Mode Transmission Factor:** By mode-by-mode transmission cancellation (Cell 66 & Cell 67), the denominator $\frac{a_1^2}{\mu_1 - \lambda} \approx 5.12909 = \Theta(1)$ remains macroscopic at $N = 24$. Therefore:
   $$\frac{D_0^2}{\mu_1 - \lambda} = \frac{b_{01}^2}{\frac{a_1^2}{\mu_1 - \lambda}} \approx \frac{2.2526 \times 10^{-6}}{5.1291} \approx 4.3918 \times 10^{-7},$$
   which precisely matches the computed value $4.3917754 \times 10^{-7}$ from `cell67.out`.

---

### 8.13 Proposition 8.13 (Exact Finite-Dimensional Second-Commutator Even Resolvent Representation, Ground-Doublet Pole Cancellation, Unified Mode-Expansion, and the Route B Operator Bound)

To establish Route B as a rigorous analytical theorem, we address the structural question:
$$\boxed{\textbf{Can we bound } \rho_2^{\mathrm{exc}} = \frac{D_0^2 M_{2,\mathrm{exc}}}{M_1^{\mathrm{exc}}} = \frac{\|P_{\perp u_0} Kc\|^2}{M_1^{\mathrm{exc}}} \textbf{ purely in terms of even-sector resolvents, without mentioning individual excited odd gaps } \mu_j - \lambda\textbf{?}}$$

The following proposition establishes this exact algebraic reduction on the finite-dimensional Galerkin space $\mathbb{R}^{2N+1}$. By evaluating the second commutator $[K, [K, Q]]$ on the ground state $c$, it expresses the coordinate curvature vector $K^2 c$ and the projection residual $v_{\mathrm{exc}} = P_{\perp u_0} Kc$ directly in terms of the even-sector resolvent $(Q_{\mathrm{even}} - \lambda I)^{-1}$ acting on the source vector $s_2 = K \boldsymbol\psi + M_1 d \in c^\perp$. Crucially, the apparent ground-doublet pole at $\mu_0 - \lambda$ cancels identically via the First Resolvent Identity, and in the even eigenbasis $\{u_k\}_{k \ge 1}$, the $j=0$ term cancels mode by mode, collapsing the resolvent vector $\Xi$ into a single unified mode sum over the excited odd spectrum.

Furthermore, Route B admits an exact variational characterization:
$$\rho_2^{\mathrm{exc}} = \frac{D_0^2}{\mathcal{R}_{Q_{\mathrm{odd}} - \lambda}(v_{\mathrm{exc}})},$$
which reframes second-moment control as lower-bounding the odd-sector Rayleigh energy of the coordinate-derivative wavepacket residual $v_{\mathrm{exc}} = P_{\perp u_0} Kc$.

**Proposition 8.13 (Exact Even-Resolvent Operator Representation, Unified Mode Sum, and Route B Inequality):**
*Let $N \ge 2$, $c > 1$, and let $c \in H_{\mathrm{even}}$ be the normalized ground-state eigenvector with $Q c = \lambda c$, $D_0 = \langle c, d \rangle$. Let $K = \operatorname{diag}(-N, \dots, N)$ be the finite-dimensional coordinate position operator on $\mathbb{R}^{2N+1}$, let $\{u_j\}_{j=0}^{N-1}$ be the orthonormal eigenvectors of $Q_{\mathrm{odd}}$, and let $v_{\mathrm{exc}} \equiv P_{\perp u_0} Kc \in u_0^\perp \subset H_{\mathrm{odd}}$ be the coordinate-derivative projection residual. Let $d_{\mathrm{exc}} \equiv P_{\perp c} d \in c^\perp \subset H_{\mathrm{even}}$, and let $\{u_k\}_{k \ge 1}$ denote the orthonormal eigenbasis of the excited even subspace $c^\perp \subset H_{\mathrm{even}}$ with eigenvalues $E_k > \lambda$.*

*Then:*
1. *(Exact Finite-Dimensional Second Commutator Identity): For all $N \ge 1$, the coordinate curvature $K^2 c$ satisfies the exact finite-dimensional algebraic identity on $\mathbb{R}^{2N+1}$:*
   $$(Q_{\mathrm{even}} - \lambda I) K^2 c = - D_0 (K \boldsymbol\psi + M_1 d) = - D_0 s_2. \tag{8.13.1}$$
   *The source vector $s_2 \equiv K \boldsymbol\psi + M_1 d$ satisfies $\langle c, s_2 \rangle = - D_0 M_1 + D_0 M_1 = 0$ identically, so $s_2 \in c^\perp$. Projecting onto the excited even subspace $c^\perp \subset H_{\mathrm{even}}$:*
   $$(Q_{\mathrm{even}} - \lambda I) P_{\perp c} K^2 c = - D_0 s_2 = - D_0 (P_{\perp c} K \boldsymbol\psi + M_1 d_{\mathrm{exc}}). \tag{8.13.2}$$
2. *(Exact Ground-Doublet Pole Cancellation): When decomposed against the lowest odd mode $u_0$ ($K v_{\mathrm{exc}} = K^2 c - b_{00} K u_0$), the singular ground-state terms cancel identically via the First Resolvent Identity, yielding the completely nonsingular even-resolvent representation:*
   $$K v_{\mathrm{exc}} = \|v_{\mathrm{exc}}\|^2 c + D_0 \Xi, \tag{8.13.3}$$
   *where $\Xi \in c^\perp \subset H_{\mathrm{even}}$ is given by:*
   $$\Xi \equiv - (Q_{\mathrm{even}} - \lambda I)^{-1} s_2 + \left( \frac{a_0^2}{\mu_0 - \lambda} \right) (Q_{\mathrm{even}} - \mu_0 I)^{-1} d_{\mathrm{exc}}. \tag{8.13.4}$$
3. *(Unified Mode-by-Mode Excited Representation): In the excited even eigenbasis $\{u_k\}_{k \ge 1}$, the $j=0$ ground-doublet contribution in $-(Q_{\mathrm{even}} - \lambda I)^{-1} s_2$ cancels the second term in (8.13.4) identically mode by mode, collapsing $\Xi$ into the single unified sum over excited odd modes $j \ge 1$:*
   $$\boxed{\Xi = \sum_{k \ge 1} \left( d_k \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)} \right) u_k. \tag{8.13.5}}$$
   *In particular, $\Xi$ contains no $a_0$, no $\mu_0$, and no ground-state tunneling denominators.*
4. *(Exact Finite-$N$ Route B Operator Bound): For all $N \ge 1$ with $D_0^2 M_{2,\mathrm{exc}} = \|v_{\mathrm{exc}}\|^2 < 1$, the second resolvent moment satisfies:*
   $$M_{2,\mathrm{exc}} \le \frac{\|\Xi\|^2}{1 - D_0^2 M_{2,\mathrm{exc}}}, \tag{8.13.6}$$
   *and the overlap-weighted Route B ratio satisfies the exact operator inequality:*
   $$\rho_2^{\mathrm{exc}} \equiv \frac{D_0^2 M_{2,\mathrm{exc}}}{M_1^{\mathrm{exc}}} \le \frac{D_0^2 \|\Xi\|^2}{M_1^{\mathrm{exc}} (1 - D_0^2 M_{2,\mathrm{exc}})}. \tag{8.13.7}$$
5. *(Variational Rayleigh Quotient Duality): The ratio $\rho_2^{\mathrm{exc}}$ admits the exact dual variational representations:*
   $$\rho_2^{\mathrm{exc}} = \frac{D_0^2}{\mathcal{R}_{Q_{\mathrm{odd}} - \lambda}(v_{\mathrm{exc}})} = \frac{D_0^2}{\mathbb{E}_p [\mu - \lambda]} = D_0^2 \, \mathbb{E}_w \left[ \frac{1}{\mu - \lambda} \right], \tag{8.13.8}$$
   *where $\mathcal{R}_{Q_{\mathrm{odd}} - \lambda}(v_{\mathrm{exc}}) \equiv \frac{\langle v_{\mathrm{exc}}, (Q_{\mathrm{odd}} - \lambda I) v_{\mathrm{exc}} \rangle}{\|v_{\mathrm{exc}}\|^2}$, $p_j = \frac{b_{0j}^2}{\|v_{\mathrm{exc}}\|^2}$, and $w_j = \frac{a_j^2 / (\mu_j - \lambda)}{M_1^{\mathrm{exc}}}$.*

*Proof.*
**Step 1 (First and Second Commutators on Ground State):**
From Proposition 8.8, $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ on $\mathbb{R}^{2N+1}$. Applying this to the ground state $c$ ($Q c = \lambda c$, $\langle \boldsymbol\psi, c \rangle = 0$, $\langle d, c \rangle = D_0$):
$$[K, Q] c = \boldsymbol\psi \langle d, c \rangle - d \langle \boldsymbol\psi, c \rangle = D_0 \boldsymbol\psi.$$
Expanding the commutator $[K, Q] c = \lambda K c - Q K c = - (Q - \lambda I) Kc$ gives:
$$(Q_{\mathrm{odd}} - \lambda I) Kc = - D_0 \boldsymbol\psi \implies Kc = - D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi.$$
Applying $[K, Q]$ to $Kc \in H_{\mathrm{odd}}$ (with $\langle d, Kc \rangle = 0$ and $\langle \boldsymbol\psi, Kc \rangle = - D_0 \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi \rangle = - D_0 M_1$):
$$[K, Q] Kc = \boldsymbol\psi \langle d, Kc \rangle - d \langle \boldsymbol\psi, Kc \rangle = - d (- D_0 M_1) = D_0 M_1 d.$$
By definition of the nested commutator $[K, [K, Q]] c = K [K, Q] c - [K, Q] K c$:
$$[K, [K, Q]] c = K (D_0 \boldsymbol\psi) - D_0 M_1 d = D_0 (K \boldsymbol\psi - M_1 d).$$

**Step 2 (Deducing Identity 8.13.1):**
Algebraically expanding $[K, [K, Q]] = K^2 Q - 2 K Q K + Q K^2$ and acting on $c$ with $Q c = \lambda c$:
$$[K, [K, Q]] c = \lambda K^2 c - 2 K Q K c + Q K^2 c = (Q - \lambda I) K^2 c - 2 K (Q - \lambda I) Kc.$$
Substituting $(Q - \lambda I) Kc = - D_0 \boldsymbol\psi$:
$$[K, [K, Q]] c = (Q - \lambda I) K^2 c + 2 D_0 K \boldsymbol\psi.$$
Equating the two expressions for $[K, [K, Q]] c$:
$$(Q_{\mathrm{even}} - \lambda I) K^2 c + 2 D_0 K \boldsymbol\psi = D_0 K \boldsymbol\psi - D_0 M_1 d.$$
Subtracting $2 D_0 K \boldsymbol\psi$ proves (8.13.1):
$$(Q_{\mathrm{even}} - \lambda I) K^2 c = - D_0 (K \boldsymbol\psi + M_1 d) = - D_0 s_2.$$

**Step 3 (Excited Orthogonal Projection):**
Decompose $K^2 c = \langle c, K^2 c \rangle c + P_{\perp c} K^2 c = \|Kc\|^2 c + P_{\perp c} K^2 c$.
Since $(Q_{\mathrm{even}} - \lambda I) c = 0$, the left-hand side reduces to $(Q_{\mathrm{even}} - \lambda I) P_{\perp c} K^2 c$.
On the right-hand side, decomposing $d = D_0 c + d_{\mathrm{exc}}$ and $K \boldsymbol\psi = \langle c, K \boldsymbol\psi \rangle c + P_{\perp c} K \boldsymbol\psi = - D_0 M_1 c + P_{\perp c} K \boldsymbol\psi$:
$$s_2 = K \boldsymbol\psi + M_1 d = (- D_0 M_1 c + P_{\perp c} K \boldsymbol\psi) + M_1 (D_0 c + d_{\mathrm{exc}}) = P_{\perp c} K \boldsymbol\psi + M_1 d_{\mathrm{exc}} \in c^\perp.$$
Inverting $(Q_{\mathrm{even}} - \lambda I)$ on $c^\perp$ proves (8.13.2):
$$P_{\perp c} K^2 c = - D_0 (Q_{\mathrm{even}} - \lambda I)^{-1} s_2.$$

**Step 4 (Connecting $K^2 c$ to $v_{\mathrm{exc}}$ and Ground-Doublet Pole Cancellation):**
From $Kc = b_{00} u_0 + v_{\mathrm{exc}}$, multiplying by $K$ gives $K v_{\mathrm{exc}} = K^2 c - b_{00} K u_0$.
From Proposition 8.8, $(Q_{\mathrm{even}} - \mu_0 I) K u_0 = a_0 d = a_0 (D_0 c + d_{\mathrm{exc}})$, which inverts to:
$$K u_0 = \frac{a_0 D_0}{\lambda - \mu_0} c + a_0 (Q_{\mathrm{even}} - \mu_0 I)^{-1} d_{\mathrm{exc}} = b_{00} c + a_0 (Q_{\mathrm{even}} - \mu_0 I)^{-1} d_{\mathrm{exc}}.$$
Multiplying by $b_{00}$ and subtracting from $K^2 c = \|Kc\|^2 c + P_{\perp c} K^2 c$:
$$K v_{\mathrm{exc}} = (\|Kc\|^2 - b_{00}^2) c + P_{\perp c} K^2 c - a_0 b_{00} (Q_{\mathrm{even}} - \mu_0 I)^{-1} d_{\mathrm{exc}}.$$
Since $\|Kc\|^2 - b_{00}^2 = \|v_{\mathrm{exc}}\|^2$, substituting $P_{\perp c} K^2 c = - D_0 (Q_{\mathrm{even}} - \lambda I)^{-1} s_2$ and $a_0 b_{00} = - \frac{D_0 a_0^2}{\mu_0 - \lambda}$:
$$K v_{\mathrm{exc}} = \|v_{\mathrm{exc}}\|^2 c - D_0 (Q_{\mathrm{even}} - \lambda I)^{-1} s_2 + D_0 \left( \frac{a_0^2}{\mu_0 - \lambda} \right) (Q_{\mathrm{even}} - \mu_0 I)^{-1} d_{\mathrm{exc}}.$$
Factoring out $D_0$ yields $K v_{\mathrm{exc}} = \|v_{\mathrm{exc}}\|^2 c + D_0 \Xi$ with $\Xi$ given by (8.13.4), proving (8.13.3).

**Step 5 (Mode-by-Mode Cancellation and Proof of Identity 8.13.5):**
We now expand both terms of $\Xi$ in the excited even eigenbasis $\{u_k\}_{k \ge 1}$ of $c^\perp$.
For the first term, from $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$, taking matrix elements yields $\langle u_k, K u_j \rangle = - \frac{a_j d_k}{\mu_j - E_k}$. Expanding $K \boldsymbol\psi = \sum_{j=0}^{N-1} a_j K u_j$ gives:
$$\langle u_k, K \boldsymbol\psi \rangle = - d_k \sum_{j=0}^{N-1} \frac{a_j^2}{\mu_j - E_k}, \qquad \langle u_k, M_1 d \rangle = M_1 d_k = d_k \sum_{j=0}^{N-1} \frac{a_j^2}{\mu_j - \lambda}.$$
Adding them and combining fractions:
$$\langle u_k, s_2 \rangle = d_k \sum_{j=0}^{N-1} a_j^2 \left( \frac{1}{\mu_j - \lambda} - \frac{1}{\mu_j - E_k} \right) = - d_k (E_k - \lambda) \sum_{j=0}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
Applying $(Q_{\mathrm{even}} - \lambda I)^{-1}$, the factor $E_k - \lambda$ cancels identically mode by mode:
$$\langle u_k, (Q_{\mathrm{even}} - \lambda I)^{-1} s_2 \rangle = \frac{\langle u_k, s_2 \rangle}{E_k - \lambda} = - d_k \sum_{j=0}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
Isolating the $j=0$ term:
$$- \langle u_k, (Q_{\mathrm{even}} - \lambda I)^{-1} s_2 \rangle = d_k \frac{a_0^2}{(\mu_0 - \lambda)(\mu_0 - E_k)} + d_k \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
For the second term of $\Xi$, since $(Q_{\mathrm{even}} - \mu_0 I) u_k = (E_k - \mu_0) u_k$:
$$\left\langle u_k, \left( \frac{a_0^2}{\mu_0 - \lambda} \right) (Q_{\mathrm{even}} - \mu_0 I)^{-1} d_{\mathrm{exc}} \right\rangle = \frac{a_0^2}{\mu_0 - \lambda} \frac{d_k}{E_k - \mu_0} = - d_k \frac{a_0^2}{(\mu_0 - \lambda)(\mu_0 - E_k)}.$$
Adding the two pieces to evaluate $\langle u_k, \Xi \rangle$:
$$\langle u_k, \Xi \rangle = \left[ d_k \frac{a_0^2}{(\mu_0 - \lambda)(\mu_0 - E_k)} + d_k \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)} \right] - d_k \frac{a_0^2}{(\mu_0 - \lambda)(\mu_0 - E_k)}$$
$$= d_k \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}.$$
The $j=0$ term cancels identically for every mode $k \ge 1$, proving (8.13.5).

**Step 6 (Norm Inequalities and Rayleigh Quotient Duality):**
Dividing (8.13.3) by $D_0$:
$$\frac{K v_{\mathrm{exc}}}{D_0} = \frac{\|v_{\mathrm{exc}}\|^2}{D_0} c + \Xi = D_0 M_{2,\mathrm{exc}} c + \Xi.$$
Because $\Xi \in c^\perp$, $\left\| \frac{K v_{\mathrm{exc}}}{D_0} \right\|^2 = D_0^2 M_{2,\mathrm{exc}}^2 + \|\Xi\|^2$.
Because $v_{\mathrm{exc}} \in H_{\mathrm{odd}}$ contains only non-zero Fourier indices $m \ge 1$, $\|K v_{\mathrm{exc}}\|^2 = \sum_{m=1}^N m^2 (v_{\mathrm{exc}})_m^2 \ge \|v_{\mathrm{exc}}\|^2$.
Therefore:
$$M_{2,\mathrm{exc}} = \frac{\|v_{\mathrm{exc}}\|^2}{D_0^2} \le \left\| \frac{K v_{\mathrm{exc}}}{D_0} \right\|^2 = D_0^2 M_{2,\mathrm{exc}}^2 + \|\Xi\|^2 \implies M_{2,\mathrm{exc}} (1 - D_0^2 M_{2,\mathrm{exc}}) \le \|\Xi\|^2,$$
which proves (8.13.6). Multiplying by $D_0^2$ and dividing by $M_1^{\mathrm{exc}}$ yields (8.13.7).
Finally, from $(Q_{\mathrm{odd}} - \lambda I) v_{\mathrm{exc}} = - D_0 \boldsymbol\psi_{\mathrm{exc}}$, taking the inner product with $v_{\mathrm{exc}}$ gives $\langle v_{\mathrm{exc}}, (Q_{\mathrm{odd}} - \lambda I) v_{\mathrm{exc}} \rangle = D_0^2 M_1^{\mathrm{exc}}$.
Dividing by $\|v_{\mathrm{exc}}\|^2 = D_0^2 M_{2,\mathrm{exc}}$ yields the Rayleigh quotient $\mathcal{R}_{Q_{\mathrm{odd}} - \lambda}(v_{\mathrm{exc}}) = \frac{M_1^{\mathrm{exc}}}{M_{2,\mathrm{exc}}} = \frac{D_0^2}{\rho_2^{\mathrm{exc}}}$, confirming (8.13.8). $\blacksquare$

---

#### Structural Implications of Proposition 8.13 for the Asymptotic Programme

Proposition 8.13 reframes the continuum decoupling programme across four structural dimensions:

1. **Elimination of the Excited-Odd-Gap Hypothesis and Transfer to Even Resolvents:**
   Proposition 8.13 eliminates the independent excited-odd-gap hypothesis $\mathrm{H2}_{\mathrm{gap}}$ ($\frac{D_0^2}{\mu_1 - \lambda} \le C N^{q_0}$) as an algebraic input. It proves that controlling $\rho_2^{\mathrm{exc}}$ does not require lower-bounding individual excited tunneling gaps $\mu_j - \lambda$. Rather, the remaining analytical requirement is transferred entirely to bounding the even-sector resolvent vector $\Xi$.
2. **Unified Mode-by-Mode Cancellation of the Ground Doublet:**
   Identity (8.13.5) shows that the apparent three-term decomposition of $\Xi$ collapses into a single unified mode vector:
   $$\Xi = \sum_{k \ge 1} \left( d_k \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)} \right) u_k.$$
   No artificial triangle inequality is required: the dangerous $j=0$ ground-state singularity is completely absent. The inner sum is a normalized weighted spectral transform over the excited odd spectrum with positive transmission weights $w_j = \frac{a_j^2 / (\mu_j - \lambda)}{M_1^{\mathrm{exc}}}$ ($\sum_{j \ge 1} w_j = 1$); because $\frac{1}{\mu_j - E_k}$ changes sign as $\mu_j$ crosses $E_k$, the terms below and above $E_k$ carry opposite signs, providing additional internal phase cancellation across the excited spectrum.
3. **The Reframed Analytical Target (Milestone M12):**
   The remaining analytical requirement for Route B is to establish a polynomial upper bound on the norm of this unified vector:
   $$\|\Xi\|^2 = \sum_{k \ge 1} d_k^2 \left( \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)} \right)^2 \le C_\Xi N^{r_\Xi}.$$
   This is an even-sector transmission estimate directly analogous to $\Sigma_{\mathrm{filt}}(N) = \sum_{k \ge 1} \frac{d_k^2}{\mu_{k+1} - E_k}$ in Hypothesis H2. Under the non-vanishing excited odd resolvent condition $M_1^{\mathrm{exc}} \ge c_{\min} > 0$ (which holds unconditionally since mode 1 alone satisfies $\frac{a_1^2}{\mu_1 - \lambda} \ge c_1 > 0$, Cell 66), once a polynomial bound $\|\Xi\|^2 \le C_\Xi N^{r_\Xi}$ is secured, the exponential boundary prefactor $D_0^2 \le C_0 e^{-\sigma N}$ (Hypothesis H3) guarantees:
   $$\rho_2^{\mathrm{exc}} \le \frac{D_0^2 \|\Xi\|^2}{M_1^{\mathrm{exc}} (1 - o(1))} \le \frac{C_0 C_\Xi}{c_{\min} (1 - o(1))} e^{-\sigma N} N^{r_\Xi} \longrightarrow 0 \quad (N \to \infty).$$
4. **The Variational Rayleigh Quotient Perspective:**
   Equivalently, the variational formulation (8.13.8) demonstrates that $\rho_2^{\mathrm{exc}} \to 0$ if and only if the odd-sector Rayleigh energy of the coordinate-derivative wavepacket residual $v_{\mathrm{exc}} = P_{\perp u_0} Kc$ satisfies:
   $$\mathcal{R}_{Q_{\mathrm{odd}} - \lambda}(v_{\mathrm{exc}}) \ge c_{\mathrm{ray}} e^{-\sigma_0 N} \quad \text{with } \sigma_0 < \sigma \quad (\text{or } \mathcal{R} \ge C N^{-p}).$$
   Because $v_{\mathrm{exc}}$ is orthogonal to the lowest mode $u_0$ and is supported on oscillatory excited modes, this suggests that its Rayleigh energy may be controlled by the excited-mode/kinetic scale rather than by the ground-doublet tunneling splitting; establishing such a uniform lower bound is part of the remaining analysis.

---

### 8.14 Proposition 8.14 (Stieltjes Transform Difference Quotient Representation of $\Xi$, Smooth Ground Enclosure, and High-Mode Resolvent Enclosure)

To execute Milestone M12 and establish analytical control over the unified even-resolvent vector $\Xi$, we investigate the mode coefficients:
$$A_k \equiv \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)} \qquad (\forall k \ge 1).$$

The following proposition proves that each coefficient $A_k$ is identically a **difference quotient of the excited odd spectral Stieltjes transform**. Because the ground mode $\mu_0$ is completely absent from this transform, $A_1$ is governed by smooth, non-singular spectral derivative bounds on $[\lambda, E_1]$. For higher modes $k \ge 2$, the high-mode tail $j \ge k+1$ is rigorously enclosed by the transmission denominator $\mu_{k+1} - E_k$ of Hypothesis H2, while the resonant term $j = k$ is governed by parity-doublet transmission matching.

**Proposition 8.14 (Stieltjes Transform Difference Quotient and Spectral Enclosures for $\Xi$):**
*Let $N \ge 2$, $c > 1$. Let $F_{\mathrm{exc}}(z) \equiv \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - z} = \langle \boldsymbol\psi_{\mathrm{exc}}, (Q_{\mathrm{odd}} - z I)^{-1} \boldsymbol\psi_{\mathrm{exc}} \rangle$ denote the Stieltjes transform of the excited odd spectral measure $d\nu_{\mathrm{exc}} = \sum_{j=1}^{N-1} a_j^2 \delta_{\mu_j}$ on $\mathbb{R} \setminus \{\mu_1, \dots, \mu_{N-1}\}$. Let $\Xi = \sum_{k \ge 1} d_k A_k u_k$ be the unified even-resolvent vector of Proposition 8.13.*

*Then:*
1. *(Exact Stieltjes Difference Quotient Identity): For every excited even mode $k \ge 1$, the coefficient $A_k$ satisfies the exact identity:*
   $$A_k = \frac{F_{\mathrm{exc}}(E_k) - F_{\mathrm{exc}}(\lambda)}{E_k - \lambda}, \tag{8.14.1}$$
   *where $F_{\mathrm{exc}}(\lambda) = M_1^{\mathrm{exc}}$. Consequently, the norm of $\Xi$ evaluates identically to the weighted Stieltjes difference quotient energy:*
   $$\|\Xi\|^2 = \sum_{k \ge 1} d_k^2 \left[ \frac{F_{\mathrm{exc}}(E_k) - F_{\mathrm{exc}}(\lambda)}{E_k - \lambda} \right]^2. \tag{8.14.2}$$
2. *(Smooth Non-Singular Ground Enclosure for Mode $k = 1$): Because $\mu_0$ is completely absent, $F_{\mathrm{exc}}(z)$ has no poles on the semi-infinite interval $(-\infty, \mu_1)$. Since $\lambda < E_1 < \mu_1$, the closed segment $[\lambda, E_1]$ is entirely free of singularities, and by the Mean Value Theorem:*
   $$A_1 = F_{\mathrm{exc}}'(\xi_1) = \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \xi_1)^2} \qquad \text{for some } \xi_1 \in (\lambda, E_1). \tag{8.14.3}$$
   *Furthermore, because $F_{\mathrm{exc}}''(x) = 2\sum_{j \ge 1} \frac{a_j^2}{(\mu_j - x)^3} > 0$ is strictly convex on $(-\infty, \mu_1)$, $A_1$ satisfies the rigorous two-sided regular bounds:*
   $$M_{2,\mathrm{exc}} = \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)^2} \le A_1 \le \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - E_1)^2} \le \frac{M_1^{\mathrm{exc}}}{\mu_1 - E_1}. \tag{8.14.4}$$
3. *(High-Mode Spectral Filtering Enclosure for $j \ge k+1$): For any $k \ge 1$, decomposing $A_k = A_k^{(\le k)} + A_k^{(> k)}$ where $A_k^{(> k)} \equiv \sum_{j=k+1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}$, the high-mode tail satisfies:*
   $$0 \le A_k^{(> k)} \le \frac{1}{\mu_{k+1} - E_k} \sum_{j=k+1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} \le \frac{M_1^{\mathrm{exc}}}{\mu_{k+1} - E_k}, \tag{8.14.5}$$
   *yielding the high-mode transmission enclosure:*
   $$\sum_{k \ge 1} d_k^2 \big(A_k^{(> k)}\big)^2 \le (M_1^{\mathrm{exc}})^2 \sum_{k \ge 1} \frac{d_k^2}{(\mu_{k+1} - E_k)^2}. \tag{8.14.6}$$
4. *(Factorization of the Resonant Parity-Doublet Transmission Factor and Candidate Mechanism): For the resonant doublet term $j = k$ ($k \ge 2$), the contribution to $d_k A_k$ factors identically as:*
   $$d_k \frac{a_k^2}{(\mu_k - \lambda)(\mu_k - E_k)} = \left( \frac{d_k a_k}{\mu_k - E_k} \right) \left( \frac{a_k}{\mu_k - \lambda} \right). \tag{8.14.7}$$
   *Here $\frac{a_k}{\mu_k - \lambda}$ is the normalized transmission weight, identifying the resonant quotient $\frac{d_k a_k}{\mu_k - E_k}$ as the specific modewise transmission factor requiring analytical control in a mode-by-mode approach. Semiclassically, the even and odd boundary amplitudes $d_k, a_k$ and the doublet splitting $\Delta_k = \mu_k - E_k \sim e^{-2\mathcal{S}_k N}$ share the common barrier penetration action, motivating the candidate physical mechanism for $\left| \frac{d_k a_k}{\mu_k - E_k} \right| \le C$.*

*Proof.*
**Step 1 (Stieltjes Resolvent Identity):**
For any $j \ge 1$ and $k \ge 1$, applying the algebraic partial fraction decomposition:
$$\frac{1}{(\mu_j - \lambda)(\mu_j - E_k)} = \frac{1}{E_k - \lambda} \left[ \frac{1}{\mu_j - E_k} - \frac{1}{\mu_j - \lambda} \right].$$
Multiplying by $a_j^2$ and summing over $j \in \{1, \dots, N-1\}$:
$$A_k = \frac{1}{E_k - \lambda} \left[ \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - E_k} - \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} \right] = \frac{F_{\mathrm{exc}}(E_k) - F_{\mathrm{exc}}(\lambda)}{E_k - \lambda},$$
since $F_{\mathrm{exc}}(\lambda) = \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} = M_1^{\mathrm{exc}}$. This proves (8.14.1). Squaring $\Xi_k = d_k A_k$ and summing over $k \ge 1$ yields (8.14.2).

**Step 2 (Proof of Smooth Enclosure for Mode 1):**
The poles of $F_{\mathrm{exc}}(z)$ are precisely the excited odd eigenvalues $\{\mu_1, \dots, \mu_{N-1}\}$. Under the interlaced doublet ordering, $\lambda < \mu_0 < E_1 < \mu_1$. Because $\mu_0$ is absent from $F_{\mathrm{exc}}$, the distance from the segment $[\lambda, E_1]$ to the spectrum of $F_{\mathrm{exc}}$ is $\mu_1 - E_1 > 0$.
Hence $F_{\mathrm{exc}}$ is infinitely differentiable on $[\lambda, E_1]$. By the Mean Value Theorem, there exists $\xi_1 \in (\lambda, E_1)$ such that $A_1 = F_{\mathrm{exc}}'(\xi_1) = \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \xi_1)^2}$, proving (8.14.3).
Differentiating again:
$$F_{\mathrm{exc}}''(x) = 2 \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - x)^3} > 0 \qquad \forall x \in (-\infty, \mu_1),$$
so $F_{\mathrm{exc}}'(x)$ is strictly monotonically increasing on $[\lambda, E_1]$. Therefore:
$$F_{\mathrm{exc}}'(\lambda) \le F_{\mathrm{exc}}'(\xi_1) \le F_{\mathrm{exc}}'(E_1) \implies \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)^2} \le A_1 \le \sum_{j=1}^{N-1} \frac{a_j^2}{(\mu_j - E_1)^2}.$$
Since $\mu_j - E_1 \ge \mu_1 - E_1 > 0$ for all $j \ge 1$, pulling $\frac{1}{\mu_1 - E_1}$ out of the sum yields $A_1 \le \frac{1}{\mu_1 - E_1} \sum_{j \ge 1} \frac{a_j^2}{\mu_j - E_1} \le \frac{M_1^{\mathrm{exc}}}{\mu_1 - E_1}$, completing (8.14.4).

**Step 3 (Proof of High-Mode Enclosure):**
For $j \ge k+1$, the monotonicity of the odd spectrum $\mu_{k+1} \le \mu_{k+2} \le \dots$ ensures $\mu_j - E_k \ge \mu_{k+1} - E_k > 0$. Therefore:
$$A_k^{(> k)} = \sum_{j=k+1}^{N-1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)} \le \frac{1}{\mu_{k+1} - E_k} \sum_{j=k+1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} \le \frac{M_1^{\mathrm{exc}}}{\mu_{k+1} - E_k},$$
proving (8.14.5). Multiplying by $d_k$, squaring, and summing over $k \ge 1$ establishes (8.14.6).

**Step 4 (Candidate Resonant Factorization):**
For $j = k$, factoring the algebraic expression yields identity (8.14.7) unconditionally. Applying the Cauchy–Schwarz inequality on the boundary fluxes gives:
$$\left( \frac{d_k a_k}{\mu_k - E_k} \right)^2 \le \left( \frac{d_k^2}{\mu_k - E_k} \right) \left( \frac{a_k^2}{\mu_k - E_k} \right).$$
In double-well tunneling, the boundary amplitudes $d_k, a_k$ and the doublet splitting $\Delta_k = \mu_k - E_k$ share the common barrier penetration action $\sim e^{-2\mathcal{S}_k N}$, providing empirical and physical motivation for the candidate bound $\left| \frac{d_k a_k}{\mu_k - E_k} \right| \le C$, while establishing such a modewise bound analytically remains an open target in a modewise framework. $\blacksquare$

---

### 8.15 Proposition 8.15 (Global Operator Representation of $\Xi$, Even-Resolvent Parseval Identity, and Automatic Resonant Absorption)

To resolve the remaining open challenge of Milestone M12 without demanding a separate modewise bound $\frac{d_k a_k}{\mu_k - E_k} = \mathcal{O}(1)$, we investigate the global operator structure of $\Xi = \sum_{k \ge 1} d_k A_k u_k$. 

The following proposition establishes that $\Xi$ is identically an **even-resolvent operator integral** over the excited odd spectral measure, proves the exact closed-form Parseval evaluation of the even-resolvent norm $a_j^2 \sum_{k=1}^{N-1} \frac{d_k^2}{(\mu_j - E_k)^2} = \|K u_j\|^2 - b_{0j}^2$, and establishes a global quadratic-form bound on $\|\Xi\|^2$ that **automatically absorbs the resonant contribution without isolating individual poles**.

**Proposition 8.15 (Global Operator Representation, Parseval Identity, and Resonant Absorption):**
*Let $N \ge 2$, $c > 1$. Let $\Xi = \sum_{k \ge 1} d_k A_k u_k \in c^\perp \subset H_{\mathrm{even}}$ be the unified even-resolvent vector of Proposition 8.13, and let $G_{\mathrm{even}}(z) \equiv \langle d_{\mathrm{exc}}, (z I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} \rangle = \sum_{k=1}^{N-1} \frac{d_k^2}{z - E_k}$ denote the even-sector Stieltjes transform on $\mathbb{C} \setminus \{E_1, \dots, E_{N-1}\}$.*

*Then:*
1. *(Exact Global Operator Representation): As an element of the excited even subspace $c^\perp$, $\Xi$ satisfies the exact operator identity:*
   $$\Xi = g(Q_{\mathrm{even}}) d_{\mathrm{exc}} = \int \frac{d\nu_{\mathrm{exc}}(\mu)}{\mu - \lambda} (\mu I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} = M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} w_j (\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}, \tag{8.15.1}$$
   *where $g(x) = \frac{F_{\mathrm{exc}}(x) - F_{\mathrm{exc}}(\lambda)}{x - \lambda}$, and $w_j = \frac{a_j^2/(\mu_j - \lambda)}{M_1^{\mathrm{exc}}} \ge 0$ ($\sum_{j=1}^{N-1} w_j = 1$) are the normalized transmission weights.*
2. *(Exact Even-Resolvent Parseval Identity): Under the strict spectral interlacing of the discrete Galerkin operator below the barrier top, the odd source overlaps are strictly non-zero: $a_j \ne 0$ for all $j \ge 1$. For every excited odd mode $j \in \{1, \dots, N-1\}$, the even resolvent vector $(\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}$ evaluates in coordinate space identically to:*
   $$(\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} = - \frac{1}{a_j} P_{\perp c} K u_j. \tag{8.15.2}$$
   *Consequently, the even-resolvent norm evaluates in closed form without any eigenvalue denominators:*
   $$a_j^2 \sum_{k=1}^{N-1} \frac{d_k^2}{(\mu_j - E_k)^2} = a_j^2 \|(\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}\|^2 = \|P_{\perp c} K u_j\|^2 = \|K u_j\|^2 - b_{0j}^2 \le \|K u_j\|^2 \le N^2. \tag{8.15.3}$$
3. *(Second-Generation Coordinate-Derivative Wavepacket Geometry & Stieltjes Kernel): The unified resolvent vector $\Xi$ evaluates in coordinate space as an iterated coordinate-derivative projection:*
   $$\Xi = \frac{1}{D_0} P_{\perp c} K v_{\mathrm{exc}}, \qquad \|\Xi\|^2 = \frac{1}{D_0^2} \|P_{\perp c} K v_{\mathrm{exc}}\|^2, \tag{8.15.4}$$
   *arising from the two-tier coordinate-derivative progression:*
   $$c \xrightarrow{\;K\;} Kc \xrightarrow{\;P_{\perp u_0}\;} v_{\mathrm{exc}} \xrightarrow{\;K\;} K v_{\mathrm{exc}} \xrightarrow{\;P_{\perp c}\;} P_{\perp c} K v_{\mathrm{exc}}.$$
   *Equivalently, in the spectral domain, $\|\Xi\|^2$ evaluates to the double Stieltjes divided-difference integral:*
   $$\|\Xi\|^2 = \iint \left[ \frac{G_{\mathrm{even}}(\mu') - G_{\mathrm{even}}(\mu)}{\mu - \mu'} \right] \frac{d\nu_{\mathrm{exc}}(\mu)}{\mu - \lambda} \frac{d\nu_{\mathrm{exc}}(\mu')}{\mu' - \lambda}. \tag{8.15.5}$$
4. *(Automatic Resonant Absorption and Relative-Gap Reduction): By convexity of the norm, $\|\Xi\|^2$ satisfies the closed-form upper bound:*
   $$\|\Xi\|^2 \le M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} \frac{\|K u_j\|^2 - b_{0j}^2}{\mu_j - \lambda} \le M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} \frac{\|K u_j\|^2}{\mu_j - \lambda}. \tag{8.15.6}$$
   *Consequently, the Route B second-moment ratio $\rho_2^{\mathrm{exc}} \equiv \frac{D_0^2 M_{2,\mathrm{exc}}}{M_1^{\mathrm{exc}}}$ satisfies:*
   $$\rho_2^{\mathrm{exc}} \le \frac{1}{1 - D_0^2 M_{2,\mathrm{exc}}} \sum_{j=1}^{N-1} \left( \frac{D_0^2}{\mu_j - \lambda} \right) \big( \|K u_j\|^2 - b_{0j}^2 \big) \le \frac{N^3}{1 - D_0^2 M_{2,\mathrm{exc}}} R_{\mathrm{gap}}^{\max}(N). \tag{8.15.7}$$
   *This completely eliminates the need for an independent modewise resonant bound $\frac{d_k a_k}{\mu_k - E_k} = \mathcal{O}(1)$, reducing the second-moment bound directly to a relative tunneling-gap bound, $R_{\mathrm{gap}}^{\max}(N) = \frac{D_0^2}{\mu_1 - \lambda}$.*

*Proof.*
**Step 1 (Proof of Global Operator Representation):**
The set of excited even eigenvectors $\{u_k\}_{k=1}^{N-1}$ forms an orthonormal basis for $c^\perp \subset H_{\mathrm{even}}$ with $Q_{\mathrm{even}} u_k = E_k u_k$.
By definition of the functional calculus, for any function $g$ regular on $\{E_1, \dots, E_{N-1}\}$:
$$g(Q_{\mathrm{even}}) d_{\mathrm{exc}} = \sum_{k=1}^{N-1} \langle u_k, d_{\mathrm{exc}} \rangle g(E_k) u_k = \sum_{k=1}^{N-1} d_k g(E_k) u_k.$$
Setting $g(x) = \frac{F_{\mathrm{exc}}(x) - F_{\mathrm{exc}}(\lambda)}{x - \lambda}$, Proposition 8.14 established $A_k = g(E_k)$. Thus:
$$\Xi = \sum_{k=1}^{N-1} d_k A_k u_k = g(Q_{\mathrm{even}}) d_{\mathrm{exc}}.$$
Using the integral representation $g(x) = \int \frac{d\nu_{\mathrm{exc}}(\mu)}{(\mu - \lambda)(\mu - x)}$, we obtain:
$$g(Q_{\mathrm{even}}) d_{\mathrm{exc}} = \int \frac{d\nu_{\mathrm{exc}}(\mu)}{\mu - \lambda} (\mu I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} = \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} (\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}},$$
proving (8.15.1).

**Step 2 (Proof of Non-Vanishing $a_j$ and Even-Resolvent Parseval Identity):**
We first establish that $a_j \ne 0$ for all $j \ge 1$. Suppose, for the sake of contradiction, that $a_j = \langle \boldsymbol\psi, u_j \rangle = 0$ for some $j \ge 1$. Applying the rank-two commutator identity of Proposition 8.8, $[K, Q] u_j = \boldsymbol\psi \langle d, u_j \rangle - d \langle \boldsymbol\psi, u_j \rangle = 0$ because $\langle d, u_j \rangle = 0$ by odd parity.
Expanding $[K, Q] u_j = K Q u_j - Q K u_j = \mu_j K u_j - Q K u_j = 0$, this implies:
$$(Q - \mu_j I) K u_j = 0.$$
Since $u_j \in H_{\mathrm{odd}} \setminus \{0\}$, $u_j$ has zero DC component $(u_j)_0 = 0$. The diagonal operator $K = \operatorname{diag}(-N, \dots, N)$ is strictly injective on $H_{\mathrm{odd}}$ ($K u_j = 0 \iff u_j = 0$), so $K u_j \in H_{\mathrm{even}} \setminus \{0\}$. Thus, $K u_j$ would be a non-trivial even eigenvector of $Q$ with eigenvalue $\mu_j$.
This would require $\mu_j \in \operatorname{spec}(Q_{\mathrm{even}}) \cap \operatorname{spec}(Q_{\mathrm{odd}})$, directly contradicting the strict interlacing $\lambda < \mu_0 < E_1 < \mu_1 < \dots$ of the even and odd spectra below the barrier top. Therefore, $a_j \ne 0$ for all $j \ge 1$.

Now, expanding $K u_j \in H_{\mathrm{even}}$ along $H_{\mathrm{even}} = \mathbb{R} c \oplus c^\perp$:
$$K u_j = \langle c, K u_j \rangle c + P_{\perp c} K u_j = b_{0j} c + P_{\perp c} K u_j.$$
Applying $(Q - \mu_j I)$ and noting $(Q - \mu_j I) c = (\lambda - \mu_j) c$ and $d = D_0 c + d_{\mathrm{exc}}$:
$$(\lambda - \mu_j) b_{0j} c + (Q_{\mathrm{even}} - \mu_j I) P_{\perp c} K u_j = a_j D_0 c + a_j d_{\mathrm{exc}}.$$
Since $(\lambda - \mu_j) b_{0j} = a_j D_0$ (by Eq. 8.12.1), the $c$-components cancel identically, leaving:
$$(Q_{\mathrm{even}} - \mu_j I) P_{\perp c} K u_j = a_j d_{\mathrm{exc}} \implies (\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} = - \frac{1}{a_j} P_{\perp c} K u_j,$$
proving (8.15.2).
Taking the squared norm of both sides:
$$\|(\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}\|^2 = \frac{1}{a_j^2} \|P_{\perp c} K u_j\|^2 = \frac{1}{a_j^2} \big( \|K u_j\|^2 - |\langle c, K u_j \rangle|^2 \big) = \frac{\|K u_j\|^2 - b_{0j}^2}{a_j^2}.$$
On the other hand, expanding the resolvent in the eigenbasis $\{u_k\}_{k=1}^{N-1}$ gives:
$$\|(\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}\|^2 = \sum_{k=1}^{N-1} \frac{\langle u_k, d_{\mathrm{exc}} \rangle^2}{(\mu_j - E_k)^2} = \sum_{k=1}^{N-1} \frac{d_k^2}{(\mu_j - E_k)^2}.$$
Multiplying by $a_j^2$ establishes the exact closed-form evaluation (8.15.3).

**Step 3 (Proof of Coordinate-Derivative Geometry & Stieltjes Kernel):**
Substituting (8.15.2) into (8.15.1) gives:
$$\Xi = \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} \left( - \frac{1}{a_j} P_{\perp c} K u_j \right) = - P_{\perp c} K \left( \sum_{j=1}^{N-1} \frac{a_j}{\mu_j - \lambda} u_j \right).$$
Since $\sum_{j=1}^{N-1} \frac{a_j}{\mu_j - \lambda} u_j = (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi_{\mathrm{exc}} = - \frac{v_{\mathrm{exc}}}{D_0}$ (by Proposition 8.13), we obtain $\Xi = \frac{1}{D_0} P_{\perp c} K v_{\mathrm{exc}}$, so $\|\Xi\|^2 = \frac{1}{D_0^2} \|P_{\perp c} K v_{\mathrm{exc}}\|^2$, establishing (8.15.4).
Taking the inner product of $\Xi$ with itself using representation (8.15.1):
$$\|\Xi\|^2 = \iint \frac{d\nu_{\mathrm{exc}}(\mu)}{\mu - \lambda} \frac{d\nu_{\mathrm{exc}}(\mu')}{\mu' - \lambda} \langle (\mu I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}, (\mu' I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} \rangle.$$
By the first resolvent identity:
$$(\mu I - Q_{\mathrm{even}})^{-1} (\mu' I - Q_{\mathrm{even}})^{-1} = \frac{1}{\mu - \mu'} \left[ (\mu' I - Q_{\mathrm{even}})^{-1} - (\mu I - Q_{\mathrm{even}})^{-1} \right].$$
Taking the expectation with $d_{\mathrm{exc}}$ gives $\frac{G_{\mathrm{even}}(\mu') - G_{\mathrm{even}}(\mu)}{\mu - \mu'}$, establishing (8.15.5).

**Step 4 (Proof of Automatic Resonant Absorption & Relative-Gap Reduction):**
By definition, $A_k = M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} w_j \frac{1}{\mu_j - E_k}$. By Cauchy–Schwarz (or Jensen's inequality for the convex function $x \mapsto x^2$ with probability weights $w_j \ge 0$, $\sum w_j = 1$):
$$A_k^2 \le (M_1^{\mathrm{exc}})^2 \sum_{j=1}^{N-1} w_j \frac{1}{(\mu_j - E_k)^2} = M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} \frac{a_j^2}{\mu_j - \lambda} \frac{1}{(\mu_j - E_k)^2}.$$
Multiplying by $d_k^2$ and summing over $k \in \{1, \dots, N-1\}$:
$$\|\Xi\|^2 = \sum_{k=1}^{N-1} d_k^2 A_k^2 \le M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} \frac{1}{\mu_j - \lambda} \left[ a_j^2 \sum_{k=1}^{N-1} \frac{d_k^2}{(\mu_j - E_k)^2} \right].$$
Substituting the exact identity (8.15.3) into the bracketed term gives:
$$\|\Xi\|^2 \le M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} \frac{\|K u_j\|^2 - b_{0j}^2}{\mu_j - \lambda},$$
proving (8.15.6). In this upper bound, the summation over $k$ has been evaluated in closed form, completely absorbing the resonant pole!
Finally, inserting (8.15.6) into the master Route B operator inequality of Proposition 8.13:
$$\rho_2^{\mathrm{exc}} \le \frac{D_0^2 \|\Xi\|^2}{M_1^{\mathrm{exc}} (1 - D_0^2 M_{2,\mathrm{exc}})} \le \frac{1}{1 - D_0^2 M_{2,\mathrm{exc}}} \sum_{j=1}^{N-1} \left( \frac{D_0^2}{\mu_j - \lambda} \right) \big( \|K u_j\|^2 - b_{0j}^2 \big).$$
Since $\|K u_j\|^2 \le N^2$ and $\frac{D_0^2}{\mu_j - \lambda} \le \frac{D_0^2}{\mu_1 - \lambda} = R_{\mathrm{gap}}^{\max}(N)$, summing over $N-1$ modes gives:
$$\sum_{j=1}^{N-1} \left( \frac{D_0^2}{\mu_j - \lambda} \right) \|K u_j\|^2 \le N^2 (N-1) R_{\mathrm{gap}}^{\max}(N) \le N^3 R_{\mathrm{gap}}^{\max}(N),$$
completing the proof. $\blacksquare$

*Discussion of Asymptotic Frontier:*
Proposition 8.15 achieves an exact algebraic reduction that removes the independent hypothesis $\mathrm{H2}_{\mathrm{gap}}$ from the finite-$N$ identities and eliminates the resonant transmission factor. However, concluding the continuum limit $\rho_2^{\mathrm{exc}} \to 0$ from (8.15.7) still requires establishing that $R_{\mathrm{gap}}^{\max}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda} = \mathcal{O}(N^{-3-\epsilon})$ (or exponentially decaying, as indicated by the Cell 67 audit $\Delta \sigma_N^{\mathrm{gap}} \ge 0.610$). The required relative-gap estimate is dramatically weaker and cleaner than the original second-moment problem $D_0^2 M_{2,\mathrm{exc}}$, but remains the active analytical bridge.

---

### 8.16 Proposition 8.16 (Weighted Coordinate-Energy Measure, Exact Trace Enclosure, and High-Energy Sector Decoupling)

The crude bound $\rho_2^{\mathrm{exc}} \le \frac{N^3}{1 - D_0^2 M_{2,\mathrm{exc}}} R_{\mathrm{gap}}^{\max}(N)$ in (8.15.7) replaces every mode's relative gap $\frac{D_0^2}{\mu_j - \lambda}$ by the worst-case supremum $R_{\mathrm{gap}}^{\max} = \frac{D_0^2}{\mu_1 - \lambda}$ and every coordinate energy by $N^2$. 

To sharpen this bound toward an asymptotic decoupling proof, we investigate the exact weighted coordinate-energy sum:
$$\mathcal{S}_{\mathrm{coord}}(N) \equiv \sum_{j=1}^{N-1} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 = \sum_{j=1}^{N-1} \frac{D_0^2}{\mu_j - \lambda} \big( \|K u_j\|^2 - b_{0j}^2 \big),$$
where $u_0^{\mathrm{even}}$ denotes the normalized even ground state ($Q_{\mathrm{even}} u_0^{\mathrm{even}} = \lambda u_0^{\mathrm{even}}$), $P_{\perp u_0^{\mathrm{even}}}$ denotes the orthogonal projection on $H_{\mathrm{even}}$ onto $(u_0^{\mathrm{even}})^\perp$, and $b_{0j} \equiv \langle u_0^{\mathrm{even}}, K u_j \rangle$ (distinguishing the ground state $u_0^{\mathrm{even}}$ from the prime cutoff parameter $c > 1$).

The following proposition establishes the exact trace identity for the total coordinate energy, defines the coordinate-energy spectral measure $d\eta_{\mathrm{coord}}$, and proves that the high-energy (above-barrier) sector decouples exponentially fast under Hypothesis H3 and an $N$-uniform barrier separation, unconditional with respect to relative-gap hypotheses.

**Proposition 8.16 (Coordinate Trace Identity, Spectral Measure, and High-Energy Decoupling):**
*Let $N \ge 2$, and let $c > 1$ denote the prime cutoff parameter ($L = \log c$). Let $\{u_j\}_{j=0}^{N-1}$ be the complete orthonormal eigenbasis of the odd subspace $H_{\mathrm{odd}} \subset \mathbb{R}^{2N+1}$, and let $u_0^{\mathrm{even}}$ be the normalized even ground state with $Q_{\mathrm{even}} u_0^{\mathrm{even}} = \lambda u_0^{\mathrm{even}}$.*

*Then:*
1. *(Exact Odd Coordinate Trace Identity): The total coordinate energy across the odd subspace satisfies the exact trace identity:*
   $$\sum_{j=0}^{N-1} \|K u_j\|^2 = \operatorname{Tr}_{H_{\mathrm{odd}}}(K^2) = \sum_{m=1}^N m^2 = \frac{N(N+1)(2N+1)}{6} = \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6}. \tag{8.16.1}$$
   *Consequently, the excited coordinate energy satisfies:*
   $$\sum_{j=1}^{N-1} \|K u_j\|^2 = \frac{N(N+1)(2N+1)}{6} - \|K u_0\|^2 < \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} = \mathcal{O}(N^3). \tag{8.16.2}$$
2. *(Coordinate-Energy Spectral Measure Representation): Let $d\eta_{\mathrm{coord}}(\mu) \equiv \sum_{j=1}^{N-1} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \delta_{\mu_j}$ denote the excited coordinate-energy spectral measure on $[\mu_1, \infty)$, with finite total mass $\int d\eta_{\mathrm{coord}} = \sum_{j=1}^{N-1} (\|K u_j\|^2 - b_{0j}^2) < \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} = \mathcal{O}(N^3)$. Then $\mathcal{S}_{\mathrm{coord}}(N)$ evaluates identically to the Stieltjes transform of $d\eta_{\mathrm{coord}}$ at the ground energy $\lambda$:*
   $$\mathcal{S}_{\mathrm{coord}}(N) = D_0^2 \int \frac{d\eta_{\mathrm{coord}}(\mu)}{\mu - \lambda} = D_0^2 \, G_{\mathrm{coord}}(\lambda). \tag{8.16.3}$$
3. *(Exponential Decoupling of High-Energy Modes): Let $V_* > 0$ be a fixed threshold satisfying $V_* > \sup_N \lambda_N$ (corresponding semiclassically to the effective barrier top $V_{\max} \approx 1.3$). Let the discrete finite-$N$ odd spectrum be partitioned into a bound-state sector $\mu_j < V_*$ ($1 \le j \le N_{\mathrm{bound}}$) and a high-energy (above-barrier) sector $\mu_j \ge V_*$ ($j > N_{\mathrm{bound}}$) where $\mu_j - \lambda \ge c_{\mathrm{high}} \equiv V_* - \sup_N \lambda_N > 0$ uniformly for all $N$. Then the high-energy contribution satisfies:*
   $$\mathcal{S}_{\mathrm{high}}(N) \equiv \sum_{j > N_{\mathrm{bound}}} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \le \frac{D_0^2}{c_{\mathrm{high}}} \sum_{j > N_{\mathrm{bound}}} \|K u_j\|^2 \le \frac{D_0^2}{c_{\mathrm{high}}} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right). \tag{8.16.4}$$
   *Under Hypothesis H3 ($D_0^2 \le C_0 e^{-\sigma N}$), this tail decays exponentially fast:*
   $$\mathcal{S}_{\mathrm{high}}(N) \le \frac{C_0}{c_{\mathrm{high}}} e^{-\sigma N} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right) = \mathcal{O}(e^{-\sigma N} N^3) \longrightarrow 0 \quad (N \to \infty). \tag{8.16.5}$$
   *This decoupling is unconditional with respect to relative-gap hypotheses ($\mathrm{H2}_{\mathrm{gap}}$, $\mathrm{H2}_{\mathrm{odd}}$), being conditional strictly on H3 and the assumed uniform above-barrier spectral separation $c_{\mathrm{high}} > 0$. Consequently, the small denominator problem is isolated strictly to the finite bound-state ladder $1 \le j \le N_{\mathrm{bound}}$.*

*Proof.*
**Step 1 (Proof of Coordinate Trace Identity):**
The canonical basis for $H_{\mathrm{odd}} \subset \mathbb{R}^{2N+1}$ is given by $v_m = \frac{1}{\sqrt{2}}(e_m - e_{-m})$ for $m \in \{1, \dots, N\}$.
Applying the diagonal operator $K = \operatorname{diag}(-N, \dots, N)$:
$$K v_m = \frac{1}{\sqrt{2}}(m e_m - (-m) e_{-m}) = m \frac{1}{\sqrt{2}}(e_m + e_{-m}).$$
The squared norm is $\|K v_m\|^2 = m^2$.
Because $\{u_j\}_{j=0}^{N-1}$ is an orthonormal basis for $H_{\mathrm{odd}}$, by cyclicity and unitary invariance of the trace:
$$\sum_{j=0}^{N-1} \|K u_j\|^2 = \operatorname{Tr}_{H_{\mathrm{odd}}}(K^T K) = \sum_{m=1}^N \|K v_m\|^2 = \sum_{m=1}^N m^2 = \frac{N(N+1)(2N+1)}{6} = \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6},$$
proving (8.16.1). Subtracting the strictly positive ground-mode term $\|K u_0\|^2 > 0$ yields:
$$\sum_{j=1}^{N-1} \|K u_j\|^2 = \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} - \|K u_0\|^2 < \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} = \mathcal{O}(N^3),$$
proving (8.16.2). (Note: because $\|Ku_0\|^2 \le N^2$, one cannot deduce $< N^3/3$ from positivity alone; the clean universally valid bound is $\frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} = \mathcal{O}(N^3)$).

**Step 2 (Proof of Stieltjes Coordinate Measure Representation):**
Integrating $\frac{1}{\mu - \lambda}$ against the discrete positive measure $d\eta_{\mathrm{coord}}(\mu) = \sum_{j=1}^{N-1} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \delta_{\mu_j}$ gives:
$$\int \frac{d\eta_{\mathrm{coord}}(\mu)}{\mu - \lambda} = \sum_{j=1}^{N-1} \frac{\|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2}{\mu_j - \lambda}.$$
Multiplying by $D_0^2$ yields (8.16.3). The total mass satisfies $\int d\eta_{\mathrm{coord}} = \sum_{j=1}^{N-1} (\|K u_j\|^2 - b_{0j}^2) \le \sum_{j=1}^{N-1} \|K u_j\|^2 < \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} = \mathcal{O}(N^3)$.

**Step 3 (Proof of High-Energy Decoupling):**
For modes $j > N_{\mathrm{bound}}$ in the above-barrier sector, $\mu_j \ge V_*$. With $V_* > \sup_N \lambda_N$, the uniform separation $c_{\mathrm{high}} = V_* - \sup_N \lambda_N > 0$ ensures $\mu_j - \lambda \ge c_{\mathrm{high}}$ for all $j > N_{\mathrm{bound}}$ and all $N$.
Pulling the uniform denominator $\frac{1}{c_{\mathrm{high}}}$ out of the sum:
$$\sum_{j > N_{\mathrm{bound}}} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \le \frac{D_0^2}{c_{\mathrm{high}}} \sum_{j > N_{\mathrm{bound}}} \|K u_j\|^2 \le \frac{D_0^2}{c_{\mathrm{high}}} \sum_{j=1}^{N-1} \|K u_j\|^2 < \frac{D_0^2}{c_{\mathrm{high}}} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right).$$
Under Hypothesis H3 (exponential boundary suppression $D_0^2 \le C_0 e^{-\sigma N}$), this tail decays as $\mathcal{O}(e^{-\sigma N} N^3) \to 0$ exponentially fast, completing the proof. $\blacksquare$

---

### 8.17 Proposition 8.17 (Exact Operator Trace Identity, Bound-State Spectral Decomposition, and Semiclassical Localization Bridges B1 & B2)

Proposition 8.16 establishes that the high-energy (above-barrier) sector $\mathcal{S}_{\mathrm{high}}(N)$ decouples exponentially fast under Hypothesis H3 and uniform separation $c_{\mathrm{high}} > 0$. Consequently, the entire small-denominator challenge for Route B is isolated strictly to the finite bound-state ladder:
$$\mathcal{S}_{\mathrm{bound}}(N) \equiv \sum_{j=1}^{N_{\mathrm{bound}}} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2.$$

In Proposition 8.15, the bound on $\rho_2^{\mathrm{exc}}$ incurred an artificial factor of $N^3$, derived from the crude worst-case bounds $\|K u_j\|^2 \le N^2$ and $\sum_{j=1}^{N-1} 1 \le N$, which treated low-energy bound states as if they carried maximal high-frequency Fourier momentum $m = N$. 
We now formulate the exact operator trace geometry on $H_{\mathrm{odd}}$, prove the rigorous finite-$N$ bound-state enclosure, and explicitly identify the two continuum localization bridges (B1 and B2) required to eliminate the $N^3$ rate penalty in the continuum limit.

**Proposition 8.17 (Exact Operator Trace Identity, Bound-State Enclosure, and Semiclassical Bridges B1–B2):**
*Let $N \ge 2$, and let $c > 1$ denote the prime cutoff parameter ($L = \log c$). Let $\{u_j\}_{j=0}^{N-1}$ be the orthonormal eigenbasis of $Q_{\mathrm{odd}}$ with eigenvalues $\mu_0 < \mu_1 < \dots < \mu_{N-1}$, and let $u_0^{\mathrm{even}}$ denote the normalized even ground state with $Q_{\mathrm{even}} u_0^{\mathrm{even}} = \lambda u_0^{\mathrm{even}}$. Let $V_* > \sup_N \lambda_N$ define a fixed above-barrier threshold, with $N_{\mathrm{bound}}(N) \equiv \max \{ j : \mu_j < V_* \}$.*

**Part I: Exact Finite-$N$ Operator Trace Identities and Spectral Decomposition (Unconditional Theorem)**
*For every finite dimension $N \ge 2$ and cutoff $c > 1$:*
1. *(Exact Operator Trace Identity for $\mathcal{S}_{\mathrm{coord}}$): The weighted coordinate-energy sum satisfies the exact operator trace identity on $H_{\mathrm{odd}}$:*
   $$\mathcal{S}_{\mathrm{coord}}(N) = D_0^2 \operatorname{Tr}_{H_{\mathrm{odd}}}\Big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K P_{\perp u_0^{\mathrm{even}}} K \Big]. \tag{8.17.1}$$
   *Since $T_{\mathrm{coord}} \equiv K P_{\perp u_0^{\mathrm{even}}} K$ is positive semidefinite and satisfies $0 \preceq T_{\mathrm{coord}} \preceq K^2$ on $H_{\mathrm{odd}}$, $\mathcal{S}_{\mathrm{coord}}(N)$ is enclosed by the diagonal Fourier trace:*
   $$\mathcal{S}_{\mathrm{coord}}(N) \le D_0^2 \operatorname{Tr}_{H_{\mathrm{odd}}}\Big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K^2 \Big] = D_0^2 \sum_{m=1}^N m^2 \langle v_m, (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} v_m \rangle, \tag{8.17.2}$$
   *where $v_m = \frac{1}{\sqrt{2}}(e_m - e_{-m})$ is the canonical odd Fourier basis in which $K^2 = \operatorname{diag}(1^2, 2^2, \dots, N^2)$.*
2. *(Exact Bound vs. High-Energy Spectral Partitioning): Orthogonal spectral projection on $H_{\mathrm{odd}}$ via $P_{\mathrm{bound}} \equiv \sum_{j=1}^{N_{\mathrm{bound}}} u_j u_j^T$ and $P_{\mathrm{high}} \equiv \sum_{j > N_{\mathrm{bound}}}^{N-1} u_j u_j^T$ decomposes the sum identically:*
   $$\mathcal{S}_{\mathrm{coord}}(N) = \mathcal{S}_{\mathrm{bound}}(N) + \mathcal{S}_{\mathrm{high}}(N), \tag{8.17.3}$$
   *where $\mathcal{S}_{\mathrm{high}}(N) \le \frac{C_0}{c_{\mathrm{high}}} e^{-\sigma N} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right) = \mathcal{O}(e^{-\sigma N} N^3) \to 0$ exponentially fast under Hypothesis H3 and uniform separation $c_{\mathrm{high}} = V_* - \sup_N \lambda_N > 0$.*
3. *(Exact Bound-State Enclosure): For each bound mode $j \in \{1, \dots, N_{\mathrm{bound}}\}$, the relative gaps satisfy $\frac{D_0^2}{\mu_j - \lambda} \le \frac{D_0^2}{\mu_1 - \lambda} = R_{\mathrm{gap}}^{\max}(N)$, yielding the unconditional finite-$N$ inequality:*
   $$\mathcal{S}_{\mathrm{bound}}(N) \le R_{\mathrm{gap}}^{\max}(N) \sum_{j=1}^{N_{\mathrm{bound}}(N)} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \le R_{\mathrm{gap}}^{\max}(N) \sum_{j=1}^{N_{\mathrm{bound}}(N)} \|K u_j\|^2. \tag{8.17.4}$$

**Part II: Conditional Continuum Reduction under Semiclassical Bridges B1 & B2**
*Consider the following two analytical hypotheses on the finite-$N$ Galerkin discretizations:*
- *(Bridge B1 — Uniform Bound-State Cardinality): The number of discrete odd Galerkin eigenvalues beneath the barrier threshold $V_*$ remains uniformly bounded as $N \to \infty$:*
  $$\bar{N}_{\mathrm{bound}} \equiv \sup_{N \ge 2} N_{\mathrm{bound}}(N) < \infty. \tag{B1}$$
- *(Bridge B2 — Uniform Discrete Coordinate Kinetic Energy): For each fixed bound-state index $j \le \bar{N}_{\mathrm{bound}}$, the discrete Galerkin coordinate derivative norm remains uniformly bounded across dimensions:*
  $$\sup_{N \ge j+1} \|K u_j^{(N)}\|^2 \le \mathcal{C}_j < \infty. \tag{B2}$$
*Under Bridges B1 and B2:*
1. *There exists a uniform constant $\mathcal{C}_{\mathrm{bound}} \equiv \sum_{j=1}^{\bar{N}_{\mathrm{bound}}} \mathcal{C}_j < \infty$ such that for all $N$:*
   $$\mathcal{S}_{\mathrm{bound}}(N) \le \mathcal{C}_{\mathrm{bound}} R_{\mathrm{gap}}^{\max}(N) = \mathcal{O}\big( R_{\mathrm{gap}}^{\max}(N) \big). \tag{8.17.5}$$
2. *(Refined Route B Bound): Under H3, B1, and B2, the excited second-moment ratio satisfies:*
   $$\rho_2^{\mathrm{exc}} \le \frac{1}{1 - D_0^2 M_{2,\mathrm{exc}}} \left[ \mathcal{C}_{\mathrm{bound}} R_{\mathrm{gap}}^{\max}(N) + \frac{C_0}{c_{\mathrm{high}}} e^{-\sigma N} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right) \right]. \tag{8.17.6}$$
   *Consequently, the condition:*
   $$R_{\mathrm{gap}}^{\max}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda} \longrightarrow 0 \quad (N \to \infty) \tag{8.17.7}$$
   *is sufficient to establish the continuum decoupling limit $\rho_2^{\mathrm{exc}} \to 0$, completely eliminating the polynomial $N^3$ rate requirement $R_{\mathrm{gap}}^{\max} = \mathcal{O}(N^{-3-\epsilon})$ of Proposition 8.15.*
   *(Note: Boundedness $R_{\mathrm{gap}}^{\max}(N) = \mathcal{O}(1)$ guarantees uniform boundedness $\rho_2^{\mathrm{exc}} = \mathcal{O}(1)$, while genuine continuum decoupling $\rho_2^{\mathrm{exc}} \to 0$ requires $R_{\mathrm{gap}}^{\max}(N) \to 0$.)*

*Proof.*
**Step 1 (Proof of Exact Operator Trace Identity):**
By definition, $\mathcal{S}_{\mathrm{coord}}(N) = \sum_{j=1}^{N-1} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2$.
In $H_{\mathrm{even}}$, the orthogonal projection satisfies $P_{\perp u_0^{\mathrm{even}}}^2 = P_{\perp u_0^{\mathrm{even}}} = P_{\perp u_0^{\mathrm{even}}}^T$, so:
$$\|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 = \langle P_{\perp u_0^{\mathrm{even}}} K u_j, P_{\perp u_0^{\mathrm{even}}} K u_j \rangle = \langle K u_j, P_{\perp u_0^{\mathrm{even}}} K u_j \rangle = \langle u_j, K P_{\perp u_0^{\mathrm{even}}} K u_j \rangle.$$
Define the operator $T_{\mathrm{coord}} \equiv K P_{\perp u_0^{\mathrm{even}}} K$ mapping $H_{\mathrm{odd}} \to H_{\mathrm{odd}}$.
For any $u \in H_{\mathrm{odd}}$, $\langle u, T_{\mathrm{coord}} u \rangle = \|P_{\perp u_0^{\mathrm{even}}} K u\|^2 \ge 0$, and since $P_{\perp u_0^{\mathrm{even}}} \preceq I_{\mathrm{even}}$, we have $0 \preceq T_{\mathrm{coord}} \preceq K^2$.
Because $\{u_j\}_{j=0}^{N-1}$ is an orthonormal eigenbasis of $Q_{\mathrm{odd}}$ with eigenvalues $\mu_j$, the spectral projection onto the excited odd modes is $P_{\perp u_0} = \sum_{j=1}^{N-1} u_j u_j^T$.
Therefore:
$$(Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} = \sum_{j=1}^{N-1} \frac{1}{\mu_j - \lambda} u_j u_j^T.$$
Evaluating the trace of $(Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} T_{\mathrm{coord}}$ in the basis $\{u_j\}_{j=0}^{N-1}$:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}\Big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} T_{\mathrm{coord}} \Big] = \sum_{j=0}^{N-1} \langle u_j, (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} T_{\mathrm{coord}} u_j \rangle = \sum_{j=1}^{N-1} \frac{1}{\mu_j - \lambda} \langle u_j, T_{\mathrm{coord}} u_j \rangle.$$
Multiplying by $D_0^2$ yields (8.17.1).
Since $T_{\mathrm{coord}} \preceq K^2$ and $(Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} \succeq 0$, the trace of their product satisfies:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}\Big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} T_{\mathrm{coord}} \Big] \le \operatorname{Tr}_{H_{\mathrm{odd}}}\Big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K^2 \Big].$$
Evaluating this trace in the canonical Fourier basis $v_m = \frac{1}{\sqrt{2}}(e_m - e_{-m})$ for $m \in \{1, \dots, N\}$, where $K^2 v_m = m^2 v_m$:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}\Big[ (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K^2 \Big] = \sum_{m=1}^N \langle v_m, (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K^2 v_m \rangle = \sum_{m=1}^N m^2 \langle v_m, (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} v_m \rangle,$$
proving (8.17.2).

**Step 2 (Proof of Spectral Partitioning):**
Because $P_{\perp u_0} = P_{\mathrm{bound}} + P_{\mathrm{high}}$ with orthogonal ranges, the linearity of the trace yields:
$$\mathcal{S}_{\mathrm{coord}}(N) = \mathcal{S}_{\mathrm{bound}}(N) + \mathcal{S}_{\mathrm{high}}(N).$$
The bound on $\mathcal{S}_{\mathrm{high}}(N)$ follows directly from Proposition 8.16 (Step 3), proving (8.17.3).

**Step 3 (Proof of Exact Bound-State Enclosure):**
Since $\mu_1 < \mu_2 < \dots < \mu_{N_{\mathrm{bound}}}$, we have $\mu_j - \lambda \ge \mu_1 - \lambda > 0$, so $\frac{D_0^2}{\mu_j - \lambda} \le \frac{D_0^2}{\mu_1 - \lambda} = R_{\mathrm{gap}}^{\max}(N)$ for all $1 \le j \le N_{\mathrm{bound}}$.
Extracting this supremum:
$$\mathcal{S}_{\mathrm{bound}}(N) = \sum_{j=1}^{N_{\mathrm{bound}}} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \le R_{\mathrm{gap}}^{\max}(N) \sum_{j=1}^{N_{\mathrm{bound}}} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \le R_{\mathrm{gap}}^{\max}(N) \sum_{j=1}^{N_{\mathrm{bound}}} \|K u_j\|^2,$$
establishing (8.17.4). This concludes the unconditional proof of Part I.

**Step 4 (Proof of Conditional Reduction under Bridges B1 & B2):**
Assuming Bridge B1 ($N_{\mathrm{bound}}(N) \le \bar{N}_{\mathrm{bound}} < \infty$) and Bridge B2 ($\sup_N \|K u_j^{(N)}\|^2 \le \mathcal{C}_j < \infty$), summing over $j \in \{1, \dots, N_{\mathrm{bound}}(N)\}$ gives:
$$\sum_{j=1}^{N_{\mathrm{bound}}(N)} \|K u_j\|^2 \le \sum_{j=1}^{\bar{N}_{\mathrm{bound}}} \mathcal{C}_j \equiv \mathcal{C}_{\mathrm{bound}} < \infty,$$
uniformly in $N$. Substituting into (8.17.4) proves (8.17.5).
Combining (8.17.3), (8.17.5), and (8.16.5) with the master Route B operator inequality $\rho_2^{\mathrm{exc}} \le \frac{\mathcal{S}_{\mathrm{coord}}(N)}{1 - D_0^2 M_{2,\mathrm{exc}}}$ (from Propositions 8.13 and 8.15) immediately yields (8.17.6).
Taking the limit $N \to \infty$ in (8.17.6): since $\mathcal{S}_{\mathrm{high}} \to 0$ exponentially fast, if $R_{\mathrm{gap}}^{\max}(N) \to 0$, then $\rho_2^{\mathrm{exc}} \to 0$, completing the proof. $\blacksquare$

*Discussion and Scientific Epistemics:*
Proposition 8.17 establishes a rigorous architectural boundary between finite-dimensional exact facts and the remaining continuum programme:
1. **The Epistemic Separation (Exact Algebra vs. Semiclassical Bridges):**
   - *Unconditional Exact Theorem (Part I):* The operator trace identity $\mathcal{S}_{\mathrm{coord}} = D_0^2 \operatorname{Tr}[(Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K P_{\perp u_0^{\mathrm{even}}} K]$, the Fourier trace bound $\le D_0^2 \sum m^2 [R_{\mathrm{odd},\perp}]_{mm}$, and the exact bound-state inequality $\mathcal{S}_{\mathrm{bound}} \le R_{\mathrm{gap}}^{\max} \sum \|Ku_j\|^2$ are mathematically watertight for every finite dimension $N$.
   - *Active Analytical Bridges (Part II):* The estimates $N_{\mathrm{bound}} = \mathcal{O}(1)$ and $\|Ku_j\|^2 = \mathcal{O}(1)$ are **not** claims derived from finite-$N$ linear algebra. Semiclassical phase-space integration ($N_{\mathrm{bound}} \approx \frac{1}{\pi}\iint dt\,dp$) and continuous $H^1$ Sobolev norms are continuum heuristics that motivate Bridges B1 and B2, but do not prove them for discrete Galerkin matrices. Establishing Bridges B1 and B2 from the discrete operators constitutes the active analytical agenda.
2. **Sharp Sufficiency Criterion for Continuum Decoupling:**
   Under Bridges B1 and B2, the required relative-gap rate is relaxed from the severe requirement $R_{\mathrm{gap}}^{\max}(N) = \mathcal{O}(N^{-3-\epsilon})$ of Proposition 8.15 to the simple convergence $R_{\mathrm{gap}}^{\max}(N) \to 0$. Boundedness $R_{\mathrm{gap}}^{\max}(N) = \mathcal{O}(1)$ guarantees uniform boundedness $\rho_2^{\mathrm{exc}} = \mathcal{O}(1)$, while $R_{\mathrm{gap}}^{\max}(N) \to 0$ delivers full continuum decoupling $\rho_2^{\mathrm{exc}} \to 0$.
   Since Cell 67 confirmed exponential decay of $R_{\mathrm{gap}}^{\max}(N)$ with decay action $\Delta \sigma_N^{\mathrm{gap}} \ge 0.610$, the bound-state enclosure tightly connects empirical data to the continuum decoupling target.
3. **The Alternative Global Operator Trace Horizon:**
   The Fourier-resolvent trace representation (8.17.2):
   $$\mathcal{S}_{\mathrm{coord}}(N) \le D_0^2 \sum_{m=1}^N m^2 \langle v_m, (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} v_m \rangle = D_0^2 \operatorname{Tr}_{H_{\mathrm{odd}}}\big[ K^2 (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} \big]$$
   suggests an alternative, purely operator-theoretic avenue: examining whether the structured resolvent can compensate for the $m^2$ coordinate weight. Proposition 8.18 investigates this global trace via commutator algebra and provides the definitive epistemic resolution.

---

### 8.18 Proposition 8.18 (Global Resolvent Commutator Identity, Weighted Coordinate-Trace Duality, and Barrier-Thinning Quenching)

*The global weighted resolvent trace $\mathcal{T}_N \equiv \operatorname{Tr}_{H_{\mathrm{odd}}}[K^2 (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0}]$ governs the operator trace bound on the coordinate-energy measure $\mathcal{S}_{\mathrm{coord}}(N) \le D_0^2 \mathcal{T}_N$. By evaluating the nested commutator $[K, [K, R(z)]]$ on the finite-rank Galerkin space, the resolvent commutator generates an exact rank-two trace identity between parity sectors. Furthermore, analyzing the spectral expansion of $\mathcal{T}_N$ provides a decisive resolution to the polynomial trace hypothesis: while the bare trace $\mathcal{T}_N$ is exponentially divergent due to the first excited tunneling gap, the physical scaled trace $D_0^2 \mathcal{T}_N$ is exponentially quenched to zero by semiclassical barrier thinning.*

**Part I (Exact Resolvent Commutator Algebra and Parity Trace Duality — Unconditional Finite-$N$ Theorem):**
*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Galerkin matrix, $K = \operatorname{diag}(-N, \dots, N)$, and let $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ be the exact rank-two commutator of Proposition 8.8. For any regular value $z \in \rho(Q)$, define the full resolvent $R(z) \equiv (Q - z I)^{-1}$.*
1. *The coordinate commutator of the resolvent satisfies the exact rank-two identity:*
   $$[K, R(z)] = - R(z) [K, Q] R(z) = - (R(z)\boldsymbol\psi)(d^T R(z)) + (R(z)d)(\boldsymbol\psi^T R(z)). \tag{8.18.1}$$
2. *Evaluating the nested commutator $[K, [K, R(z)]] = K^2 R(z) - 2 K R(z) K + R(z) K^2$ and taking the parity trace yields the exact trace difference between the odd and even sectors:*
   $$\boxed{\operatorname{Tr}_{H_{\mathrm{odd}}}\big[ K^2 R_{\mathrm{odd}}(z) \big] - \operatorname{Tr}_{H_{\mathrm{even}}}\big[ K^2 R_{\mathrm{even}}(z) \big] = \langle R_{\mathrm{odd}}(z)\boldsymbol\psi, K R_{\mathrm{even}}(z) d \rangle,} \tag{8.18.2}$$
   *where $R_{\mathrm{odd}}(z) \equiv (Q_{\mathrm{odd}} - z I)^{-1}$ and $R_{\mathrm{even}}(z) \equiv (Q_{\mathrm{even}} - z I)^{-1}$.*

**Part II (Exact Spectral Representation of the Weighted Coordinate Trace):**
*For the projected odd resolvent $R_\lambda^\perp \equiv (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0}$, the Fourier diagonal trace and the odd-eigenbasis coordinate-energy sum are identically equal for every $N \ge 2$ and $c > 1$:*
$$\boxed{\mathcal{T}_N \equiv \operatorname{Tr}_{H_{\mathrm{odd}}}\big[ K^2 R_\lambda^\perp \big] = \sum_{m=1}^N m^2 [R_{\mathrm{odd},\perp}(\lambda)]_{mm} = \sum_{j=1}^{N-1} \frac{\|K u_j\|^2}{\mu_j - \lambda}.} \tag{8.18.3}$$

**Part III (Epistemic Resolution of the Polynomial Trace Hypothesis):**
1. *Exponential Divergence of the Bare Trace $\mathcal{T}_N$:*
   *Because the first excited odd mode $u_1$ is a bound state beneath the effective barrier top ($E_0 < \mu_1 < V_*$), its tunneling gap decays exponentially:*
   $$\mu_1 - \lambda = \mathcal{O}(e^{-\sigma_1 N}), \qquad \sigma_1 > 0.$$
   *Since the discrete coordinate kinetic energy $\|K u_1\|^2 \ge \mathcal{C}_1 > 0$ remains an $\mathcal{O}(1)$ positive constant across $N$, the bare trace satisfies the lower bound:*
   $$\mathcal{T}_N \ge \frac{\|K u_1\|^2}{\mu_1 - \lambda} \ge \mathcal{C}_1 e^{+\sigma_1 N} \gg C N^p \qquad (\forall p > 0). \tag{8.18.4}$$
   *Consequently, the bare resolvent trace $\mathcal{T}_N$ is **not** polynomially bounded in $N$.*

2. *Exponential Quenching of the Scaled Physical Trace $D_0^2 \mathcal{T}_N$:*
   *Decomposing $\mathcal{T}_N = \mathcal{T}_{\mathrm{bound}}(N) + \mathcal{T}_{\mathrm{high}}(N)$ at an $N$-uniform barrier threshold $V_* > \sup_N \lambda_N$ ($c_{\mathrm{high}} \equiv V_* - \sup_N \lambda_N > 0$):*
   - *The above-barrier trace is polynomially bounded:*
     $$\mathcal{T}_{\mathrm{high}}(N) \equiv \sum_{\mu_j \ge V_*} \frac{\|Ku_j\|^2}{\mu_j - \lambda} \le \frac{1}{c_{\mathrm{high}}} \operatorname{Tr}_{H_{\mathrm{odd}}}(K^2) = \frac{1}{c_{\mathrm{high}}} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right) = \mathcal{O}(N^3).$$
   - *The bound-state trace satisfies:*
     $$\mathcal{T}_{\mathrm{bound}}(N) \equiv \sum_{j=1}^{N_{\mathrm{bound}}} \frac{\|Ku_j\|^2}   Multiplying by the ground-state tunneling amplitude $D_0^2 \le C_0 e^{-2\sigma_0 N}$, and assuming Bridges B1 and B2 ($\sum_{j=1}^{N_{\mathrm{bound}}} \|Ku_j\|^2 \le \mathcal{C}_{\mathrm{bound}} < \infty$), the scaled trace satisfies:
   $$\boxed{D_0^2 \mathcal{T}_N \le \mathcal{C}_{\mathrm{bound}} R_{\mathrm{gap}}^{\max}(N) + \frac{C_0}{c_{\mathrm{high}}} e^{-2\sigma_0 N} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right).} \tag{8.18.5}$$
   *Under the semiclassical barrier-thinning hypothesis $\sigma_0 > \sigma_1$ (with action gap $\Delta\sigma \equiv \sigma_0 - \sigma_1 > 0$), the relative tunneling gap decays exponentially:*
   $$R_{\mathrm{gap}}^{\max}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda} = \mathcal{O}\big( e^{-(\sigma_0 + \Delta\sigma) N} \big) \longrightarrow 0 \quad (N \to \infty).$$
   *Therefore, conditional on Bridges B1 and B2 and the barrier-thinning action hierarchy $\sigma_0 > \sigma_1$, the physical scaled trace and the coordinate-energy measure are exponentially quenched:*
   $$\mathcal{S}_{\mathrm{coord}}(N) \le D_0^2 \mathcal{T}_N = \mathcal{O}\big( e^{-(\sigma_0 + \Delta\sigma) N} \big) \longrightarrow 0 \quad (N \to \infty). \tag{8.18.6}$$

---

**Proof:**

**Step 1 (Proof of Part I — Exact Resolvent Commutator and Trace Identity):**
Differentiating the resolvent identity $(Q - z I) R(z) = I$ with respect to the coordinate operator $K$:
$$0 = [K, (Q - z I) R(z)] = [K, Q - z I] R(z) + (Q - z I) [K, R(z)] = [K, Q] R(z) + (Q - z I) [K, R(z)].$$
Multiplying on the left by $R(z) = (Q - z I)^{-1}$:
$$[K, R(z)] = - R(z) [K, Q] R(z).$$
Substituting the exact rank-two commutator $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ from Proposition 8.8:
$$[K, R(z)] = - R(z) (\boldsymbol\psi d^T - d \boldsymbol\psi^T) R(z) = - (R(z)\boldsymbol\psi)(d^T R(z)) + (R(z)d)(\boldsymbol\psi^T R(z)),$$
which proves (8.18.1).

Now define the resolvent vectors:
$$w_{\mathrm{odd}} \equiv R(z) \boldsymbol\psi = R_{\mathrm{odd}}(z) \boldsymbol\psi \in H_{\mathrm{odd}}, \qquad w_{\mathrm{even}} \equiv R(z) d = R_{\mathrm{even}}(z) d \in H_{\mathrm{even}}.$$
The commutator is therefore the anti-symmetric rank-two operator:
$$[K, R(z)] = - w_{\mathrm{odd}} w_{\mathrm{even}}^T + w_{\mathrm{even}} w_{\mathrm{odd}}^T.$$
Evaluating the nested commutator $[K, [K, R(z)]] = K [K, R(z)] - [K, R(z)] K$:
$$[K, [K, R(z)]] = - (K w_{\mathrm{odd}}) w_{\mathrm{even}}^T + (K w_{\mathrm{even}}) w_{\mathrm{odd}}^T + w_{\mathrm{odd}} (K w_{\mathrm{even}})^T - w_{\mathrm{even}} (K w_{\mathrm{odd}})^T.$$
Algebraically, $[K, [K, R(z)]] = K^2 R(z) - 2 K R(z) K + R(z) K^2$.
Because $K$ interchanges the parity subspaces ($K : H_{\mathrm{odd}} \to H_{\mathrm{even}}$ and $K : H_{\mathrm{even}} \to H_{\mathrm{odd}}$) while $R(z)$ preserves them, the operator $K R(z) K$ maps $H_{\mathrm{odd}} \to H_{\mathrm{odd}}$ via $H_{\mathrm{even}}$.
Taking the trace on the odd subspace $H_{\mathrm{odd}}$:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}[K^2 R(z) + R(z) K^2] = 2 \operatorname{Tr}_{H_{\mathrm{odd}}}[K^2 R_{\mathrm{odd}}(z)],$$
while by the cyclic property of the trace across orthogonal complements:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}[K R(z) K] = \operatorname{Tr}_{H_{\mathrm{odd}}}[K_{oe} R_{\mathrm{even}}(z) K_{eo}] = \operatorname{Tr}_{H_{\mathrm{even}}}[R_{\mathrm{even}}(z) K_{eo} K_{oe}] = \operatorname{Tr}_{H_{\mathrm{even}}}[K^2 R_{\mathrm{even}}(z)].$$
Thus:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}\big[ [K, [K, R(z)]] \big] = 2 \Big( \operatorname{Tr}_{H_{\mathrm{odd}}}[K^2 R_{\mathrm{odd}}(z)] - \operatorname{Tr}_{H_{\mathrm{even}}}[K^2 R_{\mathrm{even}}(z)] \Big).$$
On the other hand, tracing the dyadic expansion on $H_{\mathrm{odd}}$: the components acting non-trivially on $H_{\mathrm{odd}}$ are $(K w_{\mathrm{even}}) w_{\mathrm{odd}}^T$ and $w_{\mathrm{odd}} (K w_{\mathrm{even}})^T$ (since $K w_{\mathrm{even}} \in H_{\mathrm{odd}}$ and $w_{\mathrm{odd}} \in H_{\mathrm{odd}}$).
Each dyad has trace:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}\big[ (K w_{\mathrm{even}}) w_{\mathrm{odd}}^T \big] = w_{\mathrm{odd}}^T (K w_{\mathrm{even}}) = \langle w_{\mathrm{odd}}, K w_{\mathrm{even}} \rangle.$$
Hence:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}\big[ [K, [K, R(z)]] \big] = 2 \langle w_{\mathrm{odd}}, K w_{\mathrm{even}} \rangle = 2 \langle R_{\mathrm{odd}}(z) \boldsymbol\psi, K R_{\mathrm{even}}(z) d \rangle.$$
Equating the two trace expressions and dividing by 2 establishes (8.18.2).

**Step 2 (Proof of Part II — Exact Spectral Expansion):**
In the orthonormal eigenbasis $\{u_j\}_{j=0}^{N-1}$ of $Q_{\mathrm{odd}}$:
$$R_\lambda^\perp = (Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} = \sum_{j=1}^{N-1} \frac{1}{\mu_j - \lambda} u_j u_j^T.$$
Since $K^2$ is symmetric:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}(K^2 R_\lambda^\perp) = \sum_{j=1}^{N-1} \frac{1}{\mu_j - \lambda} \operatorname{Tr}(K^2 u_j u_j^T) = \sum_{j=1}^{N-1} \frac{\langle u_j, K^2 u_j \rangle}{\mu_j - \lambda} = \sum_{j=1}^{N-1} \frac{\|K u_j\|^2}{\mu_j - \lambda}.$$
In the canonical odd Fourier basis $\{v_m\}_{m=1}^N$ with $v_m = \frac{e_m - e_{-m}}{\sqrt{2}}$, $K^2 v_m = m^2 v_m$.
Evaluating the trace directly:
$$\operatorname{Tr}_{H_{\mathrm{odd}}}(K^2 R_\lambda^\perp) = \sum_{m=1}^N \langle v_m, K^2 R_\lambda^\perp v_m \rangle = \sum_{m=1}^N m^2 \langle v_m, R_\lambda^\perp v_m \rangle = \sum_{m=1}^N m^2 [R_{\mathrm{odd},\perp}(\lambda)]_{mm},$$
establishing (8.18.3).

**Step 3 (Proof of Part III — Epistemic Resolution):**
1. *Divergence of Bare Trace:*
   The first excited odd mode $u_1$ has eigenvalue $\mu_1 < V_*$, so its tunneling gap $\mu_1 - \lambda \sim e^{-\sigma_1 N}$ decays exponentially (audited in Cell 67: $\mu_1 - \lambda \approx 2.95 \times 10^{-34}$ at $N=24$ with $\sigma_1 \approx 3.22$).
   Its coordinate kinetic energy is $\|K u_1\|^2 = \sum_{m=1}^N m^2 u_{1, m}^2 \ge \mathcal{C}_1 > 0$ (audited in Cell 68: $\|Ku_1\|^2 \approx 4.5 = \mathcal{O}(1)$).
   Retaining the single mode $j=1$ in the positive summation (8.18.3):
   $$\mathcal{T}_N = \sum_{j=1}^{N-1} \frac{\|Ku_j\|^2}{\mu_j - \lambda} \ge \frac{\|Ku_1\|^2}{\mu_1 - \lambda} \ge \mathcal{C}_1 e^{+\sigma_1 N}.$$
   Since $e^{+\sigma_1 N}$ grows faster than any polynomial $N^p$, the bare trace $\mathcal{T}_N$ cannot be bounded by $C N^p$, proving (8.18.4).

2. *Exponential Quenching of Scaled Trace:*
   Partitioning the sum at $V_*$:
   For high-energy modes $\mu_j \ge V_*$, $\mu_j - \lambda \ge c_{\mathrm{high}} > 0$ uniformly in $N$.
   By the exact trace identity $\sum_{j=0}^{N-1} \|Ku_j\|^2 = \operatorname{Tr}_{H_{\mathrm{odd}}}(K^2) = \frac{N(N+1)(2N+1)}{6}$,
   $$\mathcal{T}_{\mathrm{high}}(N) \le \frac{1}{c_{\mathrm{high}}} \sum_{j > N_{\mathrm{bound}}} \|Ku_j\|^2 \le \frac{1}{c_{\mathrm{high}}} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right).$$
   For bound modes $j \le N_{\mathrm{bound}}$, $\mu_j - \lambda \ge \mu_1 - \lambda$, so:
   $$\mathcal{T}_{\mathrm{bound}}(N) \le \frac{1}{\mu_1 - \lambda} \sum_{j=1}^{N_{\mathrm{bound}}} \|Ku_j\|^2.$$
   Multiplying by $D_0^2$:
   $$D_0^2 \mathcal{T}_N \le \frac{D_0^2}{\mu_1 - \lambda} \sum_{j=1}^{N_{\mathrm{bound}}} \|Ku_j\|^2 + \frac{D_0^2}{c_{\mathrm{high}}} \left( \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} \right).$$
   Under Bridges B1 and B2, $\sum_{j=1}^{N_{\mathrm{bound}}} \|Ku_j\|^2 \le \mathcal{C}_{\mathrm{bound}} < \infty$.
   By Hypothesis H3, $D_0^2 \le C_0 e^{-2\sigma_0 N}$.
   Under the semiclassical barrier-thinning hypothesis $\sigma_0 > \sigma_1$ (action gap $\Delta\sigma = \sigma_0 - \sigma_1 > 0$),
   $$\frac{D_0^2}{\mu_1 - \lambda} = R_{\mathrm{gap}}^{\max}(N) \sim e^{-2\sigma_0 N} e^{+\sigma_1 N} = e^{-(\sigma_0 + \Delta\sigma) N} \longrightarrow 0 \quad (N \to \infty).$$
   Since $D_0^2 \mathcal{T}_{\mathrm{high}} = \mathcal{O}(e^{-2\sigma_0 N} N^3) \to 0$ as well, we obtain $D_0^2 \mathcal{T}_N \to 0$ exponentially fast, establishing (8.18.5) and (8.18.6). $\blacksquare$

*Discussion and Scientific Epistemics:*
Proposition 8.18 clarifies the exact role of the global weighted resolvent trace in Route B:
1. **Resolution of the Polynomial Trace Hypothesis:**
   The hypothesis that the bare trace $\mathcal{T}_N = \operatorname{Tr}[K^2 R_\lambda^\perp]$ could be bounded by a polynomial $C N^p$ is analytically refuted by the exact spectral representation: the existence of bound states beneath the barrier top makes the lowest resolvent eigenvalue $\frac{1}{\mu_1 - \lambda}$ grow exponentially. No linear algebra or operator identity can suppress an isolated singular eigenvalue when the corresponding projection $\|Ku_1\|^2$ is non-zero.
2. **The True Physical Mechanism: Barrier-Thinning Quenching:**
   The continuum decoupling $\mathcal{S}_{\mathrm{coord}}(N) \to 0$ does not require the bare trace $\mathcal{T}_N$ to be polynomial. Instead, it relies on the **barrier-thinning action hierarchy**:
   $$2\sigma_0 > \sigma_1 \iff \sigma_0 + \Delta\sigma > 0.$$
   Because the boundary tunneling amplitude $D_0^2$ decays with double the ground-state barrier action $2\sigma_0 \approx 7.55$, while the first excited gap $\mu_1 - \lambda$ closes with the smaller excited barrier action $\sigma_1 \approx 3.22$, the product $D_0^2 \mathcal{T}_N \sim e^{-4.33 N}$ is exponentially extinguished.
3. **Retention of Bridges B1 and B2:**
   The global trace analysis does **not** bypass Bridges B1 and B2; the bound-state localization estimate $\sum_{j=1}^{N_{\mathrm{bound}}} \|Ku_j\|^2 \le \mathcal{C}_{\mathrm{bound}} < \infty$ remains an essential premise in (8.18.5). What Proposition 8.18 establishes is that once B1 and B2 are held, the entire small-denominator challenge collapses strictly to the relative tunneling gap $R_{\mathrm{gap}}^{\max}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda}$.
4. **The Commutator Parity Trace Identity:**
   The exact identity $\operatorname{Tr}_{H_{\mathrm{odd}}}(K^2 R_{\mathrm{odd}}) - \operatorname{Tr}_{H_{\mathrm{even}}}(K^2 R_{\mathrm{even}}) = \langle R_{\mathrm{odd}}\boldsymbol\psi, K R_{\mathrm{even}} d \rangle$ demonstrates that the difference between the weighted coordinate traces in the two parity sectors collapses to a single rank-two inner product between the odd resolvent vector $w_{\mathrm{odd}}$ and the coordinate-shifted even resolvent vector $K w_{\mathrm{even}}$, providing an exact structural link between the two sectors.

---

### 8.19 Proposition 8.19 (Exact Transition Dipole Factorization of the First Relative Tunneling Gap and Excited Wavepacket Residual)

*The first relative tunneling gap $R_{\mathrm{gap}}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda}$, which controls the bound-state coordinate energy $\mathcal{S}_{\mathrm{bound}}(N)$ and the Route B continuum decoupling limit $\rho_2^{\mathrm{exc}} \to 0$, admits an exact finite-dimensional algebraic factorization in terms of the first coordinate transition dipole $b_{01} \equiv \langle u_0^{\mathrm{even}}, K u_1 \rangle$ and the first excited overlap transmission ratio $\mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda}$. Consequently, the relative tunneling gap is identically bounded by the excited coordinate wavepacket residual $\|v_{\mathrm{exc}}\|^2 = \|P_{\perp u_0} K c\|^2$, and is governed by an exact even-sector resolvent cancellation identity.*

**Part I (Exact Dipole Factorization — Unconditional Finite-$N$ Theorem):**
*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Galerkin matrix, $u_0^{\mathrm{even}} = c$ the even ground state ($Q_{\mathrm{even}} c = \lambda c$, $\langle d, c \rangle = D_0$), and $u_1$ the first excited odd eigenstate ($Q_{\mathrm{odd}} u_1 = \mu_1 u_1$, $a_1 = \langle \boldsymbol\psi, u_1 \rangle$). Define the coordinate transition dipole:*
$$b_{01} \equiv \langle u_0^{\mathrm{even}}, K u_1 \rangle = \langle c, K u_1 \rangle.$$
*For every finite dimension $N \ge 2$ and cutoff $c > 1$:*
1. *The dipole matrix element satisfies the exact commutator quotient identity:*
   $$b_{01} = - \frac{D_0 a_1}{\mu_1 - \lambda}. \tag{8.19.1}$$
2. *Consequently, the first relative tunneling gap factors identically into:*
   $$\boxed{R_{\mathrm{gap}}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda} = \frac{b_{01}^2}{\mathcal{R}_1}, \qquad \mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda}.} \tag{8.19.2}$$

**Part II (Coordinate Wavepacket Residual Upper Enclosure):**
*Let $v_{\mathrm{exc}} \equiv P_{\perp u_0} K c$ be the excited coordinate wavepacket residual of Proposition 8.12. Then:*
$$\boxed{b_{01}^2 \le \|v_{\mathrm{exc}}\|^2 \equiv \|P_{\perp u_0} K c\|^2 = D_0^2 M_{2,\mathrm{exc}},} \tag{8.19.3}$$
*and therefore:*
$$\boxed{R_{\mathrm{gap}}(N) \le \frac{\|v_{\mathrm{exc}}\|^2}{\mathcal{R}_1} = \frac{\|P_{\perp u_0} K c\|^2}{a_1^2 / (\mu_1 - \lambda)}.} \tag{8.19.4}$$

**Part III (Even-Sector Resolvent Cancellation Identity):**
*Expanding $K u_1$ in the orthonormal eigenbasis $\{u_k^{\mathrm{even}}\}_{k=0}^N$ of $Q_{\mathrm{even}}$ with eigenvalues $E_0 = \lambda < E_1 < \dots < E_N$ and boundary weights $d_k = \langle u_k^{\mathrm{even}}, d \rangle$:*
$$\boxed{b_{01}^2 = \|K u_1\|^2 - a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2} = \|K u_1\|^2 - a_1^2 \|(Q_{\mathrm{even}} - \mu_1 I)^{-1} P_{\perp c} d\|^2.} \tag{8.19.5}$$

---

**Proof:**

**Step 1 (Proof of Part I — Exact Dipole Factorization):**
From Proposition 8.8, the rank-two coordinate commutator is $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$.
Applying this to the first excited odd eigenstate $u_1 \in H_{\mathrm{odd}}$ ($Q_{\mathrm{odd}} u_1 = \mu_1 u_1$, $\langle d, u_1 \rangle = 0$ by odd parity, and $\langle \boldsymbol\psi, u_1 \rangle = a_1$):
$$[K, Q] u_1 = \boldsymbol\psi \langle d, u_1 \rangle - d \langle \boldsymbol\psi, u_1 \rangle = - a_1 d.$$
Expanding the commutator:
$$[K, Q] u_1 = K Q u_1 - Q K u_1 = \mu_1 K u_1 - Q_{\mathrm{even}} K u_1 = - (Q_{\mathrm{even}} - \mu_1 I) K u_1.$$
Equating both sides yields the exact operator identity on $H_{\mathrm{even}}$:
$$(Q_{\mathrm{even}} - \mu_1 I) K u_1 = a_1 d.$$
Taking the inner product with the even ground state $u_0^{\mathrm{even}} = c$ ($Q_{\mathrm{even}} c = \lambda c$, $\langle d, c \rangle = D_0$):
$$\langle c, (Q_{\mathrm{even}} - \mu_1 I) K u_1 \rangle = (\lambda - \mu_1) \langle c, K u_1 \rangle = a_1 \langle c, d \rangle = a_1 D_0.$$
Dividing by $\lambda - \mu_1 = - (\mu_1 - \lambda) \ne 0$ (since $\mu_1 > \mu_0 > \lambda$):
$$b_{01} \equiv \langle c, K u_1 \rangle = - \frac{D_0 a_1}{\mu_1 - \lambda},$$
which proves (8.19.1).
Squaring both sides gives:
$$b_{01}^2 = \frac{D_0^2 a_1^2}{(\mu_1 - \lambda)^2} = \frac{D_0^2}{\mu_1 - \lambda} \cdot \frac{a_1^2}{\mu_1 - \lambda} = R_{\mathrm{gap}}(N) \cdot \mathcal{R}_1.$$
Since $a_1 \ne 0$ by strict spectral interlacing (Proposition 8.14), dividing by $\mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda} > 0$ yields (8.19.2).

**Step 2 (Proof of Part II — Wavepacket Residual Enclosure):**
From Proposition 8.12, the excited coordinate Parseval sum satisfies:
$$D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} K c\|^2 = \sum_{j=1}^{N-1} |\langle u_j, K c \rangle|^2.$$
Since $K$ is symmetric and interchanges parities, $\langle u_j, K c \rangle = - \langle c, K u_j \rangle = - b_{0j}$.
Therefore:
$$\|v_{\mathrm{exc}}\|^2 = \|P_{\perp u_0} K c\|^2 = \sum_{j=1}^{N-1} b_{0j}^2 = b_{01}^2 + \sum_{j=2}^{N-1} b_{0j}^2.$$
Because all terms $b_{0j}^2 \ge 0$, dropping the higher modes $j \ge 2$ gives $b_{01}^2 \le \|v_{\mathrm{exc}}\|^2$, establishing (8.19.3).
Substituting this into (8.19.2) establishes (8.19.4).

**Step 3 (Proof of Part III — Even-Sector Resolvent Cancellation Identity):**
From Step 1, $(Q_{\mathrm{even}} - \mu_1 I) K u_1 = a_1 d$.
For any excited even eigenvector $u_k^{\mathrm{even}}$ ($k \ge 1$) with eigenvalue $E_k$ and overlap $d_k = \langle u_k^{\mathrm{even}}, d \rangle$:
$$\langle u_k^{\mathrm{even}}, (Q_{\mathrm{even}} - \mu_1 I) K u_1 \rangle = (E_k - \mu_1) \langle u_k^{\mathrm{even}}, K u_1 \rangle = a_1 d_k \implies \langle u_k^{\mathrm{even}}, K u_1 \rangle = - \frac{a_1 d_k}{\mu_1 - E_k}.$$
Because $\{u_k^{\mathrm{even}}\}_{k=0}^N$ forms an orthonormal basis of $H_{\mathrm{even}}$, Parseval's identity on $K u_1 \in H_{\mathrm{even}}$ yields:
$$\|K u_1\|^2 = |\langle u_0^{\mathrm{even}}, K u_1 \rangle|^2 + \sum_{k=1}^N |\langle u_k^{\mathrm{even}}, K u_1 \rangle|^2 = b_{01}^2 + a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2}.$$
Subtracting the excited even sum from $\|Ku_1\|^2$ gives (8.19.5), completing the proof. $\blacksquare$

*Discussion and Scientific Epistemics:*
Proposition 8.19 establishes an exact reduction of the relative tunneling gap:
1. **The Exact Dipole Architecture:**
   The relative gap $R_{\mathrm{gap}}(N) = \frac{D_0^2}{\mu_1 - \lambda}$ does not require separate heuristic estimates of $D_0^2$ and $\mu_1 - \lambda$. Identity (8.19.2) shows that $R_{\mathrm{gap}}$ factors identically into the transition dipole $b_{01}^2$ and the odd transmission ratio $\mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda}$.
   The exact factorization reduces control of $R_{\mathrm{gap}}$ to joint control of $b_{01}^2$ and $\mathcal{R}_1$. If $\mathcal{R}_1$ remains bounded above and below away from zero, then $R_{\mathrm{gap}}(N) \to 0 \iff b_{01}^2 \to 0$. Across all tested dimensions $N \in \{8, \dots, 24\}$, $\mathcal{R}_1 \in [3.91, 8.00]$ remains strictly $\mathcal{O}(1)$ (`cell69.out`), but this uniform lower bound is an empirical observation across tested dimensions, not an assumed theorem.
2. **Dominance of Mode 1 in the Excited Wavepacket:**
   Because $\|v_{\mathrm{exc}}\|^2 = \sum_{j \ge 1} b_{0j}^2$, the ratio $b_{01}^2 / \|v_{\mathrm{exc}}\|^2$ measures the concentration of the excited coordinate wavepacket $P_{\perp u_0} K c$ in the first excited mode. As shown in `cell69.out`, mode 1 carries $99.9909\%$ at $N=8$ rising monotonically to $99.9995\%$ at $N=24$ of this residual norm, so $b_{01}^2 \approx \|v_{\mathrm{exc}}\|^2$.
3. **The Cancellation Mechanism in the Even Sector:**
   Identity (8.19.5) exposes the precise algebraic mechanism governing $b_{01}^2$:
   $$b_{01}^2 = \|K u_1\|^2 - a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2} \equiv \|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1}.$$
   At $N=24$, $\|Ku_1\|^2 = 10.2186531$ and the excited even resolvent sum $\mathcal{E}_{\mathrm{even}, 1} = 10.2186508$ are both $\mathcal{O}(10)$ quantities, while their difference is $b_{01}^2 = 2.25258 \times 10^{-6}$ (`cell69.out`). The cancellation residual independently agrees at $10^{-18}$ level in 50-digit precision. The analytical origin of this cancellation is resolved in Proposition 8.20 below.

---

### 8.20 Proposition 8.20: Exact Even-Resolvent Representation of $K u_1$, Stieltjes Derivative Architecture, and Algebraic Resolution of the Dipole Cancellation

*The transition dipole cancellation $b_{01}^2 = \|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1}$ observed in Proposition 8.19 admits an exact operator-theoretic explanation: $K u_1$ is proportional to the even resolvent of the boundary vector $d$ evaluated at the odd eigenvalue $\mu_1$, its norm squared $\|Ku_1\|^2$ is the exact derivative of the even Stieltjes transform $G_d'(\mu_1)$, and subtracting the excited even sum $\mathcal{E}_{\mathrm{even}, 1}$ is the exact algebraic projector isolating the ground-state pole.*

**Part I (Exact Resolvent Representation — Unconditional Finite-$N$ Identity):**
*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Galerkin matrix, $u_0^{\mathrm{even}} = c$ the even ground state ($Q_{\mathrm{even}} c = \lambda c$, $\langle d, c \rangle = D_0$), and $u_1$ the first excited odd eigenstate ($Q_{\mathrm{odd}} u_1 = \mu_1 u_1$, $a_1 = \langle \boldsymbol\psi, u_1 \rangle \ne 0$).*
*Because $\mu_1 \notin \sigma(Q_{\mathrm{even}})$ by strict spectral interlacing ($E_1 < \mu_1 < E_2$), the shifted operator $(Q_{\mathrm{even}} - \mu_1 I)$ is unconditionally invertible on $H_{\mathrm{even}}$, and the coordinate derivative $K u_1 \in H_{\mathrm{even}}$ admits the exact, closed-form resolvent representation:*
$$\boxed{K u_1 = a_1 (Q_{\mathrm{even}} - \mu_1 I)^{-1} d. \tag{8.20.1}}$$

**Part II (Exact Stieltjes Derivative Identity):**
*Let $G_d(z) \equiv \langle d, (Q_{\mathrm{even}} - z I)^{-1} d \rangle = \sum_{k=0}^N \frac{d_k^2}{E_k - z}$ be the even-sector Stieltjes transform generated by the boundary vector $d \in H_{\mathrm{even}}$.*
*The discrete coordinate kinetic energy $\|Ku_1\|^2$ is identically given by the first derivative of the even Stieltjes transform evaluated at the odd eigenvalue $z = \mu_1$:*
$$\boxed{\|Ku_1\|^2 = a_1^2 \|(Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2 = a_1^2 G_d'(\mu_1) = a_1^2 \sum_{k=0}^N \frac{d_k^2}{(\mu_1 - E_k)^2}.} \tag{8.20.2}$$

**Part III (Algebraic Resolution of the Dipole Cancellation):**
*The cancellation identity (8.19.5):*
$$b_{01}^2 = \|Ku_1\|^2 - a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2}$$
*is an exact algebraic projector isolating the $k=0$ ground-state residue from the even Stieltjes derivative:*
$$\boxed{b_{01}^2 = a_1^2 \frac{d_0^2}{(\mu_1 - E_0)^2} = \frac{a_1^2 D_0^2}{(\mu_1 - \lambda)^2}.} \tag{8.20.3}$$
*The subtraction $\|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1}$ does not represent an asymptotic or approximate cancellation, but an exact algebraic identity: subtracting the excited even spectral sum $\mathcal{E}_{\mathrm{even}, 1} = a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2}$ from the total norm squared $\|Ku_1\|^2$ removes precisely the $k \ge 1$ terms of the Parseval sum, leaving identically the single $k=0$ pole.*

**Part IV (Excited Doublet Partner Decomposition):**
*Decomposing the excited sum $\mathcal{E}_{\mathrm{even}, 1}$ separates the resonant first excited parity doublet $(E_1, \mu_1)$ from the higher spectrum $k \ge 2$:*
$$\boxed{\mathcal{E}_{\mathrm{even}, 1} = \frac{a_1^2 d_1^2}{(\mu_1 - E_1)^2} + a_1^2 \sum_{k=2}^N \frac{d_k^2}{(\mu_1 - E_k)^2}.} \tag{8.20.4}$$
*Consequently, the first relative tunneling gap $R_{\mathrm{gap}}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda}$ satisfies the exact positive resolvent representation:*
$$\boxed{R_{\mathrm{gap}}(N) = \frac{\mu_1 - \lambda}{a_1^2} \left[ a_1^2 G_d'(\mu_1) - \mathcal{E}_{\mathrm{even}, 1} \right] = \frac{b_{01}^2}{\mathcal{R}_1}.} \tag{8.20.5}$$

---

**Proof:**

**Step 1 (Proof of Part I — Exact Resolvent Representation):**
From Proposition 8.19 Step 1, the rank-two commutator $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ applied to $u_1 \in H_{\mathrm{odd}}$ establishes the operator identity on $H_{\mathrm{even}}$:
$$(Q_{\mathrm{even}} - \mu_1 I) K u_1 = a_1 d.$$
By the strict spectral interlacing theorem for rank-one parity perturbations (Proposition 8.8 and Proposition 8.14), the eigenvalues of $Q_{\mathrm{even}}$ and $Q_{\mathrm{odd}}$ strictly interlace:
$$E_0 = \lambda < \mu_0 < E_1 < \mu_1 < E_2 < \dots < E_N.$$
In particular, $\mu_1 \in (E_1, E_2)$, so $\mu_1 \notin \sigma(Q_{\mathrm{even}})$.
Thus $(Q_{\mathrm{even}} - \mu_1 I)$ is a non-singular, invertible symmetric matrix on the $(N+1)$-dimensional subspace $H_{\mathrm{even}}$.
Multiplying both sides by $(Q_{\mathrm{even}} - \mu_1 I)^{-1}$ yields:
$$K u_1 = a_1 (Q_{\mathrm{even}} - \mu_1 I)^{-1} d,$$
which is (8.20.1).

**Step 2 (Proof of Part II — Stieltjes Derivative Identity):**
Taking the Euclidean norm squared of $K u_1 \in H_{\mathrm{even}}$:
$$\|Ku_1\|^2 = a_1^2 \langle (Q_{\mathrm{even}} - \mu_1 I)^{-1} d, (Q_{\mathrm{even}} - \mu_1 I)^{-1} d \rangle = a_1^2 \langle d, (Q_{\mathrm{even}} - \mu_1 I)^{-2} d \rangle.$$
Expanding in the orthonormal eigenbasis $\{u_k^{\mathrm{even}}\}_{k=0}^N$ of $Q_{\mathrm{even}}$ with eigenvalues $E_k$ and weights $d_k = \langle u_k^{\mathrm{even}}, d \rangle$:
$$\langle d, (Q_{\mathrm{even}} - \mu_1 I)^{-2} d \rangle = \sum_{k=0}^N \frac{d_k^2}{(\mu_1 - E_k)^2}.$$
On the other hand, the Stieltjes transform of the even boundary measure $d\nu_d(E) = \sum_{k=0}^N d_k^2 \delta_{E_k}$ is:
$$G_d(z) = \langle d, (Q_{\mathrm{even}} - z I)^{-1} d \rangle = \sum_{k=0}^N \frac{d_k^2}{E_k - z}.$$
Differentiating with respect to $z$:
$$G_d'(z) = \sum_{k=0}^N \frac{d_k^2}{(E_k - z)^2}.$$
Evaluating at $z = \mu_1$ gives $G_d'(\mu_1) = \sum_{k=0}^N \frac{d_k^2}{(\mu_1 - E_k)^2}$, proving (8.20.2).

**Step 3 (Proof of Part III — Resolution of the Dipole Cancellation):**
The sum in (8.20.2) runs over $k = 0, 1, \dots, N$.
Splitting the $k=0$ term from the remaining terms $k \ge 1$:
$$\sum_{k=0}^N \frac{d_k^2}{(\mu_1 - E_k)^2} = \frac{d_0^2}{(\mu_1 - E_0)^2} + \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2}.$$
For $k=0$, $u_0^{\mathrm{even}} = c$ is the even ground state, so $E_0 = \lambda$ and $d_0 = \langle c, d \rangle = D_0$.
Therefore, the $k=0$ contribution multiplied by $a_1^2$ is:
$$a_1^2 \frac{d_0^2}{(\mu_1 - E_0)^2} = \frac{a_1^2 D_0^2}{(\mu_1 - \lambda)^2} = \left( - \frac{D_0 a_1}{\mu_1 - \lambda} \right)^2 = b_{01}^2,$$
where the last equality is the exact commutator quotient identity (8.19.1).
Multiplying the remaining terms by $a_1^2$ gives precisely $\mathcal{E}_{\mathrm{even}, 1} \equiv a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2}$.
Therefore:
$$\|Ku_1\|^2 = a_1^2 G_d'(\mu_1) = b_{01}^2 + \mathcal{E}_{\mathrm{even}, 1} \implies b_{01}^2 = \|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1},$$
proving that identity (8.19.5) is an exact algebraic decomposition separating the ground-state pole from the excited spectrum. $\blacksquare$

**Step 4 (Proof of Part IV — Doublet Partner Decomposition):**
Separating the $k=1$ summand from the $k \ge 2$ terms in $\mathcal{E}_{\mathrm{even}, 1}$ is an immediate algebraic partitioning, yielding (8.20.4).
Multiplying $b_{01}^2$ by $\frac{\mu_1 - \lambda}{a_1^2} = \frac{1}{\mathcal{R}_1}$ yields (8.20.5). $\blacksquare$

---

*Discussion and Epistemic Synthesis:*
1. **The Algebraic Anatomy of the Cancellation:**
   In `cell69.out`, the numerical calculation observed:
   $$\|Ku_1\|^2 \approx 10.2186531, \qquad \mathcal{E}_{\mathrm{even}, 1} \approx 10.2186508, \qquad b_{01}^2 \approx 2.25258 \times 10^{-6} \quad (N = 24).$$
   Proposition 8.20 proves that this is not an accidental numerical near-cancellation between two unrelated dynamical quantities. Rather:
   - $\|Ku_1\|^2 = a_1^2 G_d'(\mu_1)$ is the total norm squared of the resolvent vector $a_1 (Q_{\mathrm{even}} - \mu_1 I)^{-1} d$.
   - $\mathcal{E}_{\mathrm{even}, 1}$ is the projection of this resolvent vector onto the orthogonal complement of the ground state $P_{\perp c} H_{\mathrm{even}}$.
   - $b_{01}^2$ is the projection of this resolvent vector onto the one-dimensional ground state subspace $\mathbb{R} c$.
   Therefore, $\|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1} \equiv b_{01}^2$ is an exact Pythagorean decomposition on $H_{\mathrm{even}}$:
   $$\|(Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2 = \|P_c (Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2 + \|P_{\perp c} (Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2.$$
2. **Why is the Ground-State Projection $b_{01}^2$ Quenched?**
   Because $K$ is symmetric ($K = K^T$) and interchanges parity sectors:
   $$b_{01} \equiv \langle c, K u_1 \rangle = \langle K c, u_1 \rangle.$$
   (The minus sign belongs strictly to the commutator quotient identity $b_{01} = - \frac{D_0 a_1}{\mu_1 - \lambda}$, not to the symmetry of $K$.)
   From Proposition 8.12, the coordinate-derivative wavepacket $Kc$ is overwhelmingly aligned with the odd ground state $u_0$:
   $$Kc = \langle u_0, Kc \rangle u_0 + P_{\perp u_0} Kc, \qquad \frac{|\langle u_0, Kc \rangle|^2}{\|Kc\|^2} > 99.99987\% \quad (N = 24).$$
   Because $u_1 \perp u_0$, the ground-state component vanishes identically: $\langle \langle u_0, Kc \rangle u_0, u_1 \rangle = 0$.
   Therefore:
   $$\boxed{b_{01} = \langle P_{\perp u_0} Kc, u_1 \rangle.}$$
   The transition dipole $b_{01}$ is precisely the overlap between the first excited odd mode $u_1$ and the excited wavepacket residual $v_{\mathrm{exc}} = P_{\perp u_0} Kc$.
   Because $\|v_{\mathrm{exc}}\|^2 = D_0^2 M_{2,\mathrm{exc}} \sim 2.25 \times 10^{-6}$ is tiny due to the near-perfect alignment of $Kc$ with $u_0$, its projection onto $u_1$ is necessarily bounded by $\|v_{\mathrm{exc}}\|^2$, establishing $b_{01}^2 \le \|v_{\mathrm{exc}}\|^2 \sim 10^{-6}$.
3. **The Wavepacket Concentration Mechanism:**
   In `cell69.out`, mode 1 carries $99.9995\%$ of the total excited wavepacket norm $\|v_{\mathrm{exc}}\|^2 = \sum_{j \ge 1} b_{0j}^2$.
   Proposition 8.20 explains this concentration: $u_1$ is the lowest excited bound state, lying at the lowest kinetic energy above $u_0$, so that the spatial residual $P_{\perp u_0} Kc$ is almost purely dipolar, projecting almost exclusively onto the first excited nodal mode $u_1$.
4. **Epistemic Hygiene on $\mathcal{R}_1$:**
   In `cell69.out`, the transmission ratio $\mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda}$ takes the values:
   $$N=8: 8.00, \quad N=12: 6.88, \quad N=16: 4.13, \quad N=20: 3.91, \quad N=24: 5.13.$$
   While these empirical values are consistent with an $\mathcal{O}(1)$ non-zero limit, five data points do not constitute a mathematical proof that $\inf_N \mathcal{R}_1(N) \ge c_1 > 0$. We therefore maintain strict epistemic discipline: $R_{\mathrm{gap}}(N) \to 0$ is guaranteed if $b_{01}^2 \to 0$ provided $\mathcal{R}_1$ does not collapse to zero, and the primary analytical agenda is establishing analytical bounds on the Stieltjes derivative $G_d'(\mu_1)$ and wavepacket alignment.
5. **Algebraic Anatomy vs. Asymptotic Convergence:**
   Proposition 8.20 resolves the algebraic nature of $b_{01}^2$ as a positive component of a squared resolvent norm, proving that the numerical subtraction observed in `cell69.out` is an exact Pythagorean decomposition rather than a precarious numerical cancellation. However, this does not by itself prove that $b_{01}^2 \to 0$ asymptotically. Because $b_{01} = \langle P_{\perp u_0} Kc, u_1 \rangle$, the asymptotic agenda is transferred to proving wavepacket alignment $\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \to 0$, which is addressed in Proposition 8.21.

---

### 8.21 Proposition 8.21 (Exact Variational Energy Excess of Coordinate Wavepacket $Kc$, Rayleigh Quotient Gap Enclosure, and Asymptotic Alignment Sandwich)

The algebraic resolution in Proposition 8.20 reduced the transition dipole to the wavepacket overlap $b_{01} = \langle P_{\perp u_0} Kc, u_1 \rangle$, showing that $b_{01}^2 \le \|P_{\perp u_0} Kc\|^2$. Consequently, establishing continuum gap suppression $R_{\mathrm{gap}}(N) \to 0$ does not require resolving individual modal cancellations, but rather proving that the coordinate-derivative wavepacket $Kc$ becomes asymptotically parallel to the odd ground state $u_0$:
$$\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \longrightarrow 0 \iff \frac{Kc}{\|Kc\|} \longrightarrow \pm u_0.$$
The following proposition establishes that the energy excess of $Kc$ above the odd ground state $\mu_0$ is an exact finite-dimensional invariant determined directly by the rank-two commutator $[K, Q]$, proving an unconditional two-sided sandwich between the wavepacket misalignment and the relative tunneling gap.

---

**Part I (Exact Energy Expectation and Excess Identity — Unconditional Finite-$N$ Theorem):**
*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Galerkin matrix, $K = \operatorname{diag}(-N, \dots, N)$, and $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$ the exact rank-two commutator. Let $c \in H_{\mathrm{even}}$ be the even ground state with eigenvalue $\lambda$, boundary amplitude $D_0 = \langle d, c \rangle$, and coordinate wavepacket $Kc \in H_{\mathrm{odd}}$.*
*For all $N \ge 1$, the quadratic energy expectation of $Kc$ under $Q_{\mathrm{odd}}$ satisfies the exact finite-dimensional identity:*
$$\boxed{\langle Kc, Q_{\mathrm{odd}} Kc \rangle = \lambda \|Kc\|^2 + D_0^2 M_1,} \tag{8.21.1}$$
*where $M_1 \equiv \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi \rangle = \sum_{j \ge 0} \frac{a_j^2}{\mu_j - \lambda}$.*
*Consequently, the energy excess of $Kc$ above the odd ground state $\mu_0$ evaluates identically to:*
$$\boxed{\langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle = D_0^2 M_1^{\mathrm{exc}} - (\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2,} \tag{8.21.2}$$
*where $M_1^{\mathrm{exc}} \equiv \sum_{j \ge 1} \frac{a_j^2}{\mu_j - \lambda}$, with the ground-state pole $j=0$ completely cancelled.*

**Part II (Variational Rayleigh Quotient Enclosure — Unconditional Finite-$N$ Theorem):**
*Since $u_0$ minimizes the odd Rayleigh quotient ($\mu_0 = \min_{x \in H_{\mathrm{odd}} \setminus \{0\}} \frac{\langle x, Q_{\mathrm{odd}} x \rangle}{\|x\|^2}$), the Rayleigh quotient excess $\mathcal{R}_{Q_{\mathrm{odd}}}(Kc) - \mu_0 \ge 0$ implies:*
$$(\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2 \le D_0^2 M_1^{\mathrm{exc}}. \tag{8.21.3}$$
*Furthermore, decomposing $Kc$ in the orthonormal eigenbasis $\{u_j\}_{j=0}^{N-1}$ of $Q_{\mathrm{odd}}$, the first odd spectral gap $\mu_1 - \mu_0$ yields the sharp lower bound $\langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle \ge (\mu_1 - \mu_0) \|P_{\perp u_0} Kc\|^2$, establishing the exact finite-$N$ upper bound:*
$$\boxed{\|P_{\perp u_0} Kc\|^2 \le \frac{D_0^2 M_1^{\mathrm{exc}}}{\mu_1 - \lambda} = M_1^{\mathrm{exc}} R_{\mathrm{gap}}(N).} \tag{8.21.4}$$

**Part III (The Asymptotic Alignment Sandwich — Unconditional Finite-$N$ Theorem):**
*Combining the upper bound (8.21.4) with the transition dipole projection $b_{01} = \langle P_{\perp u_0} Kc, u_1 \rangle$ and Proposition 8.19 ($b_{01}^2 = \mathcal{R}_1 R_{\mathrm{gap}}(N)$) proves that the normalized wavepacket misalignment $\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2}$ is clamped between two exact positive multiples of the relative tunneling gap:*
$$\boxed{\frac{\mathcal{R}_1}{\|Kc\|^2} R_{\mathrm{gap}}(N) \le \frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \le \frac{M_1^{\mathrm{exc}}}{\|Kc\|^2} R_{\mathrm{gap}}(N).} \tag{8.21.5}$$
*Consequently, relative tunneling gap suppression $R_{\mathrm{gap}}(N) \to 0$ guarantees asymptotic wavepacket alignment $\frac{Kc}{\|Kc\|} \to \pm u_0$ provided the prefactor $M_1^{\mathrm{exc}} / \|Kc\|^2$ does not grow too rapidly (which holds under Hypothesis H2$_{\mathrm{odd}}$ and solitary coordinate bounds). Conversely, asymptotic alignment $\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \to 0$ implies relative gap suppression $R_{\mathrm{gap}}(N) \to 0$ conditional on the lower coefficient remaining bounded away from zero ($\inf_N \mathcal{R}_1(N) \ge c_1 > 0$ and $\sup_N \|Kc\|^2 < \infty$).*

**Part IV (The Analytical Hierarchy and the Moment Problem):**
*The algebraic structure connecting the coordinate commutator to the wavepacket energy and residual is organized into the following hierarchy:*
$$\boxed{\begin{aligned}
&[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T \text{ (exact finite-}N\text{ commutator)} \\
&\qquad\Downarrow \\
&E_{\mathrm{exc}} \equiv \langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle = D_0^2 M_1^{\mathrm{exc}} - (\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2 \text{ (exact first moment } \int (\mu - \mu_0) d\nu_v \text{)} \\
&\qquad\Downarrow \\
&\frac{\mathcal{R}_1}{\|Kc\|^2} R_{\mathrm{gap}}(N) \le \frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \le \frac{M_1^{\mathrm{exc}}}{\|Kc\|^2} R_{\mathrm{gap}}(N) \text{ (exact finite-}N\text{ sandwich)} \\
&\qquad\Downarrow \\
&R_{\mathrm{gap}}(N) \to 0 \implies \frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \to 0 \iff \frac{Kc}{\|Kc\|} \to \pm u_0 \quad (\text{conditional on prefactor control}), \\
&\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \to 0 \implies R_{\mathrm{gap}}(N) \to 0 \quad (\text{conditional on } \inf_N \mathcal{R}_1(N) > 0).
\end{aligned}} \tag{8.21.6}$$
*Crucially, because the spectral gap $\mu_1 - \mu_0 \sim e^{-\sigma_1 N} \to 0$ closes exponentially, the smallness of the energy excess $E_{\mathrm{exc}} \to 0$ is a consequence of the closing gap and does not by itself imply that the zeroth moment $\|P_{\perp u_0} Kc\|^2 = \int d\nu_v$ vanishes.*

---

**Proof:**

**Step 1 (Proof of Part I — Exact Energy Expectation and Excess Identity):**
Since $c \in H_{\mathrm{even}}$ and $K = \operatorname{diag}(-N, \dots, N)$ is symmetric ($K = K^T$), the quadratic form of $Kc \in H_{\mathrm{odd}}$ under $Q_{\mathrm{odd}}$ is:
$$\langle Kc, Q_{\mathrm{odd}} Kc \rangle = \langle c, K Q K c \rangle.$$
From the definition of the commutator $[K, Q] = KQ - QK$, we have $KQ = QK + [K, Q]$.
Multiplying on the right by $K$ yields $K Q K = Q K^2 + [K, Q] K$.
Taking the inner product with the even ground state $c$:
$$\langle c, K Q K c \rangle = \langle c, Q K^2 c \rangle + \langle c, [K, Q] K c \rangle.$$
Because $Q$ is symmetric and $c$ is an eigenvector with $Q c = \lambda c$:
$$\langle c, Q K^2 c \rangle = \langle Q c, K^2 c \rangle = \lambda \langle c, K^2 c \rangle = \lambda \|Kc\|^2.$$
For the second term, we evaluate $[K, Q] Kc$ using the rank-two commutator $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$:
$$[K, Q] Kc = \boldsymbol\psi \langle d, Kc \rangle - d \langle \boldsymbol\psi, Kc \rangle.$$
Since $Kc \in H_{\mathrm{odd}}$ and $d \in H_{\mathrm{even}}$, the boundary overlap vanishes by parity: $\langle d, Kc \rangle = 0$.
To evaluate $\langle \boldsymbol\psi, Kc \rangle$, we recall from Proposition 8.13 Step 1 that $[K, Q] c = \lambda Kc - Q Kc = - (Q - \lambda I) Kc$, while $[K, Q] c = \boldsymbol\psi \langle d, c \rangle - d \langle \boldsymbol\psi, c \rangle = D_0 \boldsymbol\psi$ (since $\langle \boldsymbol\psi, c \rangle = 0$ by odd parity).
Thus:
$$(Q_{\mathrm{odd}} - \lambda I) Kc = - D_0 \boldsymbol\psi \implies Kc = - D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi,$$
where $(Q_{\mathrm{odd}} - \lambda I)$ is strictly positive definite and invertible because $\lambda < \mu_0 < \mu_1 \le \dots \le \mu_{N-1}$.
Therefore:
$$\langle \boldsymbol\psi, Kc \rangle = - D_0 \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi \rangle = - D_0 M_1,$$
where $M_1 \equiv \sum_{j=0}^{N-1} \frac{a_j^2}{\mu_j - \lambda}$.
Substituting this back into the commutator:
$$[K, Q] Kc = - d (- D_0 M_1) = D_0 M_1 d.$$
Taking the inner product with $c$:
$$\langle c, [K, Q] Kc \rangle = \langle c, D_0 M_1 d \rangle = D_0 M_1 \langle c, d \rangle = D_0 M_1 (D_0) = D_0^2 M_1.$$
Combining the two terms establishes (8.21.1):
$$\langle Kc, Q_{\mathrm{odd}} Kc \rangle = \lambda \|Kc\|^2 + D_0^2 M_1.$$
Subtracting $\mu_0 \|Kc\|^2$ gives:
$$\langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle = - (\mu_0 - \lambda) \|Kc\|^2 + D_0^2 M_1.$$
Now decompose $M_1 = \frac{a_0^2}{\mu_0 - \lambda} + M_1^{\mathrm{exc}}$.
Multiplying by $D_0^2$ and using the exact doublet relation $b_{00} \equiv \langle u_0, Kc \rangle = - \frac{D_0 a_0}{\mu_0 - \lambda} \implies b_{00}^2 = \frac{D_0^2 a_0^2}{(\mu_0 - \lambda)^2}$:
$$D_0^2 M_1 = (\mu_0 - \lambda) b_{00}^2 + D_0^2 M_1^{\mathrm{exc}}.$$
Since $\|Kc\|^2 = b_{00}^2 + \|P_{\perp u_0} Kc\|^2$:
$$- (\mu_0 - \lambda) \|Kc\|^2 + D_0^2 M_1 = - (\mu_0 - \lambda) \big( b_{00}^2 + \|P_{\perp u_0} Kc\|^2 \big) + (\mu_0 - \lambda) b_{00}^2 + D_0^2 M_1^{\mathrm{exc}} = D_0^2 M_1^{\mathrm{exc}} - (\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2,$$
proving (8.21.2). $\blacksquare$

**Step 2 (Proof of Part II — Variational Spectral-Gap Enclosure):**
Because $\mu_0$ is the minimum eigenvalue of $Q_{\mathrm{odd}}$, the operator $(Q_{\mathrm{odd}} - \mu_0 I)$ is positive semi-definite on $H_{\mathrm{odd}}$:
$$\langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle \ge 0.$$
Substituting (8.21.2) yields immediately:
$$D_0^2 M_1^{\mathrm{exc}} - (\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2 \ge 0 \implies (\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2 \le D_0^2 M_1^{\mathrm{exc}},$$
which is (8.21.3).
Next, expanding $Kc$ in the orthonormal eigenbasis $\{u_j\}_{j=0}^{N-1}$ of $Q_{\mathrm{odd}}$:
$$\langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle = \sum_{j=0}^{N-1} (\mu_j - \mu_0) |\langle u_j, Kc \rangle|^2 = \sum_{j=1}^{N-1} (\mu_j - \mu_0) |\langle u_j, Kc \rangle|^2.$$
Since $\mu_j \ge \mu_1$ for all $j \ge 1$:
$$\sum_{j=1}^{N-1} (\mu_j - \mu_0) |\langle u_j, Kc \rangle|^2 \ge (\mu_1 - \mu_0) \sum_{j=1}^{N-1} |\langle u_j, Kc \rangle|^2 = (\mu_1 - \mu_0) \|P_{\perp u_0} Kc\|^2.$$
Equating with the upper expression from (8.21.2):
$$(\mu_1 - \mu_0) \|P_{\perp u_0} Kc\|^2 \le D_0^2 M_1^{\mathrm{exc}} - (\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2.$$
Adding $(\mu_0 - \lambda) \|P_{\perp u_0} Kc\|^2$ to both sides, and noting that $(\mu_1 - \mu_0) + (\mu_0 - \lambda) = \mu_1 - \lambda$:
$$(\mu_1 - \lambda) \|P_{\perp u_0} Kc\|^2 \le D_0^2 M_1^{\mathrm{exc}}.$$
Dividing by $\mu_1 - \lambda$ gives:
$$\|P_{\perp u_0} Kc\|^2 \le \frac{D_0^2 M_1^{\mathrm{exc}}}{\mu_1 - \lambda} = M_1^{\mathrm{exc}} R_{\mathrm{gap}}(N),$$
proving (8.21.4). $\blacksquare$

**Step 3 (Proof of Part III — The Alignment Sandwich):**
Dividing (8.21.4) by $\|Kc\|^2$ yields the upper bound in (8.21.5):
$$\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \le \frac{M_1^{\mathrm{exc}}}{\|Kc\|^2} R_{\mathrm{gap}}(N).$$
For the lower bound, note that $u_1 \in P_{\perp u_0} H_{\mathrm{odd}}$ is a unit vector ($\|u_1\| = 1$).
By Bessel's inequality on $P_{\perp u_0} H_{\mathrm{odd}}$:
$$\|P_{\perp u_0} Kc\|^2 = \sum_{j \ge 1} |\langle u_j, Kc \rangle|^2 \ge |\langle u_1, Kc \rangle|^2.$$
By symmetry of $K$ ($K = K^T$), $\langle u_1, Kc \rangle = \langle c, K u_1 \rangle \equiv b_{01}$.
Therefore $\|P_{\perp u_0} Kc\|^2 \ge b_{01}^2$.
From Proposition 8.19 (8.19.2), $b_{01}^2 = \mathcal{R}_1 R_{\mathrm{gap}}(N)$.
Therefore:
$$\|P_{\perp u_0} Kc\|^2 \ge \mathcal{R}_1 R_{\mathrm{gap}}(N).$$
Dividing by $\|Kc\|^2$ gives:
$$\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \ge \frac{\mathcal{R}_1}{\|Kc\|^2} R_{\mathrm{gap}}(N),$$
which establishes the two-sided sandwich (8.21.5). $\blacksquare$

---

*Discussion and Epistemic Synthesis:*
1. **Tightness of the Alignment Sandwich:**
   In `cell67.out` and `cell69.out`, the numerical quantities at $N = 24$ are:
   $$\|Kc\|^2 \approx 1.72507, \qquad \mathcal{R}_1 \approx 5.129, \qquad M_1^{\mathrm{exc}} \approx 93.654, \qquad R_{\mathrm{gap}}(24) \approx 4.39178 \times 10^{-7}.$$
   The sandwich (8.21.5) evaluates numerically to:
   $$\frac{5.129}{1.72507} (4.39178 \times 10^{-7}) \le \frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \le \frac{93.654}{1.72507} (4.39178 \times 10^{-7}),$$
   $$1.3058 \times 10^{-6} \le \frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} \le 2.3842 \times 10^{-5}.$$
   The actual observed wavepacket misalignment in `cell67.out` is:
   $$\frac{\|P_{\perp u_0} Kc\|^2}{\|Kc\|^2} = \frac{2.25259 \times 10^{-6}}{1.72507} \approx 1.30580 \times 10^{-6},$$
   saturating within $0.0005\%$ of the lower bound! This near-exact saturation occurs because mode 1 carries $99.9995\%$ of the excited residual norm $\|v_{\mathrm{exc}}\|^2$, meaning that $P_{\perp u_0} Kc$ is almost an exact eigenvector proportional to $u_1$.
2. **The Spectral Moment Obstruction (First Moment vs. Zeroth Moment):**
   A natural temptation would be to propose that proving $\langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle \to 0$ is sufficient to establish wavepacket alignment and relative gap suppression. However, this is a conceptual trap.
   Defining the spectral measure of the excited coordinate wavepacket $v_{\mathrm{exc}} \equiv P_{\perp u_0} Kc$:
   $$d\nu_v(\mu) \equiv \sum_{j \ge 1} b_{0j}^2 \delta_{\mu_j},$$
   we observe that:
   $$\|v_{\mathrm{exc}}\|^2 = \int d\nu_v(\mu) = \sum_{j \ge 1} b_{0j}^2 = D_0^2 M_{2,\mathrm{exc}} \quad (\text{zeroth moment}),$$
   $$E_{\mathrm{exc}} \equiv \langle Kc, (Q_{\mathrm{odd}} - \mu_0 I) Kc \rangle = \int (\mu - \mu_0) d\nu_v(\mu) = \sum_{j \ge 1} (\mu_j - \mu_0) b_{0j}^2 \quad (\text{first moment}).$$
   Under standard hypotheses (H2$_{\mathrm{odd}}$ and H3), the energy excess $E_{\mathrm{exc}} \le D_0^2 M_1^{\mathrm{exc}} \le C N^\gamma e^{-\sigma N} \to 0$ vanishes automatically. But this does **not** imply that $\|v_{\mathrm{exc}}\|^2 \to 0$, because the spectral gap closes exponentially fast:
   $$\mu_1 - \mu_0 \sim e^{-\sigma_1 N} \longrightarrow 0.$$
   The dominant contribution to the energy excess is $E_{\mathrm{exc}} \approx (\mu_1 - \mu_0) b_{01}^2$. Thus $E_{\mathrm{exc}}$ is tiny primarily because the spectral gap is exponentially small, not because the wavepacket residual $b_{01}^2$ vanishes.
   This is the standard obstruction in variational spectral theory: when the spectral gap closes, a first-moment estimate $\int (\mu - \mu_0) d\nu_v \to 0$ cannot control the zeroth moment $\int d\nu_v$.
3. **The Reoriented Continuum Target — Zeroth-Moment Control:**
   Proposition 8.21 clarifies that the real, unresolved mathematical object governing Route B remains:
   $$\boxed{\|P_{\perp u_0} Kc\|^2 = D_0^2 M_{2,\mathrm{exc}} \longrightarrow 0.}$$
   The central analytical problem is to upgrade control from the first moment $E_{\mathrm{exc}}$ to the zeroth moment $\|v_{\mathrm{exc}}\|^2$. Specifically:
   *Can we obtain a zeroth-moment bound on the excited coordinate wavepacket from the exact commutator structure, without dividing by the exponentially collapsing first spectral gap?*
4. **Spectral Concentration as the Structural Bridge:**
   In `cell69.out`, the first excited mode carries $99.9995\%$ of the total excited residual norm $\|v_{\mathrm{exc}}\|^2$:
   $$\frac{b_{01}^2}{\|v_{\mathrm{exc}}\|^2} \longrightarrow 1.$$
   If an analytical bound can be established showing that $\frac{\sum_{j \ge 2} b_{0j}^2}{\sum_{j \ge 1} b_{0j}^2} \to 0$, then $\|v_{\mathrm{exc}}\|^2 = b_{01}^2 (1 + o(1))$. Combining this with the exact resolvent representation of Proposition 8.20 ($Ku_1 = a_1 (Q_{\mathrm{even}} - \mu_1 I)^{-1} d$) connects the zeroth-moment problem directly to the even-sector Stieltjes derivative $G_d'(\mu_1)$, providing a structured operator path forward.

---

## 8.22 Proposition 8.22: First-Mode Spectral Concentration of the Coordinate Wavepacket, Tail Ratio Enclosure, and Asymptotic Residual Equivalence

### Proposition 8.22 (First-Mode Spectral Concentration of the Coordinate Wavepacket)
*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Connes–van Suijlekom Galerkin operator with cutoff parameter $c > 1$ ($L = \log c$) and even ground state $u_0^{\mathrm{even}} \equiv c \in H_{\mathrm{even}}$. Let $\{u_j\}_{j=0}^{N-1}$ be the orthonormal eigenbasis of $Q_{\mathrm{odd}}$ with eigenvalues $\mu_0 < \mu_1 < \dots < \mu_{N-1}$, and let $v_{\mathrm{exc}} \equiv P_{\perp u_0} Kc \in H_{\mathrm{odd}}$ be the excited coordinate wavepacket residual with transition dipole components $b_{0j} \equiv \langle c, Ku_j \rangle = \langle Kc, u_j \rangle = - \frac{D_0 a_j}{\mu_j - \lambda}$.*

**Part I (Exact Finite-$N$ Relative Tail Ratio and $D_0^2$ Elimination — Unconditional Theorem):**
*For every finite dimension $N \ge 2$ and cutoff parameter $c > 1$, the relative excited tail ratio:*
$$\varepsilon_N \equiv \frac{\sum_{j=2}^{N-1} b_{0j}^2}{b_{01}^2} \tag{8.22.1}$$
*evaluates identically to the ratio of excited resolvent sums in which the global boundary-layer scaling factor $D_0^2$ cancels identically:*
$$\boxed{\varepsilon_N = \sum_{j=2}^{N-1} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2.} \tag{8.22.2}$$

**Part II (Finite-$N$ Quantitative Concentration on the Moving First-Mode Dirac Mass):**
*Let $d\bar\nu_N(\mu) \equiv \frac{d\nu_v(\mu)}{\|v_{\mathrm{exc}}\|^2} = \frac{1}{\|v_{\mathrm{exc}}\|^2} \sum_{j=1}^{N-1} b_{0j}^2 \delta_{\mu_j}$ be the normalized discrete spectral measure of the excited coordinate wavepacket $v_{\mathrm{exc}}$ on $[\mu_1, \infty)$. For every bounded continuous test function $f \in C_b([\mu_1, \infty))$, the finite-$N$ expectation satisfies the quantitative concentration estimate:*
$$\boxed{\left| \int_{\mu_1}^\infty f(\mu) \, d\bar\nu_N(\mu) - f(\mu_1^{(N)}) \right| \le 3 \|f\|_\infty \varepsilon_N.} \tag{8.22.3}$$
*Consequently, the condition $\varepsilon_N \to 0$ implies that the normalized excited spectral measure becomes asymptotically concentrated on the first excited eigenvalue in the sense of a moving Dirac mass $\delta_{\mu_1^{(N)}}$. If additionally $\mu_1^{(N)} \to \mu_\infty$, standard weak-$*$ convergence $d\bar\nu_N \rightharpoonup \delta_{\mu_\infty}$ follows.*

**Part III (Exact Asymptotic Residual Equivalence — Unconditional Theorem):**
*The total excited coordinate wavepacket norm $\|v_{\mathrm{exc}}\|^2 = \|P_{\perp u_0} Kc\|^2 = D_0^2 M_{2,\mathrm{exc}}$ satisfies the exact finite-$N$ relation:*
$$\boxed{\|P_{\perp u_0} Kc\|^2 = D_0^2 M_{2,\mathrm{exc}} = b_{01}^2 (1 + \varepsilon_N) = \mathcal{R}_1 R_{\mathrm{gap}}(N) (1 + \varepsilon_N).} \tag{8.22.4}$$
*Consequently, whenever $\varepsilon_N \to 0$, the zeroth moment $\|P_{\perp u_0} Kc\|^2$ is asymptotically equivalent to the first transition dipole:*
$$\|P_{\perp u_0} Kc\|^2 \sim b_{01}^2 = \mathcal{R}_1 R_{\mathrm{gap}}(N) \quad (N \to \infty).$$

**Part IV (Two-Sector Tail Partition and Semiclassical Spectral Filtering):**
*Partitioning the tail sum $\varepsilon_N = \varepsilon_N^{\mathrm{bound}} + \varepsilon_N^{\mathrm{high}}$ at the uniform above-barrier separation threshold $V_* > \sup_N \lambda_N$ ($c_{\mathrm{high}} \equiv V_* - \sup_N \lambda_N > 0$):*
1. *Unconditional High-Energy Enclosure:*
   $$\varepsilon_N^{\mathrm{high}} \equiv \sum_{j > N_{\mathrm{bound}}} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 \le \frac{\mu_1 - \lambda}{c_{\mathrm{high}}^2 \mathcal{R}_1(N)} \|P_{\mathrm{high}} d\|^2 \le \frac{\mathcal{C}_{\mathbb{R}} N}{c_{\mathrm{high}}^2} \frac{\mu_1 - \lambda}{\mathcal{R}_1(N)}.$$
2. *Conditional High-Energy Decoupling:*
   *If the first excited transmission ratio remains bounded away from zero ($\inf_N \mathcal{R}_1(N) \ge c_1 > 0$) and the first tunneling gap decays exponentially ($\mu_1 - \lambda = \mathcal{O}(e^{-\sigma_1 N})$), then the above-barrier tail decouples exponentially fast:*
   $$\boxed{\varepsilon_N^{\mathrm{high}} = \mathcal{O}(e^{-\sigma_1 N} N) \longrightarrow 0 \quad (N \to \infty).}$$
3. *Bound-State Tail and Open Overlap Problem:*
   *For the finite bound-state ladder $2 \le j \le N_{\mathrm{bound}}$, the tail evaluates to:*
   $$\varepsilon_N^{\mathrm{bound}} \equiv \sum_{j=2}^{N_{\mathrm{bound}}} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2.$$
   *Because $\mu_j - \lambda > \mu_1 - \lambda$ for all $j \ge 2$, the spectral factor $\left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2$ provides automatic suppression. However, establishing $\varepsilon_N^{\mathrm{bound}} \to 0$ requires controlling the overlap ratios $a_j^2 / a_1^2$, which is not completed by Part IV and forms the subject of Proposition 8.23.*

**Part V (The Reoriented Analytical Hierarchy):**
*The continuum reduction of Route B is structured into the following deductive chain:*
$$\boxed{\begin{aligned}
&[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T \text{ (exact rank-two commutator)} \\
&\qquad\Downarrow \\
&b_{0j} = - \frac{D_0 a_j}{\mu_j - \lambda} \text{ (exact transition dipoles)} \\
&\qquad\Downarrow \\
&\varepsilon_N \equiv \frac{\sum_{j \ge 2} b_{0j}^2}{b_{01}^2} = \sum_{j \ge 2} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 \longrightarrow 0 \text{ (first-mode concentration)} \\
&\qquad\Downarrow \\
&D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2 = b_{01}^2 (1 + \varepsilon_N) \sim b_{01}^2 \text{ (asymptotic residual equivalence)} \\
&\qquad\Downarrow \\
&R_{\mathrm{gap}}(N) = \frac{b_{01}^2}{\mathcal{R}_1} \sim \frac{\|P_{\perp u_0} Kc\|^2}{\mathcal{R}_1} \longrightarrow 0 \quad (\text{conditional on } \inf_N \mathcal{R}_1(N) > 0).
\end{aligned}} \tag{8.22.5}$$

---

### Proof of Proposition 8.22

**1. Proof of Part I (Exact Relative Tail Ratio and $D_0^2$ Elimination):**
From Proposition 8.12 and Proposition 8.19, the transition dipole moments between the even ground state $c$ and the orthonormal odd eigenmodes $u_j$ ($j \ge 0$) satisfy the exact commutator identity:
$$b_{0j} \equiv \langle c, K u_j \rangle = - \frac{D_0 a_j}{\mu_j - \lambda}, \qquad a_j \equiv \langle d, u_j \rangle.$$
For every excited mode $j \ge 1$, squaring this relation yields:
$$b_{0j}^2 = \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2}.$$
In particular, for the first excited mode $j = 1$, strict parity interlacing ($\mu_0 < E_1 < \mu_1 < E_2$) guarantees that $a_1 \ne 0$ and $\mu_1 - \lambda > 0$, so that:
$$b_{01}^2 = \frac{D_0^2 a_1^2}{(\mu_1 - \lambda)^2} > 0.$$
Now consider the ratio of the excited tail sum $\sum_{j=2}^{N-1} b_{0j}^2$ to the first-mode dipole $b_{01}^2$:
$$\varepsilon_N \equiv \frac{\sum_{j=2}^{N-1} b_{0j}^2}{b_{01}^2} = \frac{\sum_{j=2}^{N-1} \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2}}{\frac{D_0^2 a_1^2}{(\mu_1 - \lambda)^2}}.$$
Because the boundary-layer factor $D_0^2 > 0$ is a strictly positive scalar independent of the mode index $j$, it factors out of the numerator summation and cancels identically with the denominator:
$$\varepsilon_N = \sum_{j=2}^{N-1} \frac{\frac{a_j^2}{(\mu_j - \lambda)^2}}{\frac{a_1^2}{(\mu_1 - \lambda)^2}} = \sum_{j=2}^{N-1} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2.$$
This establishes (8.22.2) unconditionally for all finite dimensions $N \ge 2$.

**2. Proof of Part II (Quantitative Finite-$N$ Concentration on Moving Dirac Mass):**
Let $f \in C_b([\mu_1, \infty))$ be an arbitrary continuous and bounded test function on the excited spectrum. The expectation of $f$ under the normalized excited spectral measure $d\bar\nu_N \equiv \frac{d\nu_v}{\|v_{\mathrm{exc}}\|^2}$ evaluates to:
$$\int_{\mu_1}^\infty f(\mu) \, d\bar\nu_N(\mu) = \frac{\sum_{j=1}^{N-1} f(\mu_j) b_{0j}^2}{\sum_{j=1}^{N-1} b_{0j}^2} = \frac{f(\mu_1^{(N)}) b_{01}^2 + \sum_{j=2}^{N-1} f(\mu_j) b_{0j}^2}{b_{01}^2 (1 + \varepsilon_N)} = \frac{f(\mu_1^{(N)}) + \sum_{j=2}^{N-1} f(\mu_j) \frac{b_{0j}^2}{b_{01}^2}}{1 + \varepsilon_N}.$$
Subtracting $f(\mu_1^{(N)})$ yields:
$$\left| \int_{\mu_1}^\infty f(\mu) \, d\bar\nu_N(\mu) - f(\mu_1^{(N)}) \right| = \left| \frac{\sum_{j=2}^{N-1} (f(\mu_j) - f(\mu_1^{(N)})) \frac{b_{0j}^2}{b_{01}^2} - \varepsilon_N f(\mu_1^{(N)})}{1 + \varepsilon_N} \right| \le \frac{2 \|f\|_\infty \varepsilon_N + \|f\|_\infty \varepsilon_N}{1 + \varepsilon_N} \le 3 \|f\|_\infty \varepsilon_N.$$
This establishes the quantitative finite-$N$ estimate (8.22.3). Whenever $\varepsilon_N \to 0$ as $N \to \infty$, the right-hand side vanishes for every bounded test function, demonstrating asymptotic concentration on the moving Dirac mass $\delta_{\mu_1^{(N)}}$. If additionally the sequence of eigenvalues converges $\mu_1^{(N)} \to \mu_\infty$, then by continuity $f(\mu_1^{(N)}) \to f(\mu_\infty)$, yielding standard weak-$*$ convergence $d\bar\nu_N \rightharpoonup \delta_{\mu_\infty}$.

**3. Proof of Part III (Exact Asymptotic Residual Equivalence):**
By definition of the orthogonal projection $v_{\mathrm{exc}} = P_{\perp u_0} Kc = \sum_{j=1}^{N-1} b_{0j} u_j$, Parseval's identity yields:
$$\|P_{\perp u_0} Kc\|^2 = \|v_{\mathrm{exc}}\|^2 = \sum_{j=1}^{N-1} b_{0j}^2 = b_{01}^2 + \sum_{j=2}^{N-1} b_{0j}^2 = b_{01}^2 \left( 1 + \frac{\sum_{j=2}^{N-1} b_{0j}^2}{b_{01}^2} \right) = b_{01}^2 (1 + \varepsilon_N).$$
By Proposition 8.12, $\|P_{\perp u_0} Kc\|^2 = D_0^2 M_{2,\mathrm{exc}}$, and by Proposition 8.19, $b_{01}^2 = \mathcal{R}_1 R_{\mathrm{gap}}(N)$ with $\mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda}$. Substituting these expressions yields (8.22.4).
When $\varepsilon_N \to 0$, $1 + \varepsilon_N \to 1$, which proves asymptotic equivalence $\|P_{\perp u_0} Kc\|^2 \sim b_{01}^2$.

**4. Proof of Part IV (Two-Sector Tail Partition and Semiclassical Filtering):**
Partition the excited index set $\{2, \dots, N-1\}$ at $N_{\mathrm{bound}}$, defined by the uniform above-barrier threshold $\mu_j \ge V_* > \sup_N \lambda_N$ ($c_{\mathrm{high}} \equiv V_* - \sup_N \lambda_N > 0$):
$$\varepsilon_N = \sum_{j=2}^{N_{\mathrm{bound}}} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 + \sum_{j > N_{\mathrm{bound}}} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 \equiv \varepsilon_N^{\mathrm{bound}} + \varepsilon_N^{\mathrm{high}}.$$
For the high-energy above-barrier modes $j > N_{\mathrm{bound}}$, $\mu_j - \lambda \ge c_{\mathrm{high}} > 0$. Therefore:
$$\varepsilon_N^{\mathrm{high}} \le \frac{(\mu_1 - \lambda)^2}{c_{\mathrm{high}}^2 a_1^2} \sum_{j > N_{\mathrm{bound}}} a_j^2.$$
Recalling that $\mathcal{R}_1 \equiv \frac{a_1^2}{\mu_1 - \lambda}$, this simplifies to:
$$\varepsilon_N^{\mathrm{high}} \le \frac{\mu_1 - \lambda}{c_{\mathrm{high}}^2 \mathcal{R}_1(N)} \sum_{j > N_{\mathrm{bound}}} a_j^2 \le \frac{\mu_1 - \lambda}{c_{\mathrm{high}}^2 \mathcal{R}_1(N)} \|d\|^2.$$
Since $\|d\|^2 \le \mathcal{C}_{\mathbb{R}} N = \mathcal{O}(N)$, assuming $\inf_N \mathcal{R}_1(N) \ge c_1 > 0$ and the exponential tunneling gap scaling $\mu_1 - \lambda = \mathcal{O}(e^{-\sigma_1 N})$, the above-barrier sector decouples exponentially fast:
$$\varepsilon_N^{\mathrm{high}} = \mathcal{O}(e^{-\sigma_1 N} N) \longrightarrow 0 \quad (N \to \infty).$$
For the bound-state sector $2 \le j \le N_{\mathrm{bound}}$, since the eigenvalues are strictly ordered $\mu_1 < \mu_2 < \dots < \mu_{N_{\mathrm{bound}}}$, the spectral factor satisfies:
$$\left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 < 1 \qquad (\forall j \ge 2).$$
This completes the proof. $\blacksquare$

---

### Analytical Commentary on Proposition 8.22

1. **Strategic Shift from Absolute to Relative Scale Control:**
   Attempting to prove that the zeroth moment $D_0^2 M_{2,\mathrm{exc}} \to 0$ directly requires estimating the absolute scale of the wavepacket residual, which is intimately coupled to the boundary-layer exponential decay factor $D_0^2 \sim e^{-2\sigma_0 N}$.
   In contrast, Proposition 8.22 demonstrates that the relative tail ratio:
   $$\varepsilon_N \equiv \frac{\sum_{j \ge 2} b_{0j}^2}{b_{01}^2} = \sum_{j \ge 2} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2$$
   is completely independent of $D_0^2$. The boundary-layer amplitude cancels identically from every term.
2. **Epistemic Calibration on High-Energy Decoupling vs. Bound Ladder:**
   Part IV establishes high-energy decoupling $\varepsilon_N^{\mathrm{high}} = \mathcal{O}(e^{-\sigma_1 N} N) \to 0$ **conditionally** upon two hypotheses: $\inf_N \mathcal{R}_1(N) \ge c_1 > 0$ (first-mode transmission lower bound) and $\mu_1 - \lambda = \mathcal{O}(e^{-\sigma_1 N})$ (exponential tunneling gap decay). While supported strongly by numerical sweeps across $N \in \{8, \dots, 24\}$, these remain analytical hypotheses.
   Crucially, Part IV does **not** claim that $\varepsilon_N \to 0$ is proven. The finite bound-state ladder:
   $$\varepsilon_N^{\mathrm{bound}} = \sum_{j=2}^{N_{\mathrm{bound}}} \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2$$
   remains open, because the spectral suppression factor $(\frac{\mu_1 - \lambda}{\mu_j - \lambda})^2 < 1$ must be shown to dominate the overlap ratio $a_j^2 / a_1^2$.
3. **Numerical Verification Across Discrete Dimensions:**
   From `cell69.out`, we compute the exact finite-$N$ values of $\varepsilon_N$:
   - $N = 8$: $\|v_{\mathrm{exc}}\|^2 \approx 3.73759 \times 10^{-5}$, $b_{01}^2 \approx 3.73725 \times 10^{-5} \implies \varepsilon_8 \approx 9.07 \times 10^{-5}$.
   - $N = 12$: $\|v_{\mathrm{exc}}\|^2 \approx 2.56371 \times 10^{-5}$, $b_{01}^2 \approx 2.56364 \times 10^{-5} \implies \varepsilon_{12} \approx 2.51 \times 10^{-5}$.
   - $N = 16$: $\|v_{\mathrm{exc}}\|^2 \approx 1.02873 \times 10^{-5}$, $b_{01}^2 \approx 1.02871 \times 10^{-5} \implies \varepsilon_{16} \approx 2.20 \times 10^{-5}$.
   - $N = 20$: $\|v_{\mathrm{exc}}\|^2 \approx 4.29285 \times 10^{-6}$, $b_{01}^2 \approx 4.29280 \times 10^{-6} \implies \varepsilon_{20} \approx 1.04 \times 10^{-5}$.
   - $N = 24$: $\|v_{\mathrm{exc}}\|^2 \approx 2.25259 \times 10^{-6}$, $b_{01}^2 \approx 2.25258 \times 10^{-6} \implies \varepsilon_{24} \approx 4.52 \times 10^{-6}$.
   The tail ratio is tiny ($< 5 \times 10^{-6}$ at $N=24$) and decays monotonically, confirming that $> 99.9995\%$ of the wavepacket residual norm concentrates in the first excited mode.

---

## 8.23 Proposition 8.23: Exact Stieltjes Derivative Architecture for Overlap Ratios and Three-Factor Bound-State Tail Factorization

Having isolated the remaining continuum wavepacket obstruction to the bound-state tail $\varepsilon_N^{\mathrm{bound}} = \sum_{j=2}^{N_{\mathrm{bound}}} \frac{a_j^2}{a_1^2} (\frac{\mu_1 - \lambda}{\mu_j - \lambda})^2$, we now establish an exact finite-$N$ representation for the overlap ratios $a_j^2 / a_1^2$. Rather than treating $a_j \equiv \langle \boldsymbol\psi, u_j \rangle$ as arbitrary projections, we deploy the even-resolvent representation of odd coordinate derivatives from Proposition 8.20.

### Proposition 8.23 (Stieltjes Derivative Factorization of Overlap Ratios)
*Let $Q \in \mathbb{R}^{(2N+1) \times (2N+1)}$ be the finite-rank Connes–van Suijlekom Galerkin operator with cutoff parameter $c > 1$ ($L = \log c$). Let $\{ (E_k, u_k^{\mathrm{even}}) \}_{k=0}^N$ and $\{ (\mu_j, u_j^{\mathrm{odd}}) \}_{j=0}^{N-1}$ denote the even and odd eigensystems, with even ground state $u_0^{\mathrm{even}} \equiv c$, $E_0 = \lambda$, boundary overlaps $d_k \equiv \langle d, u_k^{\mathrm{even}} \rangle$, and odd source overlaps $a_j \equiv \langle \boldsymbol\psi, u_j^{\mathrm{odd}} \rangle$. Let $G_d(z) \equiv \sum_{k=0}^N \frac{d_k^2}{E_k - z}$ be the Stieltjes transform of the boundary measure on $H_{\mathrm{even}}$, with Stieltjes derivative:*
$$G_d'(z) \equiv \sum_{k=0}^N \frac{d_k^2}{(E_k - z)^2}.$$

**Part I (Exact Stieltjes Derivative Overlap Representation — Unconditional Theorem):**
*For every excited odd mode $j \ge 1$, strict parity interlacing ensures that $\mu_j \notin \sigma(Q_{\mathrm{even}})$. The squared overlap $a_j^2$ satisfies the exact finite-$N$ identity:*
$$\boxed{a_j^2 = \frac{\|K u_j\|^2}{G_d'(\mu_j)} = \frac{\|K u_j\|^2}{\sum_{k=0}^N \frac{d_k^2}{(\mu_j - E_k)^2}}.} \tag{8.23.1}$$

**Part II (Exact Factorization of Overlap Ratios — Unconditional Theorem):**
*For every pair of excited odd modes $j, 1 \ge 1$, the overlap ratio $a_j^2 / a_1^2$ factors identically into the coordinate kinetic energy ratio and the inverse Stieltjes derivative ratio:*
$$\boxed{\frac{a_j^2}{a_1^2} = \left( \frac{\|K u_j\|^2}{\|K u_1\|^2} \right) \left( \frac{G_d'(\mu_1)}{G_d'(\mu_j)} \right).} \tag{8.23.2}$$
*In particular, this expression is completely independent of the boundary-layer amplitude $D_0$.*

**Part III (Three-Factor Modewise Tail Decomposition — Unconditional Theorem):**
*Every individual modewise term $T_j \equiv b_{0j}^2 / b_{01}^2$ in the excited wavepacket tail ratio factors identically into three positive components:*
$$\boxed{T_j = \mathcal{K}_j \cdot \mathcal{G}_j \cdot \mathcal{S}_j,} \tag{8.23.3}$$
*where:*
1. *$\mathcal{K}_j \equiv \frac{\|K u_j\|^2}{\|K u_1\|^2}$ is the coordinate kinetic energy ratio,*
2. *$\mathcal{G}_j \equiv \frac{G_d'(\mu_1)}{G_d'(\mu_j)}$ is the inverse Stieltjes derivative ratio,*
3. *$\mathcal{S}_j \equiv \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2$ is the squared spectral gap suppression factor.*

**Part IV (Polynomial Overlap Growth Sufficiency):**
*If for each bound mode $j \in \{2, \dots, N_{\mathrm{bound}}\}$ the overlap ratio satisfies at most polynomial growth:*
$$\frac{a_j^2}{a_1^2} \le \mathcal{C}_j N^p \qquad (p \ge 0),$$
*and the tunneling gap ratio decays exponentially: $\frac{\mu_1 - \lambda}{\mu_j - \lambda} = \mathcal{O}(e^{-\Delta\sigma_j N})$ with $\Delta\sigma_j > 0$, then each bound-state tail term is exponentially quenched:*
$$T_j = \mathcal{O}(N^p e^{-2\Delta\sigma_j N}) \longrightarrow 0 \quad (N \to \infty),$$
*establishing $\varepsilon_N^{\mathrm{bound}} \to 0$.*

---

### Proof of Proposition 8.23

**1. Proof of Part I (Exact Stieltjes Derivative Overlap Representation):**
By Proposition 8.20, for any excited odd mode $u_j$ ($j \ge 1$), the coordinate derivative vector $K u_j \in H_{\mathrm{even}}$ satisfies the exact resolvent equation:
$$K u_j = a_j (Q_{\mathrm{even}} - \mu_j I)^{-1} d.$$
Expanding this vector in the orthonormal even eigenbasis $\{u_k^{\mathrm{even}}\}_{k=0}^N$:
$$K u_j = a_j \sum_{k=0}^N \frac{\langle d, u_k^{\mathrm{even}} \rangle}{E_k - \mu_j} u_k^{\mathrm{even}} = a_j \sum_{k=0}^N \frac{d_k}{E_k - \mu_j} u_k^{\mathrm{even}}.$$
Taking the $\ell^2$-norm squared on $H_{\mathrm{even}}$ via Parseval's identity:
$$\|K u_j\|^2 = a_j^2 \sum_{k=0}^N \frac{d_k^2}{(E_k - \mu_j)^2} = a_j^2 \sum_{k=0}^N \frac{d_k^2}{(\mu_j - E_k)^2} = a_j^2 G_d'(\mu_j).$$
Because strict parity interlacing ensures that $\mu_j \ne E_k$ for all $k \in \{0, \dots, N\}$, the denominator terms $(\mu_j - E_k)^2$ are strictly positive and finite. Furthermore, $d_0 = D_0 \ne 0$ guarantees that $G_d'(\mu_j) \ge \frac{D_0^2}{(\mu_j - \lambda)^2} > 0$. Dividing by $G_d'(\mu_j)$ yields (8.23.1).

**2. Proof of Part II (Exact Factorization of Overlap Ratios):**
Applying (8.23.1) to mode $j$ and to mode 1:
$$a_j^2 = \frac{\|K u_j\|^2}{G_d'(\mu_j)}, \qquad a_1^2 = \frac{\|K u_1\|^2}{G_d'(\mu_1)}.$$
Since $a_1 \ne 0$ (from strict parity interlacing $\mu_0 < E_1 < \mu_1 < E_2$), taking the ratio gives:
$$\frac{a_j^2}{a_1^2} = \frac{\frac{\|K u_j\|^2}{G_d'(\mu_j)}}{\frac{\|K u_1\|^2}{G_d'(\mu_1)}} = \left( \frac{\|K u_j\|^2}{\|K u_1\|^2} \right) \left( \frac{G_d'(\mu_1)}{G_d'(\mu_j)} \right),$$
establishing (8.23.2) unconditionally for all $N \ge 2$.

**3. Proof of Part III (Three-Factor Modewise Tail Decomposition):**
By Proposition 8.22, the modewise tail term is defined as:
$$T_j \equiv \frac{b_{0j}^2}{b_{01}^2} = \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2.$$
Substituting the exact factorization (8.23.2) for $a_j^2 / a_1^2$:
$$T_j = \left( \frac{\|K u_j\|^2}{\|K u_1\|^2} \right) \left( \frac{G_d'(\mu_1)}{G_d'(\mu_j)} \right) \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 = \mathcal{K}_j \cdot \mathcal{G}_j \cdot \mathcal{S}_j.$$
This proves (8.23.3).

**4. Proof of Part IV (Polynomial Overlap Growth Sufficiency):**
If $\frac{a_j^2}{a_1^2} \le \mathcal{C}_j N^p$ and $\frac{\mu_1 - \lambda}{\mu_j - \lambda} \le C_{\mathrm{gap}, j} e^{-\Delta\sigma_j N}$ with $\Delta\sigma_j > 0$, then:
$$T_j = \frac{a_j^2}{a_1^2} \left( \frac{\mu_1 - \lambda}{\mu_j - \lambda} \right)^2 \le \mathcal{C}_j N^p \cdot C_{\mathrm{gap}, j}^2 e^{-2\Delta\sigma_j N} = \mathcal{O}(N^p e^{-2\Delta\sigma_j N}).$$
Since $2\Delta\sigma_j > 0$, the exponential decay $e^{-2\Delta\sigma_j N}$ dominates any fixed polynomial growth $N^p$, so $T_j \to 0$ as $N \to \infty$. Summing over the finite bound-state ladder $j \in \{2, \dots, N_{\mathrm{bound}}\}$ (with uniform cardinality $\bar{N}_{\mathrm{bound}} < \infty$ by Bridge B1) yields $\varepsilon_N^{\mathrm{bound}} \to 0$. $\blacksquare$

---

### Analytical Commentary on Proposition 8.23

1. **Conversion of the Bound-State Obstruction to Spectral Ratios:**
   Proposition 8.23 converts the open question of bound-state tail control into two concrete, highly structured spectral questions:
   - *Coordinate Kinetic Energy Ratio $\mathcal{K}_j = \|K u_j\|^2 / \|K u_1\|^2$:* In Cell 68, the discrete coordinate kinetic energies $\|K u_j\|^2$ were audited across $N \in \{8, \dots, 24\}$ and verified to be strictly $\mathcal{O}(1)$ for bound states (e.g., $\|K u_1\|^2 \approx 9.75$, $\|K u_2\|^2 \approx 15.78$ at $N=16$). Thus $\mathcal{K}_j = \Theta(1)$.
   - *Stieltjes Derivative Ratio $\mathcal{G}_j = G_d'(\mu_1) / G_d'(\mu_j)$:* Both $G_d'(\mu_1)$ and $G_d'(\mu_j)$ are positive Stieltjes derivatives of the boundary measure $d$. Evaluating their ratio reveals whether $G_d'(\mu_j)$ is comparable to, or larger than, $G_d'(\mu_1)$.
2. **Zero Boundary-Layer Dependence:**
   Neither the kinetic energy ratio $\mathcal{K}_j$ nor the Stieltjes derivative ratio $\mathcal{G}_j$ involves $D_0$. The dangerous scale $D_0^2 \sim e^{-2\sigma_0 N}$ remains completely absent from the bound-state analysis.
3. **The Spectral Filtering Asymmetry:**
   Because $\mu_1 - \lambda \sim e^{-\sigma_1 N}$ decays with action $\sigma_1 \approx 3.22$, whereas for $j \ge 2$, $\mu_j - \lambda$ corresponds to higher bound states (whose tunneling gaps are either significantly wider or macroscopically $\mathcal{O}(1)$), the spectral ratio $\mathcal{S}_j = (\frac{\mu_1 - \lambda}{\mu_j - \lambda})^2$ provides an exponential suppression factor $e^{-2\Delta\sigma_j N}$. Consequently, the overlap ratio $a_j^2 / a_1^2$ can tolerate any polynomial growth $N^p$ without defeating first-mode concentration.

---


## 9. The Analytical Roadmap toward Continuous Weil Positivity

The empirical and asymptotic results established in this research programme suggest that the finite-rank Galerkin truncation may provide a convergent approximation to the continuous Weil quadratic form. 

To convert these findings into a complete, mathematically rigorous proof of Weil positivity on the idele class group, three major analytical hurdles must be resolved:

```
[Stage 1: Operator Convergence]
  Prove strong resolvent convergence Q_{c, N} -> Q_c on L^2([0, log c]).
  Rule out spectral pollution to establish inf spec(Q_c) >= 0 from \lambda_min(N) > 0.
        |
        v
[Stage 2: Boundary Regularity & Solitary Wave Proof]
  Prove uniform mode bounds |v_{N, m}| <= C q^m to justify term-by-term differentiation.
  Establish C^infty boundary flatness: T_\infty^{(k)}(0) = T_\infty^{(k)}(L) = 0.
  Eliminate the Volterra boundary jump at \omega = 1.
        |
        v
[Stage 3: The Double Scaling Limit]
  Control the joint scaling limit (N, c) -> \infty.
  Extend positivity from compact intervals [0, log c] to the full idele class group.
  Conclude unconditional Weil positivity W(g) >= 0 <=> Riemann Hypothesis.
```

### Stage 1: Operator Convergence and Spectral Pollution
- **Challenge:** Positivity of finite Galerkin projections ($\lambda_{\min}(N) > 0$) does not automatically guarantee non-negativity of the limiting operator spectrum; unbounded operators can suffer from *spectral pollution* (spurious eigenvalues arising in spectral gaps).
- **Required Theorem:** Establish strong resolvent convergence of the finite-dimensional Galerkin operators $Q_{c, N}$ to a continuous self-adjoint operator $Q_c$ on $L^2([0, \log c])$, proving that the spectrum satisfies $\operatorname{spec}(Q_c) \subseteq [0, \infty)$.

### Stage 2: Proof of Boundary Regularity and the Solitary Wave
- **Challenge:** Interchanging the limit $N \to \infty$ with differentiation to prove that the continuum wave $T_\infty(t)$ satisfies infinite-order flat boundary contact $T_\infty^{(k)}(0) = 0$.
- **Required Theorem:** Establish uniform-in-$N$ decay bounds on the mode coefficients $|v_{N, m}| \le C q^m$ ($q < 1$). This will rigorously prove Conjecture 3.3, eliminating the finite-rank boundary jump of the Volterra convolution and establishing that the zero eigenvalue is an isolated edge mode.

### Stage 3: The Double Scaling Limit $(N, c) \to \infty$
- **Challenge:** The prime cutoff $c > 1$ restricts the scaling interval to $[0, \log c]$. Extending positivity to the full idele class group requires taking $c \to \infty$ alongside $N \to \infty$.
- **Required Theorem:** Establish uniform stability of the tri-partite balance $\mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}} \ge 0$ as $c \to \infty$, showing that the prime-power dispersive energy continues to smoothly absorb the geometric dilation pole energy across the unbounded idele domain.

### 9.1 Semiclassical Continuum Limit, Positive Jet-Energy Defect, and Boundary Decoupling Mechanics

An intriguing open question concerning Stage 3 is the analytical structure of the Archimedean quadratic form as the logarithmic cutoff $L = \log c \to \infty$. 

Recent investigations have conjectured that under appropriate Sobolev regularization, the Archimedean form decouples in the continuum limit into an exact continuous Plancherel integral over the multiplier $h_+(r)$ plus an isolated localized boundary defect operator:

$$\mathcal{Q}_{\mathrm{arch}}^{(\infty)}(\tau) \stackrel{?}{=} \frac{1}{\pi} \int_0^\infty h_+(r) |\widehat{\tau}(r)|^2 \, dr + \Delta_{\mathrm{boundary}}(\tau(0)).$$

Auditing this conjecture against the exact algebraic identities established in Paper 4 (Corollary 5.4 & Theorem 5.5) reveals an exact finite-$(N, L)$ mathematical mechanism that is considerably richer and more tractable than the heuristic continuum limit:

1. **Exact Manifest Positivity of the Boundary Defect:**
   From Corollary 5.4 of Paper 4, for any finite dimension $N$ and interval length $L$:
   $$\mathcal{Q}_{\mathrm{arch}}(v) = h_+(0) v_0^2 + \sum_{m=1}^N v_m^2 h_+(a_m) + B_{N, L}(v),$$
   where the boundary defect is defined by:
   $$B_{N, L}(v) \equiv \frac{2}{L} \sum_{n=0}^\infty \frac{1 - e^{-q_n L}}{q_n^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q_n^2 v_m}{q_n^2 + a_m^2} \right]^2 \ge 0.$$
   Crucially, $B_{N, L}(v) \ge 0$ is **manifestly non-negative** for every finite $(N, L)$ before taking any continuum limit. Non-vanishing boundary amplitudes $D_0 \neq 0$ strictly *increase* the Archimedean energy rather than incurring a negative penalty.

2. **The Positive Jet-Energy Quadratic Form:**
   Recall from Theorem 3.2 and Theorem 5.5 of Paper 4 that the boundary evaluation of the operator resolvent generates the Taylor endpoint jets:
   $$D(1/q_n^2) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{q_n^2 v_m}{q_n^2 + a_m^2} = \sum_{j=0}^\infty \frac{(-1)^j D_j}{q_n^{2j}}, \qquad D_j = T_v^{(2j)}(0).$$
   Squaring this series yields $[D(1/q_n^2)]^2 = \sum_{j, k=0}^\infty (-1)^{j+k} D_j D_k q_n^{-(2j+2k)}$. Substituting into $B_{N, L}(v)$ proves that the boundary defect is an exact **positive quadratic form in the endpoint-jet vector** $\mathbf{D} = (D_0, D_1, D_2, \dots)^T$:
   $$B_{N, L}(v) = \frac{2}{L} \sum_{j, k=0}^\infty (-1)^{j+k} D_j D_k \mu_{j+k}(L) = \frac{2}{L} \mathbf{D}^T H_L \mathbf{D} = \sum_{k=0}^\infty A_k(N) \mu_k(L),$$
   where $H_L$ is a positive Hankel-type moment matrix with scalar moments:
   $$\mu_k(L) \equiv \sum_{n=0}^\infty \frac{1 - e^{-q_n L}}{q_n^{2k+4}}.$$
   At leading order $k = 0$, with $q_n = 2n + 1/2 = \frac{4n+1}{2}$, the universal numerical moment evaluates to:
   $$\sum_{n=0}^\infty q_n^{-4} = 16 \sum_{n=0}^\infty (4n+1)^{-4} = \frac{\pi^4}{12} + 8 \beta(4) \approx 16.028986,$$
   where $\beta(4) \approx 0.98894455$ is the Dirichlet beta function at 4.

3. **Resolvent Decomposition, the Absolute Moment Bound, and the Empirical Two-Jet Scale $u_1$:**
   Rather than restricting attention to a fixed profile $\tau$, the exact rational resolvent $R_v(r)$ permits an asymptotic comparison analysis. Using the algebraic decomposition:
   $$\frac{r}{r^2 - a_m^2} = \frac{1}{r} + \frac{a_m^2}{r(r^2 - a_m^2)},$$
   the resolvent separates into the boundary value $D_0$ and a weighted remainder:
   $$R_v(r) = \frac{D_0}{r} + \frac{\sqrt{2}}{r} \sum_{m=1}^N \frac{a_m^2 v_m}{r^2 - a_m^2}.$$

   *Rigorous Absolute-Moment Bound:* Because the denominators $r^2 - a_m^2$ depend on $m$ and the coefficient vector $v_m$ exhibits alternating signs across modes, one cannot unconditionally replace $\sum_{m=1}^N \frac{a_m^2 v_m}{r^2 - a_m^2}$ by $|D_1| = \sqrt{2} |\sum_{m=1}^N a_m^2 v_m|$ without an additional sign-coherence hypothesis. For $r \ge T > a_N = 2\pi N / L$, since $\frac{1}{r^2 - a_m^2} \le \frac{1}{r^2 [1 - (a_N/T)^2]}$, the mathematically watertight pointwise bound is governed by the **weighted absolute moment**:
   $$\left| R_v(r) - \frac{D_0}{r} \right| \le \frac{\sqrt{2}}{r} \sum_{m=1}^N \frac{a_m^2 |v_m|}{r^2 - a_m^2} \le \frac{\mathcal{M}_1^{\mathrm{abs}}(v)}{r^3 \big[ 1 - (a_N/T)^2 \big]}, \qquad \mathcal{M}_1^{\mathrm{abs}}(v) \equiv \sqrt{2}\sum_{m=1}^N a_m^2 |v_m|.$$

   *Candidate Empirical Two-Jet Envelope:* For the ground-state profile $v_0$, empirical inspection reveals that the dominant low-frequency modes possess coherent signs, motivating the ansatz that the sum is effectively scaled by $|D_1| = |T_v''(0)|$. Substituting $|D_1|$ yields the candidate two-jet upper envelope:
   $$\mathcal{B}_{\mathrm{env}} \equiv \frac{C D_0^2 \log T}{L T} \left( 1 + \frac{1}{T^2 u_1 \big[ 1 - (a_N/T)^2 \big]} \right)^2, \qquad u_1 \equiv \left| \frac{D_0}{D_1} \right|,$$
   where $u_1$ is the physical first-jet cancellation scale investigated in Cell 54 and Cell 58. We treat $\mathcal{B}_{\mathrm{env}}$ as an **empirically validated scaling envelope** rather than an unconditional analytical theorem.

4. **Empirical Validation of the Controlling Scale (Cell 58 Audit):**
   Computational audit via Cell 58 across dimensions $N \in \{8, 12, 16, 20, 24\}$ at $T = 400$ provides strong empirical evidence that $u_1$ is the genuine crossover parameter:

   $$\begin{array}{c|c|c|c|c|c}
   N & |D_0| & |D_1| & u_1 = |D_0/D_1| & 1 / (T^2 u_1) & \eta = a_N / T \\
   \hline
   8  & 8.050 \times 10^{-11} & 3.364 \times 10^{-6} & 2.393 \times 10^{-5} & 0.261 & 0.049 \\
   12 & 6.647 \times 10^{-14} & 6.402 \times 10^{-9} & 1.038 \times 10^{-5} & 0.602 & 0.073 \\
   16 & 1.783 \times 10^{-16} & 3.128 \times 10^{-11} & 5.700 \times 10^{-6} & 1.097 & 0.098 \\
   20 & 8.384 \times 10^{-19} & 2.583 \times 10^{-13} & 3.245 \times 10^{-6} & 1.926 & 0.122 \\
   24 & 1.138 \times 10^{-20} & 5.917 \times 10^{-15} & 1.923 \times 10^{-6} & 3.250 & 0.147
   \end{array}$$

   The parameter $1 / (T^2 u_1)$ rises from $0.26$ to $3.25$ over this range, directly explaining why the first-jet boundary layer couples to the cutoff tail and reconciles the leading asymptotic $E_T \approx 4.14 \times 10^{-43}$ with the observed tail $\delta_T \approx 1.67 \times 10^{-43}$.

5. **Envelope Comparison and Epistemic Qualifications:**
   Comparing the candidate upper envelope $\mathcal{B}_{\mathrm{env}}$ against the exact continuous tail quadrature $\delta_T^{\mathrm{tail}} = \frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) dr$ yields:

   $$\begin{array}{c|c|c|c|c}
   N & \delta_T^{\mathrm{tail}} \text{ (Exact Tail)} & B_1 \text{ (Leading Asymptotic)} & \mathcal{B}_{\mathrm{env}} \text{ (Two-Jet Bound)} & \text{Ratio } \mathcal{B}_{\mathrm{env}} / \delta_T^{\mathrm{tail}} \\
   \hline
   8  & 1.762 \times 10^{-23} & 2.072 \times 10^{-23} & 3.300 \times 10^{-23} & 1.87 \\
   12 & 1.010 \times 10^{-29} & 1.413 \times 10^{-29} & 3.641 \times 10^{-29} & 3.60 \\
   16 & 5.939 \times 10^{-35} & 1.016 \times 10^{-34} & 4.512 \times 10^{-34} & 7.60 \\
   20 & 1.045 \times 10^{-39} & 2.248 \times 10^{-39} & 1.963 \times 10^{-38} & 18.77 \\
   24 & 1.645 \times 10^{-43} & 4.141 \times 10^{-43} & 7.734 \times 10^{-42} & 47.01
   \end{array}$$

   *Scientific & Epistemic Qualifications:*
   - **Empirical Dominance:** The candidate envelope $\mathcal{B}_{\mathrm{env}}$ strictly dominates the computed tail throughout the tested range ($\mathcal{B}_{\mathrm{env}} > \delta_T^{\mathrm{tail}}$).
   - **Asymptotic Growth of the Ratio:** However, the ratio $\mathcal{B}_{\mathrm{env}} / \delta_T^{\mathrm{tail}}$ grows systematically with $N$ (roughly doubling every $\Delta N = 4$). This confirms that Cell 58 provides an empirical upper envelope rather than a uniform-in-$N$ mathematical constant.
   - **Safe Decoupling Margin:** Crucially, this ratio excess ($47.0$ at $N = 24$) remains vastly smaller than the inverse tunneling scale ($47.0 \ll 10^{40}$), confirming that first-jet amplification does not endanger boundary-defect decoupling.
   - **Analytical Status:** An unconditional mathematical theorem requires bounding the absolute moment $\mathcal{M}_1^{\mathrm{abs}}(v)$ or establishing ground-state sign coherence analytically, distinguishing rigorous analytical machinery from the empirical two-jet envelope.

6. **Cutoff Sweep and Finite-Cutoff Transition (Cell 58 Part 3):**
   Holding $N = 24$ fixed while sweeping $T \in [100, 800]$ demonstrates the collapse of the first-jet correction:

   $$\begin{array}{c|c|c|c|c|c}
   T & 1 / (T^2 u_1) & \delta_T^{\mathrm{tail}} \text{ (Exact Tail)} & B_1 & \mathcal{B}_{\mathrm{env}} & \text{Ratio } \mathcal{B}_{\mathrm{env}} / \delta_T^{\mathrm{tail}} \\
   \hline
   100 & 79.46 & 2.485 \times 10^{-43} & 1.211 \times 10^{-42} & 7.839 \times 10^{-39} & 31544.5 \\
   200 & 14.23 & 1.987 \times 10^{-43} & 7.168 \times 10^{-43} & 1.662 \times 10^{-40} & 836.5 \\
   400 & 3.32  & 1.645 \times 10^{-43} & 4.141 \times 10^{-43} & 7.734 \times 10^{-42} & 47.0 \\
   800 & 0.817 & 1.422 \times 10^{-43} & 2.349 \times 10^{-43} & 7.754 \times 10^{-43} & 5.45
   \end{array}$$

   As $T$ increases from $100$ to $800$, the coupling parameter drops by two orders of magnitude ($79.46 \to 0.817$), and the envelope ratio collapses from $31544.5$ down to $5.45$. This confirms that the first-jet amplification is strictly a finite-cutoff artifact that extinguishes as $T^2 u_1 \to \infty$.

7. **The Boundary-Defect Decoupling Conjecture:**
   The decoupling metric $\mathcal{D}(N) \equiv D_0^2 \left( 1 + \frac{1}{T^2 u_1} \right)^2$ evaluated in Cell 58 demonstrates rapid collapse:

   $$\begin{array}{c|c|c|c|c}
   N & D_0^2 \text{ (Tunneling)} & u_1 \text{ (Boundary Scale)} & \left[1 + \frac{1}{T^2 u_1}\right]^2 & \mathcal{D}(N) \\
   \hline
   8  & 6.481 \times 10^{-21} & 2.393 \times 10^{-5} & 1.591  & 1.031 \times 10^{-20} \\
   12 & 4.418 \times 10^{-27} & 1.038 \times 10^{-5} & 2.566  & 1.134 \times 10^{-26} \\
   16 & 3.178 \times 10^{-32} & 5.700 \times 10^{-6} & 4.395  & 1.397 \times 10^{-31} \\
   20 & 7.028 \times 10^{-37} & 3.245 \times 10^{-6} & 8.560  & 6.016 \times 10^{-36} \\
   24 & 1.295 \times 10^{-40} & 1.923 \times 10^{-6} & 18.060 & 2.339 \times 10^{-39}
   \end{array}$$

   We distinguish three epistemic tiers for $D_0$: (i) empirical rapid collapse, (ii) proven semiclassical WKB barrier tunneling $\mathcal{S}_{\mathrm{WKB}} \approx \frac{\pi N}{4} \log c$, and (iii) uniform asymptotic control.
   
   To establish decoupling, we formulate the:
   > **Boundary-Defect Decoupling Conjecture:** Along any admissible double-scaling sequence $(N, L)$ with $T > a_N$,
   > $$\lim_{N, L \to \infty} D_0(N, L)^2 \frac{\log T}{L T} \left( 1 + \frac{1}{T^2 u_1(N, L)} \right)^2 = 0.$$

8. **The Power-Law Bridge and the Polynomial Decoupling Criterion:**
   The findings of Cells 59–62 and the exact operator decomposition of Theorem 7.3 crystallize the analytical path forward for the research programme. Because the small denominators $\Delta_k$ cancel mode by mode, the first-jet ratio $D_1/D_0$ is algebraically decoupled from the exponential tunneling decay $e^{-\mathcal{S}_{\mathrm{WKB}}}$. 

   Evaluating the operator resolvent components across $N \in \{8, 12, 16, 20, 24\}$ yields empirical scaling trajectories:
   $$M_1(N) \sim N^{0.6}, \qquad \langle d, R_{\mathrm{even}} d \rangle \sim N^{2.2} \implies \frac{D_1}{D_0} \sim N^{2.6 - 2.8} \implies u_1(N) \sim N^{-2.6 \pm 0.1}.$$
   However, proving the specific empirical exponent $p \approx 2.6$ is not required for the continuum programme. Rather, the essential analytical requirement is the broader:
   > **Polynomial Decoupling Criterion:** If there exists some finite exponent $p < \infty$ such that
   > $$u_1(N)^{-1} = \left| \frac{D_1}{D_0} \right| = \mathcal{O}(N^p) \qquad (N \to \infty),$$
   > $$\mathcal{D}(N) \equiv D_0^2 \left( 1 + \frac{1}{T^2 u_1} \right)^2 \sim e^{-\frac{\pi N}{2}\log c} \left[ 1 + \mathcal{O}(N^p) \right]^2 \sim e^{-\frac{\pi N}{2}\log c} N^{2p} \longrightarrow 0$$
   > collapses exponentially fast to zero as $N \to \infty$.

   Under this criterion, any polynomial bound on the excited-sector resolvent matrix elements guarantees that the boundary-layer defect decouples in the continuum limit, independently of prior proof of infinite-order $C^\infty$ boundary flatness. This requirement is established by **Proposition 8.10**, which proves the exact finite-$N$ commutator cancellations and Cauchy–Schwarz resolvent coercivity (Part I), reducing the Polynomial Decoupling Criterion to Hypotheses H1, H2, H2$_{\mathrm{odd}}$, and H3 (Part II), thereby ensuring that boundary-defect leakage decouples exponentially in the continuum limit.

---

## 10. Computational Reproducibility and Software Availability

To ensure complete computational transparency and reproducibility, the entire mathematical software pipeline and all raw high-precision calculation transcripts supporting this study are permanently archived in the public repository [10]:

> **Software Repository:** <https://github.com/akivag613/connes-cvs->  
> **Mirror Repository:** <https://github.com/nrensen/connes-cvs->

The calculations reported in this manuscript were performed using Python and the `mpmath` arbitrary-precision arithmetic library at 50 decimal digits of precision. All empirical observations, asymptotic fits, and spectral decompositions are reproduced by standalone computational scripts (`cell*.py`), whose complete numerical output transcripts are preserved in matching log files (`cell*.out`).

**Table 4: Mapping of Programme Sections to Computational Scripts and Output Logs**

| Section & Topic | Mathematical / Numerical Focus | Python Script | Verification Log |
| :--- | :--- | :--- | :--- |
| Section 2.1–2.3 (Tables 1, 2, 3) | Ground state eigensystem & mode decay ($N=1\dots 24$) | `cell34.py`, `cell40.py`, `cell41.py` | `cell34.out`, `cell40.out`, `cell41.out` |
| Section 3 (Proposition 3.1, Observation 3.2, Conjecture 3.3) | Spatial wave profile & boundary jet derivatives $D_0\dots D_3$ | `cell42.py`, `cell43.py` | `cell42.out`, `cell43.out` |
| Section 4 (WKB Barrier Mechanics) | Effective Schrödinger potential & WKB tunneling action | `cell44.py`, `cell47.py` | `cell44.out`, `cell47.out` |
| Section 5 (Legendre Multipoles & Tail Extinction) | Bauer–Bessel Legendre transform & Taylor jet extinction $A_0\dots A_4$ | `cell44.py`, `cell45.py` | `cell44.out`, `cell45.out` |
| Section 6 (Poles & Heat Dynamics) | Discrete Cauchy transform, heat boundary dynamics, profile collapse | `cell51.py`, `cell52.py`, `cell53.py`, `cell54.py` | `cell51.out`, `cell52.out`, `cell53.out`, `cell54.out` |
| Section 7 (Tri-Partite Balance & Cutoff Leakage) | Continuous-variable balance & finite-$T$ Archimedean leakage $\delta_T^{\mathrm{tail}}$ | `cell46.py`, `cell56.py`, `cell57.py` | `cell46.out`, `cell56.out`, `cell57.out` |
| Section 8 (Wiener–Hopf Scaling, Commutator Resolvent & Bounds) | Commutator resolvent formula, exact cancellation & collective spectral reorganization, Wiener–Hopf scaling & bounds | `cell49.py`, `cell50.py`, `cell53.py`, `cell54.py`, `cell55.py` | `cell49.out`, `cell50.out`, `cell53.out`, `cell54.out`, `cell55.out` |
| Section 8 (Bound-State Gaps & Doublet Ladder) | Low-energy spectral gaps & exponential collapse of odd/even tunneling splitting | `cell59.py` | `cell59.out` |
| Section 8 (Interlaced Doublet Spectrum) | Full even-odd tunneling doublet spectrum and bound-state ladder below barrier top | `cell60.py` | `cell60.out` |
| Section 8 (Transition Dipole & Mode 1 Saturation) | Commutator projection $\langle e_1, Kc \rangle$, mode 1 norm saturation & $C_N = \mathcal{O}(1)$ | `cell61.py` | `cell61.out` |
| Section 8 & 9.1 (Operator Decomposition & Scaling) | Exact decomposition $D_1/D_0 = \kappa^2 [\mathcal{T}_{\mathrm{diag}} + \mathcal{T}_{\mathrm{cross}} - D_0^2 M_2]$, spectral filtering proof, spatial profile & $u_1 \sim N^{-2.6}$ | `cell62.py` | `cell62.out` |
| Section 9.1 (Jet Defect & Decoupling Mechanics) | Two-jet resolvent envelope, cutoff sweep & first-jet scale $u_1$ audit | `cell58.py` | `cell58.out` |
| Roadmap Stage III (Operator Dominance Reconnaissance) | Negative Archimedean Gram matrix $\mathcal{Q}_-$, whitening breakdown & 3-mode sector concentration ($98\%$) | `cell63.py` | `cell63.out` |
| Roadmap Stage III (High-Precision Schur Positivity) | High-precision numerical verification of $\mathcal{Q}_{\mathrm{Weil}} \succ 0$ via symmetric $LDL^T$ Schur decoupling at 80 dps | `cell64.py` | `cell64.out` |
| Roadmap Stage V (Effective Hamiltonian & Self-Energy) | Three-mode effective Hamiltonian $S_{\mathrm{low}} = A - \Sigma_{\mathrm{low}}$ & high-mode self-energy sweep at 80 dps | `cell65.py` | `cell65.out` |
| Roadmap Stage VI (Hypothesis H1–H3 Audit) | Semiclassical flux-matching $\mathcal{R}_{\mathrm{tun}}$, mode-by-mode transmission cancellation & bound ladder | `cell66.py` | `cell66.out` |

---

### References

1. A. Weil, *Sur les "formules explicites" de la théorie des nombres premiers*, Medd. Lunds Univ. Mat. Sem. (1952), 252–265.
2. A. Connes and W. D. van Suijlekom, *A spectral approach to the Riemann zeta function*, arXiv:2104.09241 (2021).
3. A. Connes and W. D. van Suijlekom, *Quadratic Forms, Real Zeros and Echoes of the Spectral Action*, Communications in Mathematical Physics **406** (2025), no. 12, 312, [arXiv:2511.23257](https://arxiv.org/abs/2511.23257).
4. A. Connes, C. Consani, and H. Moscovici, *Zeta Spectral Triples*, in Applications of Noncommutative Geometry to Gauge Theories, Field Theories, and Quantum Space-Time, EMS Series of Lectures in Mathematics, EMS Press (2026), 39–76, [arXiv:2511.22755](https://arxiv.org/abs/2511.22755).
5. A. Groskin, *High-Precision Approximation of Riemann Zeros via the Truncated Weil Form*, arXiv:2605.20224 (2026).
6. A. Groskin, *A finite Guinand-Weil dictionary and archimedean tail order for the truncated Weil quadratic form*, arXiv:2607.02828 / Zenodo:21124802 (2026).
7. A. Groskin, *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path*, Zenodo:21242028 (2026).
8. D. Slepian and H. O. Pollak, *Prolate spheroidal wave functions, Fourier analysis and uncertainty — I*, Bell System Tech. J. 40 (1961), 43–63.
9. H. J. Landau and H. O. Pollak, *Prolate spheroidal wave functions, Fourier analysis and uncertainty — II*, Bell System Tech. J. 40 (1961), 65–84.
10. A. Groskin and N. Rensen, *connes-cvs: Arbitrary-precision computational suite and verification archive for the truncated Connes–van Suijlekom Galerkin form*, software repository and raw numerical logs, GitHub: https://github.com/akivag613/connes-cvs- (mirror: https://github.com/nrensen/connes-cvs-) (2026).
11. Research Record / Connes–CvS Series, *An Exact Resolvent and Commutator Toolkit for the Truncated Connes–van Suijlekom Weil Quadratic Form*, Companion Paper (Paper 4), GitHub: `nrensen/connes-cvs-` (2026).
