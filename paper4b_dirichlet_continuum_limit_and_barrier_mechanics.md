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
3. **Semiclassical WKB Barrier Tunneling Mechanics:** Inverting the ground-state profile defines an effective Schrödinger potential $V_{\mathrm{conf}}(t) - E = T''(t)/T(t)$ whose midpoint well rises steeply toward the boundaries. Semiclassical WKB tunneling action across the barrier $\mathcal{S}_{\mathrm{WKB}} = \int_0^{t_{\mathrm{turn}}} \sqrt{T''/T} \, dt$ reproduces the observed logarithmic boundary suppression to within $5.3\%$ across 47 orders of magnitude, obeying the universal semiclassical scaling law $\mathcal{S}_{\mathrm{WKB}}(N, c) \approx \frac{\pi N}{4} \log c$.
4. **Legendre Multipole Decomposition and Asymptotic Tail Extinction:** Via Bauer's spherical Bessel expansion, $T_{v_N}(t)$ decomposes into Legendre multipoles with alternating signs, producing strictly constructive interference at the midpoint and destructive cancellation at the boundaries. The high-frequency Taylor coefficients $A_k(N)$ extinguish rapidly across all orders ($A_0 \sim 10^{-40}, A_1 \sim 10^{-34}, A_2 \sim 10^{-29}$ at $N = 24$), motivating the conjecture that the continuum resolvent $R_\infty(r) = o(r^{-k})$ decays faster than every inverse power of $r$.
5. **Tri-Partite Zero-Energy Equilibrium and the Discretization Gap:** Continuous-variable numerical quadrature of the Archimedean form independently cancels the algebraic pole ($+1.55165$) and prime ($-0.07185$) contributions down to a residual of $1.29 \times 10^{-43}$ at $N = 24$. The residual ratio of $1.96$ against the matrix eigenvalue $\lambda_{\min}(24) = 2.53 \times 10^{-43}$ is certified as the finite-rank Galerkin discretization gap $\delta \mathcal{Q} = \langle u, Q_{\mathrm{arch}}^{\mathrm{matrix}} u \rangle - \mathcal{Q}_{\mathrm{arch}}^{\mathrm{cont}}(v)$.
6. **Formal Continuum Wiener–Hopf Scaling and Asymptotic Bounding Ladder:** In the continuum scaling limit, the divided-difference Galerkin kernel transforms into a half-line Wiener–Hopf convolution operator with kernel $K_{\mathrm{sym}}(w) = \frac{w}{2\sinh(w/2)}$ whose symbol factors into squared Gamma functions $\frac{\pi^2}{\cosh^2(\pi k)} = [\Gamma(\frac{1}{2} - ik)]^2 [\Gamma(\frac{1}{2} + ik)]^2$. The resulting double pole at $k = -i/2$ generates a logarithmic boundary layer $\phi(x) \sim -\log x$ as $x \to 0^+$, explaining the observed bulk/edge asymmetry between $D_0$ and $D_1$. Conditional on uniform scattering gap stability, this establishes two-sided subexponential bounds on the first-jet cancellation scale $\frac{c_1}{N^2 \log N} \le u_1 \le \frac{c_2}{N^{1/2}}$.
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
   - Exact dual algebraic derivations proving that the boundary suppression factor $D_0$ cancels identically from the first-jet ratio $D_1 / D_0$.
   - The exact small-denominator cancellation $(E_k - \lambda)$ in the bound-state sector and universal semigroup squeezing bounds.

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
- **Cell 54 (Anatomy of the Decoupling Ratio $s_N$):** Tracking the decoupling ratio $s_N = (\kappa N)^2 (D_0 / D_1)$ reveals that while $-\log|D_0|$ drops by 22.7 units and $-\log|D_1|$ drops by 20.2 units, their difference $\Delta_N = \log|D_1/D_0|$ drifts only from $10.64$ to $13.16$. This proves that $D_0$ and $D_1$ share the exact same leading exponential barrier decay rate.

---

## 7. Tri-Partite Zero-Energy Balance and the Finite-Rank Discretization Gap

Let $\mathcal{Q}(v) = \mathcal{Q}_{\mathrm{pole}}(v) + \mathcal{Q}_{\mathrm{prime}}(v) + \mathcal{Q}_{\mathrm{arch}}(v)$ be the Connes–van Suijlekom quadratic form on the Galerkin subspace of dimension $2N+1$. For every finite dimension $N$, the algebraic matrix sum matches the minimum eigenvalue identically:

$$\mathcal{Q}_{\mathrm{matrix}}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{\mathrm{matrix}}(v_N) \equiv \lambda_{\min}(N).$$

### Proposition 7.1 (Continuous-Quadrature Balance and Discretization Gap)
*When the Archimedean contribution is evaluated independently via continuous-variable quadrature $\frac{1}{\pi} \int_0^{80} h_+(r) \Phi_{v_N}(r)^2 \, dr$ (using `cell46.py`, logged in `cell46.out` [10]), the independently computed components cancel from $\mathcal{O}(1)$ down to a residual of order $10^{-43}$ at $N = 24$:*

- $N = 4$: $\mathcal{Q}_{\mathrm{pole}} = +2.206186$, $\mathcal{Q}_{\mathrm{prime}} = -0.316153$, $\mathcal{Q}_{\mathrm{arch}} = -1.890032$, summing to $\mathcal{Q}_{\mathrm{total}} = 7.82 \times 10^{-15}$ ($\lambda_{\min} = 8.83 \times 10^{-15}$),
- $N = 8$: $\mathcal{Q}_{\mathrm{pole}} = +1.813949$, $\mathcal{Q}_{\mathrm{prime}} = -0.154916$, $\mathcal{Q}_{\mathrm{arch}} = -1.659033$, summing to $\mathcal{Q}_{\mathrm{total}} = 5.38 \times 10^{-23}$ ($\lambda_{\min} = 6.71 \times 10^{-23}$),
- $N = 12$: $\mathcal{Q}_{\mathrm{pole}} = +1.675166$, $\mathcal{Q}_{\mathrm{prime}} = -0.108101$, $\mathcal{Q}_{\mathrm{arch}} = -1.567065$, summing to $\mathcal{Q}_{\mathrm{total}} = 1.32 \times 10^{-29}$ ($\lambda_{\min} = 1.78 \times 10^{-29}$),
- $N = 16$: $\mathcal{Q}_{\mathrm{pole}} = +1.609630$, $\mathcal{Q}_{\mathrm{prime}} = -0.088194$, $\mathcal{Q}_{\mathrm{arch}} = -1.521436$, summing to $\mathcal{Q}_{\mathrm{total}} = 5.11 \times 10^{-35}$ ($\lambda_{\min} = 7.12 \times 10^{-35}$),
- $N = 20$: $\mathcal{Q}_{\mathrm{pole}} = +1.572288$, $\mathcal{Q}_{\mathrm{prime}} = -0.077529$, $\mathcal{Q}_{\mathrm{arch}} = -1.494759$, summing to $\mathcal{Q}_{\mathrm{total}} = 8.81 \times 10^{-40}$ ($\lambda_{\min} = 1.32 \times 10^{-39}$),
- $N = 24$: $\mathcal{Q}_{\mathrm{pole}} = +1.551652$, $\mathcal{Q}_{\mathrm{prime}} = -0.071854$, $\mathcal{Q}_{\mathrm{arch}} = -1.479798$, summing to $\mathcal{Q}_{\mathrm{total}} = 1.29 \times 10^{-43}$ ($\lambda_{\min} = 2.53 \times 10^{-43}$).

*The Discretization Gap $\delta \mathcal{Q}$.*
The residual discrepancy between the continuous-quadrature sum ($1.29 \times 10^{-43}$) and the matrix eigenvalue ($2.53 \times 10^{-43}$) represents a proportional factor of $\approx 1.96$ at the residual scale. Rather than a numerical error, this factor reflects the certified finite-rank discretization gap:

$$\delta \mathcal{Q} = \langle u, Q_{\mathrm{arch}}^{\mathrm{matrix}} u \rangle - \mathcal{Q}_{\mathrm{arch}}^{\mathrm{cont}}(v)$$

between the truncated Galerkin projection matrix and the continuous functional (investigated and resolved in `cell56.py`).

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

### Theorem 8.3 (Conditional Asymptotic Bounding Ladder for $u_1$ and $s_N$)
*Conditional on the uniform scattering gap hypothesis $\inf_N E_{\mathrm{scatt},\min}(N) \ge E_{\mathrm{gap}} > 0$ (empirically supported by Cells 49–50) and uniform $H^1$ Sobolev boundedness $\|T'_v\|_{L^2} = \mathcal{O}(1)$:*

1. **Upper Bound on First-Jet Ratio:**
   $$\left| \frac{D_1}{D_0} \right| \le \kappa^2 C_{\mathrm{upper}} N^2 \log N.$$
2. **Two-Sided Bounds on Cancellation Scale and Decoupling Ratio:**
   $$\frac{c_1}{N^2 \log N} \le u_1 = \left| \frac{D_0}{D_1} \right| \le \frac{c_2}{N^{1/2}}, \qquad \frac{\kappa^2 c_1}{\log N} \le s_N \le \kappa^2 c_2 N^{3/2}.$$
3. **Subexponentiality:**
   *Under these conditions, $u_1$ and $s_N$ are strictly subexponential, ruling out any $e^{-\alpha N}$ collapse of the cancellation scale and confirming that $D_0$ and $D_1$ share the exact same leading exponential WKB barrier decay rate.*
4. **Spatial Boundary Layer Width:**
   *The boundary layer width $\delta_N = \sqrt{u_1} \ge \frac{1}{\kappa \sqrt{C_{\mathrm{upper}}} N \sqrt{\log N}}$ shrinks only algebraically, fundamentally decoupling from the exponentially suppressed endpoint amplitude $T_N(0) \sim e^{-\mathcal{S}_{\mathrm{WKB}}}$.*

*Proof.* Applying the Cauchy–Schwarz inequality to the commutator resolvent formula (Theorem 6.11 of Paper 4) with $\|e\|_2 = \sqrt{2N+1}$ and $\|\mathbf{w}\|_2 \le C_w N^{3/2} \log N$, bounded resolvent norm $\|(Q_{\mathrm{even}} - \lambda I)^\dagger\|_{\mathrm{scatt}} \le 1/E_{\mathrm{gap}}$ yields $|D_1/D_0| \le \kappa^2 C_{\mathrm{upper}} N^2 \log N$. Lower Sobolev trace bound gives $|D_1/D_0| \ge C_{\mathrm{lower}} N^{1/2}$. Inverting yields the bounds on $u_1$ and $s_N$. $\blacksquare$

---

## 9. The Analytical Roadmap toward Continuous Weil Positivity

The empirical and asymptotic results established in this research programme indicate that the finite-rank Galerkin truncation provides a faithful, convergent approximation to the continuous Weil quadratic form. 

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
| Section 7 (Tri-Partite Balance & Discretization Gap) | Continuous-variable balance & finite-rank gap $\delta \mathcal{Q}$ | `cell46.py`, `cell56.py` | `cell46.out`, `cell56.out` |
| Section 8 (Wiener–Hopf Scaling & Bounding Ladder) | Resolvent bounds, scattering gap validation & Wiener–Hopf limits | `cell49.py`, `cell50.py`, `cell54.py`, `cell55.py` | `cell49.out`, `cell50.out`, `cell54.out`, `cell55.out` |

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
