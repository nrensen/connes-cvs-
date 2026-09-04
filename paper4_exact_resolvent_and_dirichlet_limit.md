# Exact Rational Resolvent and Pointwise Positivity of the Archimedean Kernel in the Truncated Weil Quadratic Form: Numerical Evidence for a Dirichlet Continuum Limit

**Authors:** Research Record / Connes–CvS Investigation Series  
**Date:** September 2026  
**Software & Reproducibility Suite:** `https://github.com/akivag613/connes-cvs-` (mirror: `nrensen/connes-cvs-`)  
**Status:** Standalone Manuscript

---

### Abstract

The truncated Weil quadratic form developed by Connes–van Suijlekom and Connes–Consani–Moscovici at prime cutoff $c > 1$ and band $N$ produces finite-rank Galerkin matrices whose deep spectra provide an explicit computational window into Weil positivity and the Riemann Hypothesis. The omitted Archimedean tail of this truncation has historically been treated as a difficult oscillatory numerical integration problem or as an empirical asymptotic inverse-power expansion.

In this paper, we establish the exact algebraic solution to the finite-$N$ Archimedean kernel and explore its infinite-dimensional limit $N \to \infty$ through exact theorems, empirical observations, and precise conjectures:

1. **Exact Rational Resolvent & Operator Identity (Theorem):** Starting from the four-term analytic reduction of the Archimedean Volterra integral, we prove algebraically and independently of numerical quadrature that the reduced Fourier kernel $R_v(r) = K_{\mathrm{Fourier}}(v, r, L) / (1 - \cos(rL))$ is identically equal to the squared Cauchy resolvent:
   $$R_v(r) \equiv \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^{N} \frac{r v_m}{r^2 - a_m^2} \right]^2, \qquad a_m = \frac{2\pi m}{L},$$
   on the punctured complex plane $\mathbb{C} \setminus \{0, \pm a_1, \dots, \pm a_N\}$. Furthermore, introducing the Neumann Laplacian $\mathcal{L} = -d^2/dt^2$ on $[0, L]$ with $T'(0) = T'(L) = 0$, the rational generating function $D(z)$ is the boundary evaluation of the operator resolvent:
   $$D(z) \equiv \big[(I + z\mathcal{L})^{-1} T_v\big](0) = \int_0^\infty e^{-s} \big[ e^{-sz\mathcal{L}} T_v \big](0) \, ds \quad (\operatorname{Re} z > 0),$$
   whose Taylor expansion around $z = 0$ reproduces the entire endpoint-jet hierarchy $D_k = T_v^{(2k)}(0)$.
2. **Universal Fourier Factorization and Unconditional Kernel Positivity (Theorem):** The entire Fourier-side amplitude $\Phi_v(r)$ factors directly in terms of the boundary Neumann resolvent evaluated at the inverted spectral variable $z = -1/r^2$:
   $$\Phi_v(r) \equiv \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D\left(-\frac{1}{r^2}\right), \qquad K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 = \frac{4}{L} \frac{\sin^2(rL/2)}{r^2} D\left(-\frac{1}{r^2}\right)^2 \ge 0,$$
   proving algebraically and unconditionally that the Fourier-side Archimedean kernel $K_{\mathrm{Fourier}}(v, r, L)$ is pointwise non-negative for all real $r$ and all real coefficient vectors $v \in \mathbb{R}^{N+1}$. (This pointwise non-negativity does not imply positivity of the integrated Archimedean quadratic form $Q_{\mathrm{arch}}$, whose digamma weight $h_+(r)$ changes sign).
3. **Spectral Lattice Orthogonality (Theorem):** At the lattice nodes $r = a_m$, the apparent poles cancel cleanly against the envelope zeros via removable singularities, yielding the exact sampling identity:
   $$K_{\mathrm{Fourier}}(v, 0, L) = L v_0^2 = L u_0^2, \qquad K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 = L u_m^2 \quad (m = 1, \dots, N),$$
   uncoupling the kernel into the squared Fourier coefficients.
4. **Observed Mode Concentration and Asymptotic Laws (Numerical & Conjectural):** Across 24 Galerkin dimensions ($N = 1, \dots, 24$), the computed ground states exhibit rapidly decreasing successive differences and strong concentration of their $\ell^2$ mass in the lowest modes (over $99.98\%$ in $m \le 4$). The boundary value drops from $7.52 \times 10^{-3}$ to $1.14 \times 10^{-20}$ (approximately 17.8 decimal orders of magnitude), with the effective decay exponent decreasing toward values near $1.1$:
   $$|T_{v_N}(0)| \sim C(c) \cdot \rho(c)^N \quad \text{with } 0 < \rho(c) < 1 \quad \text{(Conjectured)}.$$
   While both quantities span approximately 43 decimal orders of magnitude, the ground-state eigenvalue $\lambda_{\min}(N)$ appears asymptotically proportional to the boundary leakage energy:
   $$\lambda_{\min}(N) \sim \kappa_c \cdot [T_{v_N}(0)]^2 \longrightarrow 0 \quad \text{(Numerical Conjecture)}.$$
5. **Continuum Solitary Wave and Dirichlet Nodes (Conjectural):** Numerical evidence suggests that in the continuum limit $N \to \infty$, $T_{v_N}(t)$ appears to converge to a symmetric, strictly positive solitary wave $T_\infty(L - t) = T_\infty(t)$ with dual Dirichlet boundary vanishing $T_\infty(0) = T_\infty(L) = 0$ and conjectured infinite-order flat boundary contact $\widetilde{T}_\infty \in C_c^\infty(\mathbb{R})$ with $\operatorname{supp} \widetilde{T}_\infty \subseteq [0, L]$.
   Conditional on this flat contact, the Volterra boundary jump at $\omega = 1$ is eliminated, removing the finite-rank obstruction to Weil positivity.

---

## 1. Introduction

The explicit formula of Guinand and Weil relates the nontrivial zeros of the Riemann zeta function $\zeta(s)$ to arithmetic prime-power sums, pole contributions, and Archimedean gamma-factor terms. In André Weil's formulation (1952), the Riemann Hypothesis (RH) is equivalent to the non-negativity of the associated quadratic functional:

$$W(g) \ge 0$$

on all admissible test functions $g = f * f^*$.

In Alain Connes’ non-commutative geometry program, Weil positivity is realized through an operator-theoretic spectral framework on the prolate spheroidal wave spaces of band-limited functions. Recent work by Connes and van Suijlekom (2025) and Connes, Consani, and Moscovici (2026) discretizes the continuous Weil form using a finite-rank Galerkin projection: for a logarithmic prime cutoff $L = \log c$ (with $c > 1$) and a finite frequency band $N \ge 1$, the continuous form is projected onto an explicit $(2N+1) \times (2N+1)$ matrix $Q_{c, N}$.

The total quadratic form decomposes into prime, pole, and Archimedean components:

$$\langle v, Q v \rangle = \langle v, Q_{\mathrm{prime}} v \rangle + \langle v, Q_{\mathrm{pole}} v \rangle + \langle v, Q_{\mathrm{arch}} v \rangle.$$

While the prime and pole contributions admit exact closed finite algebraic representations, the Archimedean term involves an integral over the positive real axis:

$$\langle v, Q_{\mathrm{arch}} v \rangle = \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr,$$

where

$$h_+(r) = \operatorname{Re} \psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi$$

is the smooth Archimedean density.

### The Historical Problem of the Archimedean Tail

In finite implementations, numerical quadrature must truncate the $r$-integral at some cutoff $T$. This truncation introduced significant empirical and analytical difficulties:
1. **Numerical Quadrature Failure:** At large $T$, oscillatory quadrature routines (such as unsubdivided Gauss–Legendre or double-exponential rules) suffer severe cancellation errors and false convergence.
2. **Heuristic Asymptotics:** In earlier computational iterations of this investigation series, the Archimedean tail was approached via oscillatory numerical quadrature or by expanding $K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r)$ as an asymptotic inverse-power series $\sum_{k \ge 0} A_k / r^{2k+2}$. However, the coefficients $A_k$ appeared as highly non-trivial combinatorial sums of spectral moments, and bounding the remainder $\varepsilon_N(r)$ remained an open obstacle.
3. **Question of Positivity:** The positivity of $K_{\mathrm{Fourier}}(v, r, L)$ was not made explicit in earlier finite-$N$ formulations, leaving open whether sign-oscillations could induce negative eigenvalues at large $T$.

This paper establishes the exact closed-form algebraic solution to this problem, proves global finite-$N$ non-negativity of the Fourier-side Archimedean kernel independently of numerical quadrature, and formulates the precise conjectures governing the infinite-dimensional limit $N \to \infty$.

---

## 2. Geometric Setup and the Volterra Kernel

Let $c > 1$ and define the logarithmic interval length $L = \log c$. Let $v = (v_0, v_1, \dots, v_N)^\top \in \mathbb{R}^{N+1}$ be a canonical real-even coefficient vector, normalized such that $\|v\|_2^2 = \sum_{m=0}^N v_m^2 = 1$.

The canonical vector $v$ maps to full symmetric Fourier coefficients $u = (u_{-N}, \dots, u_N)^\top \in \mathbb{R}^{2N+1}$ via:

$$u_0 = v_0, \qquad u_{+m} = u_{-m} = \frac{v_m}{\sqrt{2}} \quad (m = 1, \dots, N).$$

### 2.1 The Trigonometric Wave and the Spatial Volterra Kernel

The vector $v$ generates an even trigonometric polynomial on the physical interval $[0, L]$:

$$T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^{N} v_m \cos\left(\frac{2\pi m t}{L}\right) \qquad (0 \le t \le L),$$

which in normalized coordinates $s = t/L \in [0, 1]$ is represented by:

$$\tau_v(s) := T_v(L s) = v_0 + \sqrt{2} \sum_{m=1}^{N} v_m \cos(2\pi m s) \qquad (0 \le s \le 1).$$

The quadratic spatial kernel entering the Archimedean explicit formula is the normalized Volterra auto-convolution for $\omega \in [0, 1]$:

$$K_v(\omega) = 2 \int_0^\omega \tau_v(s) \tau_v(\omega - s) \, ds = \frac{2}{L} \int_0^{L\omega} T_v(t) T_v(L\omega - t) \, dt \qquad (0 \le \omega \le 1).$$

In physical coordinates $x = L\omega \in [0, L]$, the convolution evaluates to $K_v^{\mathrm{phys}}(x) = 2 \int_0^x T_v(t) T_v(x - t) \, dt = L K_v(x/L)$.

### 2.2 The Fourier-Side Representation

Transforming $K_v(\omega)$ to the spectral variable $r \in \mathbb{R}$ against $\cos(r L \omega)$ on $[0, 1]$ produces the Fourier-side kernel $K_{\mathrm{Fourier}}(v, r, L)$. 

By direct integration of the Volterra convolution against $\cos(r L \omega)$, the boundary terms at the endpoint $\omega = 1$ ($x = L$) factor out cleanly (implemented and verified in the companion calculation script `cell32.py` and output log `cell32.out` [10]), isolating the common oscillatory factor $1 - \cos(rL)$:

$$K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r),$$

where $R_v(r)$ is a rational function of $r^2$. Defining the fundamental lattice scale:

$$\kappa = \frac{2\pi}{L}, \qquad a_m = \kappa m = \frac{2\pi m}{L},$$

and the elementary mode integrals:

$$S_{\mathrm{bar}}(m) = \frac{a_m}{a_m^2 - r^2}, \qquad C_{\mathrm{bar}}(m) = \frac{r^2 + a_m^2}{L (r^2 - a_m^2)^2},$$

the rational kernel was originally defined by the four-part sum:

$$\begin{aligned}
R_v(r) &= \frac{2 v_0^2}{L r^2} + 2 \sum_{m=1}^{N} v_m^2 C_{\mathrm{bar}}(m) - \frac{2\sqrt{2} v_0}{\pi} \sum_{m=1}^{N} \frac{v_m S_{\mathrm{bar}}(m)}{m} - \frac{1}{\pi} \sum_{m=1}^{N} \frac{v_m^2 S_{\mathrm{bar}}(m)}{m} \\
&\quad + \frac{4}{\pi} \sum_{1 \le m < n \le N} v_m v_n \frac{m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n)}{n^2 - m^2}.
\end{aligned}$$

---

## 3. The Exact Rational Resolvent Identity

Using the analytic expression for the reduced kernel $R_v(r)$ derived from the Archimedean Volterra integral (implemented and cross-validated in the companion calculation script `cell32.py` [10]), we now state the first central theorem of this paper. This result proves **algebraically, in closed form, and completely independently of numerical quadrature**, that the four-term interaction sum matches the square of a single Cauchy resolvent at every finite dimension $N$.

### Theorem 3.1 (Exact Rational Resolvent Identity)
*For any canonical coefficient vector $v \in \mathbb{R}^{N+1}$ and any $r \in \mathbb{C} \setminus \{0, \pm a_1, \dots, \pm a_N\}$, the reduced rational kernel $R_v(r)$ is identically equal to the square of a single rational sum:*

$$R_v(r) \equiv \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^{N} \frac{r v_m}{r^2 - a_m^2} \right]^2.$$

*Furthermore, defining the boundary-resolvent rational function:*

$$D(z) := v_0 + \sqrt{2} \sum_{m=1}^{N} \frac{v_m}{1 + a_m^2 z},$$

*and the coefficient-generating rational function:*

$$A(z) := \frac{2}{L} D(-z)^2 = \frac{2}{L} \left[ v_0 + \sqrt{2} \sum_{m=1}^{N} \frac{v_m}{1 - a_m^2 z} \right]^2,$$

*the rational kernel satisfies the identity:*

$$R_v(r) \equiv \frac{1}{r^2} A\left(\frac{1}{r^2}\right) \qquad \forall r \in \mathbb{C} \setminus \{0, \pm a_1, \dots, \pm a_N\}.$$

### Proof of Theorem 3.1
We expand the square in the proposed formula:

$$\frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^{N} \frac{r v_m}{r^2 - a_m^2} \right]^2 = \frac{2}{L} \left[ \frac{v_0^2}{r^2} + 2\sqrt{2} \frac{v_0}{r} \sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} + 2 \left(\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2}\right)^2 \right],$$

and verify equivalence with the four distinct interaction blocks of $R_v(r)$ term-by-term.

**1. The $v_0^2$ Block:**
The leading term is:

$$\frac{2}{L} \left(\frac{v_0}{r}\right)^2 = \frac{2 v_0^2}{L r^2},$$

which matches the first term of $R_v(r)$ identically.

**2. The $v_0 v_m$ Cross-Term Block:**
In the square expansion, the cross-term between $v_0$ and $v_m$ is:

$$\frac{2}{L} \left[ 2 \frac{v_0}{r} \cdot \sqrt{2} \frac{r v_m}{r^2 - a_m^2} \right] = \frac{4\sqrt{2} v_0 v_m}{L (r^2 - a_m^2)}.$$

In $R_v(r)$, this term is given by $-\frac{2\sqrt{2}}{\pi} v_0 v_m \frac{S_{\mathrm{bar}}(m)}{m}$. Recalling that $S_{\mathrm{bar}}(m) = \frac{a_m}{a_m^2 - r^2}$ and $a_m = \frac{2\pi m}{L}$, we have:

$$-\frac{2\sqrt{2}}{\pi m} S_{\mathrm{bar}}(m) = -\frac{2\sqrt{2}}{\pi m} \frac{2\pi m / L}{a_m^2 - r^2} = \frac{4\sqrt{2}}{L (r^2 - a_m^2)}.$$

Multiplying by $v_0 v_m$ shows that this block matches identically.

**3. The $v_m^2$ Diagonal Block:**
In the square expansion, the diagonal term for mode $m$ is:

$$\frac{2}{L} \cdot 2 \cdot \left(\frac{r v_m}{r^2 - a_m^2}\right)^2 = \frac{4 v_m^2}{L} \frac{r^2}{(r^2 - a_m^2)^2}.$$

In $R_v(r)$, the coefficient of $v_m^2$ is:

$$2 C_{\mathrm{bar}}(m) - \frac{1}{\pi m} S_{\mathrm{bar}}(m) = \frac{2(r^2 + a_m^2)}{L (r^2 - a_m^2)^2} - \frac{1}{\pi m} \frac{a_m}{a_m^2 - r^2}.$$

Using $a_m / (\pi m) = 2/L$, the second term becomes:

$$-\frac{1}{\pi m} \frac{a_m}{a_m^2 - r^2} = \frac{2}{L (r^2 - a_m^2)} = \frac{2(r^2 - a_m^2)}{L (r^2 - a_m^2)^2}.$$

Adding the two fractions:

$$\frac{2(r^2 + a_m^2) + 2(r^2 - a_m^2)}{L (r^2 - a_m^2)^2} = \frac{4 r^2}{L (r^2 - a_m^2)^2}.$$

This matches the diagonal term of the square identically.

**4. The $v_m v_n$ ($m < n$) Off-Diagonal Block:**
In the square expansion, the pairwise cross-term between mode $m$ and mode $n$ ($m \ne n$) is:

$$\frac{2}{L} \cdot 2 \cdot \sqrt{2} \frac{r v_m}{r^2 - a_m^2} \cdot \sqrt{2} \frac{r v_n}{r^2 - a_n^2} = \frac{8 v_m v_n}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)}.$$

In $R_v(r)$, the corresponding term is:

$$\frac{4 v_m v_n}{\pi} \frac{m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n)}{n^2 - m^2}.$$

We perform partial fractions on the numerator using $S_{\mathrm{bar}}(m) = \frac{\kappa m}{\kappa^2 m^2 - r^2}$:

$$m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n) = \frac{\kappa m^2}{\kappa^2 m^2 - r^2} - \frac{\kappa n^2}{\kappa^2 n^2 - r^2} = \kappa \left[ \frac{n^2}{r^2 - \kappa^2 n^2} - \frac{m^2}{r^2 - \kappa^2 m^2} \right].$$

Combining over a common denominator:

$$n^2 (r^2 - \kappa^2 m^2) - m^2 (r^2 - \kappa^2 n^2) = r^2 (n^2 - m^2).$$

Therefore:

$$\frac{m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n)}{n^2 - m^2} = \frac{\kappa r^2 (n^2 - m^2)}{(n^2 - m^2)(r^2 - a_m^2)(r^2 - a_n^2)} = \frac{2\pi}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)}.$$

Multiplying by $\frac{4 v_m v_n}{\pi}$:

$$\frac{4 v_m v_n}{\pi} \cdot \frac{2\pi}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)} = \frac{8 v_m v_n}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)}.$$

This matches the cross-terms of the square identically. Since all four blocks match identically, the algebraic identity $R_v(r) \equiv \frac{1}{r^2} A(1/r^2)$ is exact. $\blacksquare$

### Theorem 3.2 (Neumann Resolvent Representation and Heat-Kernel Identity)
*Let $\mathcal{L} = -\frac{d^2}{dt^2}$ denote the Neumann Laplacian on the physical interval $[0, L]$ with domain $\mathcal{D}(\mathcal{L}) = \{f \in H^2(0, L) : f'(0) = f'(L) = 0\}$. The normalized cosine eigenbasis of $\mathcal{L}$ is given by:*

$$\phi_0(t) = 1, \qquad \phi_m(t) = \sqrt{2} \cos\left(\frac{2\pi m t}{L}\right) \quad (m \ge 1),$$

*with eigenvalues $\mathcal{L} \phi_m = a_m^2 \phi_m$, where $a_m = \frac{2\pi m}{L}$. For any canonical coefficient vector $v \in \mathbb{R}^{N+1}$, the spatial trigonometric wave $T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos(a_m t)$ expands in this basis as $T_v(t) = \sum_{m=0}^N v_m \phi_m(t)$.*

*Then:*
1. **Boundary Resolvent Identity:** *The rational function $D(z)$ is identically equal to the boundary evaluation of the operator resolvent $(I + z\mathcal{L})^{-1}$ acting on the spatial profile $T_v$:*

   $$D(z) \equiv \big[(I + z\mathcal{L})^{-1} T_v\big](0) \qquad \forall z \in \mathbb{C} \setminus \{-a_1^{-2}, \dots, -a_N^{-2}\}.$$

2. **Heat-Resolvent Integral Representation:** *For $\operatorname{Re}(z) > 0$, $D(z)$ admits the exact Laplace representation against the heat evolution of the spatial wave profile:*

   $$D(z) = \int_0^\infty e^{-s} \big[ e^{-s z \mathcal{L}} T_v \big](0) \, ds.$$

3. **Boundary Resolvent Expansion:** *The Taylor coefficients of $D(z)$ around $z = 0$ are the boundary evaluations of the iterated Neumann Laplacian, reproducing the endpoint-jet hierarchy:*

   $$D(z) = \sum_{k=0}^\infty (-1)^k \big[\mathcal{L}^k T_v\big](0) z^k = \sum_{k=0}^\infty T_v^{(2k)}(0) z^k = \sum_{k=0}^\infty D_k z^k.$$

### Proof of Theorem 3.2
By spectral decomposition of the self-adjoint operator $\mathcal{L}$, the resolvent acts diagonally on the eigenbasis:

$$(I + z\mathcal{L})^{-1} \phi_m = \frac{1}{1 + a_m^2 z} \phi_m \qquad (m = 0, 1, \dots, N),$$

with $a_0 = 0$. Evaluating at the boundary $t = 0$, where $\phi_0(0) = 1$ and $\phi_m(0) = \sqrt{2}$ for all $m \ge 1$:

$$\big[(I + z\mathcal{L})^{-1} T_v\big](0) = v_0 \phi_0(0) + \sum_{m=1}^N \frac{v_m}{1 + a_m^2 z} \phi_m(0) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 z} = D(z),$$

which proves (1). 

For $\operatorname{Re}(z) > 0$, using the elementary Laplace transform identity $(1 + a_m^2 z)^{-1} = \int_0^\infty e^{-s(1 + a_m^2 z)} \, ds = \int_0^\infty e^{-s} e^{-s z a_m^2} \, ds$, we substitute into the eigenmode expansion:

$$D(z) = \sum_{m=0}^N v_m \phi_m(0) \int_0^\infty e^{-s} e^{-s z a_m^2} \, ds = \int_0^\infty e^{-s} \left[ \sum_{m=0}^N v_m e^{-s z \mathcal{L}} \phi_m \right](0) \, ds = \int_0^\infty e^{-s} \big[ e^{-s z \mathcal{L}} T_v \big](0) \, ds,$$

which establishes (2). 

> [!NOTE]
> **Domain of the Heat Representation:** The Laplace integral representation (2) holds strictly in the right half-plane $\operatorname{Re}(z) > 0$, where $e^{-s z \mathcal{L}}$ acts as a contractive heat diffusion semigroup. On the negative spectral axis $z = -1/r^2 < 0$, substituting $z = -1/r^2$ inside the integrand produces terms of the form $e^{+a_m^2 s / r^2}$, so the integral fails to converge termwise. Thus, the heat-semigroup representation belongs naturally to the positive resolvent half-plane $\operatorname{Re}(z) > 0$; analytic continuation to the negative spectral axis is governed by the meromorphic rational resolvent $D_N(z)$ rather than direct continuation of the heat integral.

Finally, expanding the resolvent as a formal geometric series $(I + z\mathcal{L})^{-1} = \sum_{k=0}^\infty (-1)^k z^k \mathcal{L}^k$, we compute the boundary evaluations of the powers of $\mathcal{L}$. Since $\mathcal{L}^k = (-1)^k \frac{d^{2k}}{dt^{2k}}$ and $T_v^{(2k)}(0) = (-1)^k \sqrt{2} \sum_{m=1}^N a_m^{2k} v_m$, we obtain the clean identity:

$$\big[\mathcal{L}^k T_v\big](0) = \sqrt{2} \sum_{m=1}^N a_m^{2k} v_m = (-1)^k T_v^{(2k)}(0),$$

which immediately yields:

$$(-1)^k \big[\mathcal{L}^k T_v\big](0) = T_v^{(2k)}(0) = D_k,$$

reproducing the Taylor expansion (3) and the endpoint-jet coefficients $D_k = T_v^{(2k)}(0)$. $\blacksquare$

---

## 4. Unconditional Finite-$N$ Positivity and the Spectral Lattice Identity

### Theorem 4.1 (Unconditional Finite-$N$ Pointwise Kernel Positivity and Entire Amplitude)
*The Fourier-side Archimedean kernel $K_{\mathrm{Fourier}}(v, r, L)$ is unconditionally pointwise non-negative on the real line for all $v \in \mathbb{R}^{N+1}$:*

$$K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0 \qquad \forall r \in \mathbb{R},$$

*where $\Phi_v(r)$ is an entire function of exponential type at most $L/2$ given by:*

$$\Phi_v(r) = \frac{2}{\sqrt{L}} \left[ v_0 \frac{\sin(rL/2)}{r} + \sqrt{2} \sum_{m=1}^{N} v_m \frac{r \sin(rL/2)}{r^2 - a_m^2} \right].$$

*Since $\Phi_v(r)$ is real for real $r$ and real $v$, this may equivalently be written as $|\Phi_v(r)|^2$ on the real axis.*

> [!IMPORTANT]
> **Pointwise Kernel Non-Negativity vs Integrated Form Positivity:** Theorem 4.1 establishes that the Fourier-side kernel $K_{\mathrm{Fourier}}(v, r, L)$ is pointwise non-negative for every real $r$. This does **not** imply that the integrated Archimedean quadratic form $\mathcal{Q}_{\mathrm{arch}}(v) = \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$ is positive semi-definite, because the smooth Archimedean density $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$ takes negative values at low frequencies ($h_+(0) \approx -5.37$). In our numerical evaluations, the integrated Archimedean contribution is negative ($\mathcal{Q}_{\mathrm{arch}} \approx -1.48$), acting as a negative barrier that counterbalances the positive zeta pole dilation energy.

### Proof of Theorem 4.1
Using the trigonometric identity $1 - \cos(rL) = 2 \sin^2(rL/2)$, we factor the full kernel:

$$K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r) = 2 \sin^2(rL/2) \cdot \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} \right]^2.$$

Bringing the envelope $\sqrt{2}\sin(rL/2)$ inside the square yields $\Phi_v(r)^2$. 

At $r = 0$, $\lim_{r\to 0} \frac{\sin(rL/2)}{r} = \frac{L}{2}$ and $\lim_{r\to 0} \frac{r\sin(rL/2)}{r^2 - a_m^2} = 0$, giving the finite limit:

$$\Phi_v(0) = \sqrt{L} v_0.$$

At the apparent poles $r = \pm a_m$, we have $a_m L / 2 = \pi m$. Taylor expansion of $\sin(rL/2)$ around $r = a_m$ gives:

$$\sin(rL/2) = \sin(\pi m + (r - a_m)L/2) = (-1)^m \sin((r - a_m)L/2) = (-1)^m \frac{L}{2}(r - a_m) + O((r - a_m)^3).$$

Because the denominator contains $r^2 - a_m^2 = (r - a_m)(r + a_m)$, both apparent poles at $r = \pm a_m$ are removable:

$$\lim_{r\to \pm a_m} \frac{r \sin(rL/2)}{r^2 - a_m^2} = (-1)^m \frac{L}{4}.$$

Because there are only finitely many apparent singularities ($r = 0$ and $r = \pm a_m$ for $m \in \{1, \dots, N\}$) and each is removable, $\Phi_v(r)$ extends to an entire function on the complex plane $\mathbb{C}$. Its exponential type is at most $L/2$. Because $v$ is real, $\Phi_v(r) \in \mathbb{R}$ for all $r \in \mathbb{R}$, which forces $\Phi_v(r)^2 \ge 0$ unconditionally on the real axis. $\blacksquare$

### Corollary 4.2 (Universal Fourier Factorization)
*The entire Fourier amplitude $\Phi_v(r)$ factors directly into the product of the universal sinc envelope and the boundary Neumann resolvent evaluated at the inverted spectral variable $z = -1/r^2$:*

$$\Phi_v(r) \equiv \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D\left(-\frac{1}{r^2}\right),$$

*and consequently the reduced rational kernel satisfies:*

$$R_v(r) \equiv \frac{2}{L r^2} \left[ D\left(-\frac{1}{r^2}\right) \right]^2.$$

### Proof of Corollary 4.2
Factoring $1/r$ from the bracketed term in Theorem 4.1:

$$\frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} = \frac{1}{r} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{r^2 v_m}{r^2 - a_m^2} \right] = \frac{1}{r} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - a_m^2 / r^2} \right].$$

Recalling the definition $D(z) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 z}$, evaluating at $z = -1/r^2$ yields precisely $D(-1/r^2)$. Substituting this into $\Phi_v(r) = \frac{2}{\sqrt{L}} \sin(rL/2) \left[ \frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} \right]$ yields the factorized amplitude $\Phi_v(r) = \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D(-1/r^2)$. Squaring and multiplying by $1 - \cos(rL) = 2\sin^2(rL/2)$ produces $R_v(r) = \frac{2}{L r^2} D(-1/r^2)^2$. $\blacksquare$

### Theorem 4.3 (Spectral Lattice Sampling Identity)
*At the discrete Fourier frequencies $a_m = 2\pi m / L$, the Archimedean kernel samples the squared Fourier coefficients orthogonally:*

$$K_{\mathrm{Fourier}}(v, 0, L) = L v_0^2 = L u_0^2,$$

$$K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 = L u_m^2 \qquad (m = 1, \dots, N).$$

### Proof of Theorem 4.3
At $r = 0$, $\Phi_v(0) = \sqrt{L} v_0$, so $K_{\mathrm{Fourier}}(v, 0, L) = \Phi_v(0)^2 = L v_0^2$. 

For $m \in \{1, \dots, N\}$, evaluate $\Phi_v(a_m)$. For every $n \ne m$, $\sin(a_m L/2) = \sin(\pi m) = 0$ while the denominator $a_m^2 - a_n^2 \ne 0$, so all non-diagonal terms vanish identically. The leading term $v_0 \frac{\sin(a_m L/2)}{a_m}$ also vanishes. 

Only the removable diagonal limit survives:

$$\Phi_v(a_m) = \frac{2}{\sqrt{L}} \cdot \sqrt{2} v_m \cdot (-1)^m \frac{L}{4} = (-1)^m \sqrt{\frac{L}{2}} v_m.$$

Squaring this amplitude:

$$K_{\mathrm{Fourier}}(v, a_m, L) = [\Phi_v(a_m)]^2 = \frac{L}{2} v_m^2.$$

Using $u_m = v_m / \sqrt{2}$, we obtain $L u_m^2$. $\blacksquare$

---

## 5. Large-$N$ Asymptotics of the Galerkin Ground State

Having established the exact finite-$N$ algebra, we now turn to the infinite-dimensional limit $N \to \infty$. We analyze the sequence of Galerkin ground states $v_N \in \mathbb{R}^{N+1}$ defined by the eigensystem:

$$Q_{c, N} v_N = \lambda_{\min}(N) v_N, \qquad \|v_N\|_2 = 1.$$

High-precision numerical diagonalizations of the Galerkin operator $Q_{c, N}$ were performed at 50 decimal digits of precision for $c = 13$ across all dimensions $N = 1, \dots, 24$ using the companion arbitrary-precision Python analysis suite (`cell34.py`, `cell40.py`, `cell41.py` with corresponding verification logs `cell34.out`, `cell40.out`, `cell41.out`) [10].

### 5.1 Observed $\ell^2$ Concentration and Numerical Convergence

**Table 1: Ground-State Mode Convergence in $\ell^2$**

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

The computed ground states exhibit strong $\ell^2$ concentration and rapidly decreasing successive increments:
- The step difference $\|v_N - v_{N-1}\|_{\ell^2}$ decreases monotonically from $0.198$ down to $0.00199$ at $N = 24$. However, monotone decrease of successive differences does not imply that their sum is finite, and thus does not mathematically establish that $(v_N)$ is a Cauchy sequence in $\ell^2$.
- What the numerical data directly show is strong mode concentration: at $N = 24$, over $99.98\%$ of the total $\ell^2$ mass is concentrated in the first five Fourier modes ($m \le 4$).

### 5.2 Observed Geometric Boundary Suppression and Conjecture 5.1

The boundary value $D_0(N) = T_{v_N}(0) = v_0 + \sqrt{2} \sum_{m=1}^N v_m$ and the second derivative $D_1(N) = T_{v_N}''(0)$ were tracked across all dimensions:

**Table 2: Geometric Boundary Decay**

| $N$ | $|D_0(N)|$ | $|D_1(N)|$ | Step Ratio $|D_0(N)| / |D_0(N-1)|$ | Decay Exponent $\alpha_N$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | $7.52 \times 10^{-3}$ | 4.8694 | — | — |
| 2 | $2.68 \times 10^{-4}$ | 0.5326 | 0.0357 | 3.333 |
| 4 | $6.00 \times 10^{-7}$ | $5.92 \times 10^{-3}$ | 0.0488 | 3.020 |
| 8 | $8.05 \times 10^{-11}$ | $3.36 \times 10^{-6}$ | 0.1434 | 1.942 |
| 12 | $6.65 \times 10^{-14}$ | $6.40 \times 10^{-9}$ | 0.2268 | 1.484 |
| 16 | $1.78 \times 10^{-16}$ | $3.13 \times 10^{-11}$ | 0.2861 | 1.251 |
| 20 | $8.38 \times 10^{-19}$ | $2.58 \times 10^{-13}$ | 0.3105 | 1.170 |
| 24 | $1.14 \times 10^{-20}$ | $5.92 \times 10^{-15}$ | 0.3244 | 1.126 |

Between $N = 1$ and $N = 24$, the boundary value drops from $7.52 \times 10^{-3}$ to $1.14 \times 10^{-20}$, i.e., by approximately **17.8 decimal orders of magnitude**. The effective decay exponent $\alpha_N = -\frac{\log(|D_0(N)|/|D_0(N-1)|)}{\log c}$ decreases from $3.333$ toward values near $1.1$ over the computed range. Thus the data do not support identifying a limiting exponent from the present range of $N$; the weaker geometric-decay conjecture below is the only asymptotic claim made here. We formulate this behavior with appropriate generality:

### Conjecture 5.1 (Geometric Boundary Suppression)
*For fixed prime cutoff $c > 1$, the boundary values of the normalized Galerkin ground states satisfy:*

$$|T_{v_N}(0)| \sim C(c) \cdot \rho(c)^N \qquad (N \to \infty),$$

*for some positive constant $C(c)$ and decay base $0 < \rho(c) < 1$. A secondary hypothesis suggests $\rho(c) \approx c^{-1/2}$, but determining the exact asymptotic base remains an open problem.*

### 5.3 Observed Eigenvalue-to-Boundary Proportionality and Numerical Conjecture 5.2

**Table 3: Ground-State Eigenvalue vs Boundary Energy**

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

While both quantities span approximately 43 decimal orders of magnitude, the ratio $\lambda_{\min}(N) / A_0(N)$ stabilizes remarkably:

$$\frac{\lambda_{\min}(N)}{A_0(N)} \approx 0.00246 \pm 0.0001 \quad (N = 18, \dots, 24).$$

For fixed $c = 13$, the ratio appears to approach a non-zero limiting constant numerically. This shows an empirical correlation while both quantities span approximately 43 decimal orders of magnitude, but does not prove asymptotic equivalence. We formulate this asymptotic relationship as a numerical conjecture:

### Conjecture 5.2 (Numerical Conjecture: Eigenvalue Gap Law)
*For a given cutoff $c > 1$, the minimum eigenvalue of the truncated Galerkin matrix is asymptotically proportional to the boundary leakage energy:*

$$\lambda_{\min}(N) \sim \kappa_c A_0(N) \equiv \frac{2 \kappa_c}{L} [T_{v_N}(0)]^2 \longrightarrow 0 \qquad (N \to \infty),$$

*where $\kappa_c = \kappa(c) > 0$ is a cutoff-dependent constant. If Conjecture 5.1 holds with secondary exponent $\rho(c) = c^{-1/2}$, then Conjecture 5.2 further predicts $\lambda_{\min}(N) \sim \widetilde{\kappa}_c \cdot c^{-N}$.*

---

## 6. The Continuum Limit: Solitary Wave and Dual Dirichlet Boundary Conditions

The continuous spatial wave profile:

$$T_{v_N}(t) = v_0 + \sqrt{2} \sum_{m=1}^{N} v_m \cos\left(\frac{2\pi m t}{L}\right)$$

was evaluated across $[0, L]$ on a dense uniform grid of 2,000 points using the companion analysis script `cell42.py` (with full profile transcript recorded in `cell42.out`) [10]. We begin by stating the exact properties valid at every finite dimension $N$, followed by the empirical continuum observations.

### Proposition 6.1 (Finite-$N$ Symmetry and Normalization)
*For every finite dimension $N \ge 1$, the trigonometric wave profile $T_{v_N}(t)$ generated by the Galerkin ground state satisfies:*

1. **Midpoint Parity Symmetry:**

   $$T_{v_N}(L - t) = T_{v_N}(t) \qquad \forall t \in [0, L].$$

2. **Exact Energy Normalization:**

   $$\|T_{v_N}\|_{L^2([0, L])}^2 = \int_0^L T_{v_N}(t)^2 \, dt = L.$$

*Proof.* Because $Q_{c, N}$ commutes with the reflection operator, its eigenspaces are reflection-invariant; in particular, an even ground-state eigenvector may always be chosen ($u_{-m} = u_m = v_m / \sqrt{2}$). (Numerically, the ground state is simple and non-degenerate for all $N$ examined). Energy normalization follows from Fourier orthogonality on $[0, L]$:

$$\int_0^L T_{v_N}(t)^2 \, dt = L \left( v_0^2 + \sum_{m=1}^N v_m^2 \right) = L \|v_N\|_2^2 = L.$$

This holds identically for all $N$. $\blacksquare$

### Numerical Observation 6.2 (Apparent Continuum Profile and Dual Dirichlet Nodes)
*Dense grid evaluations through $N = 24$ indicate that as $N$ increases, the sequence of trigonometric profiles $T_{v_N}(t)$ appears to converge to a smooth, strictly positive solitary wave $T_\infty(t)$ satisfying:*

1. **Dual Dirichlet Boundary Nodes:**

   $$T_\infty(0) = T_\infty(L) = 0.$$

2. **Interior Positivity:**

   $$T_\infty(t) > 0 \qquad \forall t \in (0, L),$$

   *with a single central maximum at $t = L/2$ of height $T_{\max} \approx 2.5382 \approx L$.*

### Conjecture 6.3 (Conjectured $C^\infty$ Boundary Flatness)
*The limiting continuum solitary wave $T_\infty(t)$ is conjectured to satisfy infinite-order flat boundary contact at both endpoints:*

$$T_\infty^{(k)}(0) = T_\infty^{(k)}(L) = 0 \qquad \forall k \ge 0.$$

*If $T_\infty$ exists with the required boundary regularity, these vanishing jets imply that the extension of $T_\infty(t)$ by zero outside $[0, L]$, denoted $\widetilde{T}_\infty(t)$, belongs to $C_c^\infty(\mathbb{R})$ with $\operatorname{supp} \widetilde{T}_\infty = [0, L]$.*

*Discussion of Convergence and Boundary Jet Extinction.*
For each fixed mode index $m$, the numerical sequence $v_{N, m}$ appears to converge to a non-zero limiting value $v_{\infty, m}$. The rapid geometric mode decay $v_{N, m} \sim c^{-m/2}$ describes decay with respect to the mode index $m$, not convergence of a fixed mode as $N \to \infty$. Justifying the interchange of differentiation with the continuum limit ($\lim_{N\to\infty} T_{v_N}^{(2k)}(0) = T_\infty^{(2k)}(0)$) requires a uniform bound of the form $|v_{N, m}| \le C q^m$ ($q < 1$) valid uniformly in $N$, or an analogous weighted $\ell^2$ estimate controlling $\sum_{m=1}^\infty m^{2k} |v_{N, m}|$.

While such uniform bounds remain to be established analytically, numerical evaluation of the even derivatives $D_k(N) = T_{v_N}^{(2k)}(0)$ for $k \in \{0, 1, 2, 3\}$ across $N \in \{8, 16, 24\}$ (computed via the arbitrary-precision script `cell43.py` and recorded in `cell43.out` [10]) shows rapid extinction consistent with geometric decay across all computed orders:
- $D_0$: $8.05 \times 10^{-11} \longrightarrow 1.78 \times 10^{-16} \longrightarrow 1.14 \times 10^{-20}$,
- $D_1$: $3.36 \times 10^{-6} \longrightarrow 3.13 \times 10^{-11} \longrightarrow 5.92 \times 10^{-15}$,
- $D_2$: $2.63 \times 10^{-2} \longrightarrow 1.37 \times 10^{-6} \longrightarrow 7.20 \times 10^{-10}$,
- $D_3$: $71.43 \longrightarrow 2.45 \times 10^{-2} \longrightarrow 3.61 \times 10^{-5}$.

By midpoint reflection symmetry, all odd derivatives vanish identically at all finite dimensions: $T_{v_N}^{(2k+1)}(0) \equiv 0$. The simultaneous geometric decay across the first four even derivatives provides strong empirical evidence, but the extrapolation to all orders remains conjectural.

### 6.4 Numerical WKB Interpretation of Boundary Suppression
Given any positive profile $T(t)$, one can formally define an effective Schrödinger potential by:

$$V_{\mathrm{conf}}(t) - E := \frac{T''(t)}{T(t)}.$$

Under this definition, $T(t)$ formally satisfies the stationary Schrödinger equation $-T''(t) + V_{\mathrm{conf}}(t) T(t) = E T(t)$ as an identity. For the computed ground-state profile, this construction produces an effective potential whose minimum lies at the midpoint $t = L/2$ and which rises steeply toward the boundaries, but does not independently establish that the underlying mathematical system is a Schrödinger differential operator.

The boundary suppression can then be modeled semiclassically via the WKB tunneling action across the barrier $[0, t_{\mathrm{turn}}]$:

$$\mathcal{S}_{\mathrm{WKB}} = \int_0^{t_{\mathrm{turn}}} \sqrt{\frac{T''(t)}{T(t)}} \, dt,$$

where $t_{\mathrm{turn}} \approx 0.4079 L$ is the effective turning point, defined here by $T''(t_{\mathrm{turn}}) = 0$.

*Numerical Comparison and WKB Barrier Computation.*
At $N = 24$, the numerical turning point is $t_{\mathrm{turn}} \approx 1.046259$ ($0.40791 L$), computed via `cell44.py` (output log `cell44.out` [10]). The WKB barrier action evaluates to:

$$\mathcal{S}_{\mathrm{WKB}} = 44.363852.$$

Comparing this with the actual boundary suppression across 20 orders of magnitude:

$$\text{Actual Suppression} = \log\left(\frac{T(L/2)}{T(0)}\right) = \log\left(\frac{2.538158}{1.137963 \times 10^{-20}}\right) = 46.853901.$$

$$\frac{\text{Actual Suppression}}{\mathcal{S}_{\mathrm{WKB}}} = \frac{46.853901}{44.363852} = 1.05613.$$

The numerically constructed effective potential yields a WKB action whose exponential scale matches the observed boundary suppression within **$5.6\%$** for a boundary suppression corresponding to approximately 20 decimal orders of magnitude. The agreement is consistent with a WKB interpretation of the observed boundary suppression, but does not by itself establish that quantum tunneling is the underlying mathematical mechanism.

### Proposition 6.5 (Legendre Expansion via Bauer–Bessel Transform)
*In normalized coordinates $x = \frac{2t}{L} - 1 \in [-1, 1]$, the finite-$N$ normalized even wave $\psi_N(x) = T_{v_N}\left(\frac{x+1}{2} L\right)$ admits an exact Legendre expansion in $L^2([-1, 1])$ (converging uniformly on $[-1, 1]$ since $\psi_N$ is smooth):*

$$\psi_N(x) = \sum_{k=0}^\infty c_{2k}^{(N)} P_{2k}(x),$$

*whose coefficients are given in exact closed analytical form via Bauer's spherical Bessel expansion:*

$$c_0^{(N)} = v_{N, 0}, \qquad c_{2k}^{(N)} = (4k + 1) \sqrt{2} (-1)^k \sum_{m=1}^N (-1)^m v_{N, m} j_{2k}(\pi m) \quad (k \ge 1),$$

*where $j_n(z) = \sqrt{\frac{\pi}{2z}} J_{n+1/2}(z)$ is the spherical Bessel function of the first kind. Conditional on the existence of a limiting coefficient vector $v_\infty \in \ell^2$, the limiting multipoles are given by the formal series $c_{2k}^{(\infty)} = (4k + 1) \sqrt{2} (-1)^k \sum_{m=1}^\infty (-1)^m v_{\infty, m} j_{2k}(\pi m)$.*

*Observed Properties of the Legendre Multipoles.*
Numerical evaluation of the Legendre expansion coefficients via the Bauer–Bessel transform (conducted in `cell44.py`, logged in `cell44.out` [10]) reveals two striking structural features:

1. **Spectral Concentration:** Truncation at $K = 10$ ($P_{20}(x)$) captures **$99.999984\%$** of the $L^2$ norm:

   $$\sum_{k=0}^{10} \frac{2}{4k + 1} [c_{2k}^{(24)}]^2 = 1.99999968 \approx 2.00000000.$$

   Over **$93.7\%$** of the wave's total energy is concentrated in the lowest four even multipoles: $P_0$ ($29.9\%$), $P_2$ ($31.4\%$), $P_4$ ($21.2\%$), and $P_6$ ($11.1\%$).
2. **Observed Alternating Phases:** Throughout the resolved multipoles ($k \le 10$), the computed coefficients exhibit an alternating-sign pattern:

   $$c_{2k}^{(N)} = (-1)^k |c_{2k}^{(N)}|.$$

   Because $P_{2k}(0) = (-1)^k \frac{(2k)!}{2^{2k}(k!)^2}$, for the resolved multipoles $k \le 10$, the observed alternating signs make the contributions at the center $x = 0$ ($t = L/2$) strictly constructive:

   $$\psi_N(0) \approx \sum_{k=0}^{10} |c_{2k}^{(N)}| \frac{(2k)!}{2^{2k}(k!)^2} > 0.$$

   Conversely, at the boundaries $x = \pm 1$ ($t = 0, L$), $P_{2k}(\pm 1) = 1$. If the limiting expansion converges uniformly at the boundary and the alternating sign pattern continues to all orders, destructive cancellation would yield:

   $$\psi_\infty(\pm 1) = \sum_{k=0}^\infty c_{2k}^{(\infty)} = |c_0^{(\infty)}| - |c_2^{(\infty)}| + |c_4^{(\infty)}| - |c_6^{(\infty)}| + \dots = 0.$$

### Conjecture 6.6 (Extinction of the Asymptotic Tail Hierarchy)
*For every fixed $k \ge 0$, the $k$-th coefficient $A_k(N)$ in the inverse-power asymptotic expansion of the Archimedean resolvent vanishes identically in the continuum limit:*

$$A_k(\infty) = \lim_{N\to\infty} \frac{2}{L} (-1)^k \sum_{j=0}^k D_j(N) D_{k-j}(N) = 0 \qquad \forall k \ge 0.$$

*If, in addition, the asymptotic expansion is sufficiently uniform in $N$ to permit the interchange of $\lim_{N\to\infty}$ and $r \to \infty$, the limiting continuous-variable resolvent:*

$$R_\infty(r) = \lim_{N\to\infty} \frac{2}{L} \left[ \frac{v_{N, 0}}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_{N, m}}{r^2 - a_m^2} \right]^2$$

*is conjectured to decay faster than every inverse power of $r$:*

$$R_\infty(r) = o(r^{-k}) \qquad \forall k \in \mathbb{N}.$$

*Logical Hierarchy of the Tail Extinction Progression.*
The mathematical status of the tail extinction decomposes into four distinct logical levels:
1. **Proven at Finite $N$:** The high-frequency expansion coefficients $A_k(N)$ of $R_{v_N}(r) = \sum_{j=0}^M A_j(N) r^{-(2j+2)} + \mathcal{O}(r^{-(2M+4)})$ are determined algebraically by the boundary Taylor jet $D_j(N) = T_{v_N}^{(2j)}(0)$ of the wave profile:
   $$A_k(N) = \frac{2}{L} (-1)^k \sum_{j=0}^k D_j(N) D_{k-j}(N).$$
2. **Numerically Observed:** High-precision evaluation confirms $A_k(N) \to 0$ across all tested orders $k \in \{0, \dots, 4\}$ and dimensions $N \in \{4, \dots, 24\}$.
3. **Conjectured:** For each fixed $k \ge 0$, $\lim_{N\to\infty} A_k(N) = 0$.
4. **Further Conjectured (Requiring an Additional Uniformity Theorem):** Conditional on uniform control over the expansion remainder in $N$ permitting the interchange of $\lim_{N\to\infty}$ and $r \to \infty$, the limiting continuous resolvent is conjectured to decay faster than every inverse power of $r$: $R_\infty(r) = o(r^{-k})$ for all $k \in \mathbb{N}$.

*High-Frequency Resolvent Asymptotics and Remainder Extinction.*
While $A_k(N) \to 0$ reflects the extinction of each individual Taylor-jet coefficient at finite $N$, establishing $R_\infty(r) = o(r^{-k})$ requires uniform control over the expansion remainder to justify interchanging $\lim_{N\to\infty}$ and $r \to \infty$. 

Numerical evaluation across $N \in \{4, 8, 12, 16, 20, 24\}$ (implemented in `cell45.py` and logged in `cell45.out` [10]) shows rapid extinction consistent with geometric decay across all computed orders:
- $A_0$: $2.81 \times 10^{-13} \to 5.05 \times 10^{-21} \to 1.01 \times 10^{-40}$ (collapsing by 27 orders of magnitude),
- $A_1$: $5.54 \times 10^{-9} \to 4.22 \times 10^{-16} \to 1.05 \times 10^{-34}$,
- $A_2$: $3.48 \times 10^{-5} \to 1.21 \times 10^{-11} \to 4.01 \times 10^{-29}$,
- $A_3$: $7.65 \times 10^{-2} \to 1.47 \times 10^{-7} \to 7.28 \times 10^{-24}$,
- $A_4$: $73.42 \to 9.23 \times 10^{-4} \to 7.53 \times 10^{-19}$.

At high frequencies, the finite-$N$ resolvent plunges precipitously: $R_{v_{24}}(10.0) = 0.0368$, $R_{v_{24}}(15.0) = 6.30 \times 10^{-6}$, $R_{v_{24}}(20.0) = 1.10 \times 10^{-8}$, and $R_{v_{24}}(50.0) = 5.40 \times 10^{-30}$. The effective logarithmic slope $\gamma_{\mathrm{eff}}(r) = -r R'(r)/R(r)$ reaches $\gamma_{\mathrm{eff}} \approx 78.6$ at $r = 15.0$, $154.0$ at $r = 20.0$, and $270.3$ at $r = 30.0$, demonstrating extremely strong finite-$N$ suppression over the computed high-frequency range (though for any fixed finite $N$, the rational tail must eventually dominate at sufficiently large $r$).

### 6.7 The Accumulating Pole Mechanism and Boundary Flatness
The operator-resolvent representation $D_N(z) = \big[(I + z\mathcal{L})^{-1} T_{v_N}\big](0)$ established in Theorem 3.2 provides the conceptual mechanism for understanding how the high-frequency Fourier tail behaves in the continuum limit, reconciling the extinction of the inverse-power tail coefficients ($A_k(N) \to 0$ for all $k \ge 0$) with a non-trivial continuous-variable resolvent $R_\infty(r)$:

1. **Accumulation of Resolvent Poles at the Origin (Proven):** At every finite dimension $N$, $D_N(z)$ is a rational function whose poles lie on the negative real axis at:

   $$z_m = -\frac{1}{a_m^2} = -\frac{L^2}{4\pi^2 m^2} \in \left(-\frac{L^2}{4\pi^2}, 0\right) \quad (m = 1, \dots, N).$$

   As $N \to \infty$, the number of poles grows without bound and their locations accumulate at the origin:

   $$\lim_{m\to\infty} z_m = 0^-.$$

2. **Conditional Obstruction to Analyticity (Analytic Mechanism):** The residue of $D_N(z)$ at each pole $z_m$ is proportional to $\frac{\sqrt{2} v_{N, m}}{a_m^2}$. If the mode coefficients $v_{N, m}$ persist with sufficient weight in the large-$N$ limit, the infinite accumulation of poles at $z = 0^-$ obstructs analytic continuation through the origin from the negative real axis.
3. **Vanishing Taylor Jet $\not\Rightarrow$ Triviality:** The numerical extinction of every tested fixed Taylor coefficient $D_k(N) = T_{v_N}^{(2k)}(0) \to 0$ is consistent with a limiting object that is $C^\infty$-flat at $z = 0$ from the right ($\operatorname{Re}(z) > 0$). However, because $z = 0$ is an accumulation boundary of singularities, this flatness does not force $D_\infty(z)$ to vanish identically on the negative axis.
4. **Motivated Asymptotic Hypothesis:** When the accumulating-pole mechanism survives in the infinite-dimensional limit, it provides a natural obstruction to ordinary Taylor analyticity, though it does not by itself select a unique essential-singularity scale. The observed WKB quantum barrier behavior in Section 6.1 motivates testing an exponentially flat ansatz:

   $$D_\infty\left(-\frac{1}{r^2}\right) \sim e^{-C r} \qquad (r \to \infty),$$

   corresponding to $D_\infty(z) \sim e^{-C / \sqrt{-z}}$ as $z \to 0^-$. Whether the true continuum resolvent settles precisely to this exponential scale, a stretched exponential, or a broader boundary-layer scaling remains an open question under numerical investigation.

### Numerical Observations from the Operator Resolvent and Discrete Cauchy Transform (Cell 51)
Dedicated numerical investigation across dimensions $N \in \{8, 12, 16, 20, 24\}$ (implemented in `cell51.py` and logged in `cell51.out` [10]) reveals five concrete structural features of $D_N(z)$:

1. **Exact Discrete Cauchy Transform Identity:**
   On the negative real axis $z = -1/r^2$, defining the quadratic lattice variable $w = -r^2 / \kappa^2$, the generating function $D_N(-1/r^2)$ matches the discrete Cauchy transform:

   $$D_N\left(-\frac{1}{r^2}\right) \equiv v_{N, 0} + \sqrt{2} w \sum_{m=1}^N \frac{v_{N, m}}{w - m^2}$$

   identically to 51 decimal digits ($|\text{diff}| = 2.67 \times 10^{-51}$ at $N = 24$). The poles on the $w$-axis are simply the positive integers $w = m^2$.

2. **Persistent Lattice Oscillations Across Pole Cells:**
   While $|D_{24}(-1/r^2)|$ drops across 14 orders of magnitude (falling from order unity to $8.38 \times 10^{-15}$ at $r \approx 55$), the ratio $-\log|D_{24}|/r$ does not settle to a single constant $C$, but oscillates between $0.37$ and $0.59$. Probing off-lattice points $r = \kappa(m + \delta)$ across $\delta \in \{0.10, 0.25, 0.50, 0.75, 0.90\}$ confirms that suppression is universal across the entire cell between poles ($\sim 10^{-12} - 10^{-13}$ at $m = 20$ for all $\delta$), ruling out half-integer sampling artifacts while confirming that a simple monotonic $e^{-Cr}$ law is not supported at finite $N$.

3. **Irregular Mode Coefficient Signs:**
   The modulated Fourier coefficients $b_{24, m} = (-1)^m v_{24, m}$ are strictly positive for $m = 1, \dots, 5$ ($0.674 \to 0.443 \to 0.213 \to 0.069 \to 0.011$), but reverse sign at $m = 6, 7, 8$ ($-9.30 \times 10^{-4}, -8.65 \times 10^{-4}, -1.14 \times 10^{-4}$) and oscillate irregularly thereafter. This disproves simple geometric decay $v_m \sim (-1)^m C q^m$, establishing that the finite-$N$ boundary suppression arises from a delicate balance between a smooth low-frequency profile and an oscillatory edge correction near $m \sim N$.

4. **Spectral-Edge Time Scale $u_N \sim (\kappa N)^{-2}$ and Resolvent Crossover:**
   For $N = 24$, the heat boundary trace $H_N(u) = [e^{-u\mathcal{L}} T_N](0)$ plunges across 20 orders of magnitude ($0.544 \to 1.77 \times 10^{-20}$), reaching the exact boundary value $T_{24}(0)$ at $u = 10^{-6}$. While this confirms that in the continuum limit $H_\infty(u) = 0$ for all $u > 0$, the inverse spectral-edge scale $u_N = a_N^{-2} = \frac{1}{\kappa^2 N^2}$ ($u_{24} \approx 2.9 \times 10^{-4}$) acts as an exact crossover scale for the resolvent integral rather than an $N$-independent universal heat profile. Dedicated investigation (Cell 52 [10]) reveals that normalized profiles $\Theta_N(s) = H_N(s u_N) / T_N(0)$ steepen systematically with $N$ due to the rapid cancellation of the boundary jet $T_N(0) / T_N''(0)$, demonstrating that the spectral-edge scale and the endpoint-cancellation scale decouple in the large-$N$ limit. Subsequent investigation (Cell 53 [10]) shows that when heat time is scaled by the first-jet cancellation scale $u_1 = D_0 / D_1$, the normalized profiles $\Theta_N^{\mathrm{cancel}}(\theta) = H_N(\theta u_1) / D_0$ exhibit a near-perfect universal collapse across all dimensions $N \in \{8, \dots, 24\}$ ($2.12 \pm 0.02$ at $\theta = 1.0$), with stable shape invariants $\beta_N = D_0 D_2 / D_1^2 \approx 0.19 - 0.26$.

5. **Resolvent Asymmetry:**
   Along the positive real axis $x > 0$, $(I + x\mathcal{L})^{-1}$ is a bounded, positive operator, and $D_N(x)$ remains $O(1)$ across all dimensions ($D_{24}(1.0) \approx 0.431$, $D_{24}(10.0) \approx 0.533$), proving that the limiting operator resolvent cannot vanish identically. Thus, $z = 0$ serves as an asymmetric boundary separating an $O(1)$ positive resolvent from a super-suppressed boundary layer at $z \to 0^+$ and an accumulating discrete Cauchy singularity at $z = -1/r^2 < 0$.

6. **Anatomy of the First-Jet Ratio $D_1/D_0$, Subexponential Scale, and Exact Commutator Algebra (Cell 54):**
   Dedicated investigation of the first-jet cancellation scale $u_1 = D_0 / D_1$ and its decoupling ratio $s_N = (\kappa N)^2 (D_0 / D_1)$ (implemented in `cell54.py` and logged in `cell54.out` [10]) reveals three exact algebraic identities and four empirical mechanisms:
   - *Exact Archimedean Resolvent Identity:* From the tail convolution $A_0 = \frac{2}{L} D_0^2$ and $A_1 = -\frac{4}{L} D_0 D_1$, the first-jet ratio satisfies the exact algebraic identity:
     $$\frac{D_1}{D_0} \equiv -\frac{1}{2} \frac{A_1}{A_0},$$
     identifying $D_1 / D_0$ directly as the relative first correction to the large-$r$ resolvent $R_v(r) = \frac{A_0}{r^2} + \frac{A_1}{r^4} + \cdots$. In terms of the generating function $F(z) = e^T (I - z M^2)^{-1} v = \int \frac{d\mu_N(x)}{1 - z x}$ for the signed spectral measure $\mu_N = v_0 \delta_0 + \sqrt{2} \sum_{m=1}^N v_m \delta_{a_m^2}$, this ratio is the logarithmic derivative at zero:
     $$-\frac{1}{2} \frac{A_1}{A_0} = \frac{F'(0)}{F(0)} = \frac{\int x \, d\mu_N(x)}{\int d\mu_N(x)}.$$
   - *Rank-4 Quadratic Commutator & Forced Moment Balance:* In full Fourier coordinates, the quadratic commutator $[M^2, Q]$ with $M = \operatorname{diag}(n)$ has rank $\le 4$:
     $$[M^2, Q] = b e^T + a p^T - p a^T - e b^T, \qquad a_n = n, \, b_n = n\psi(n), \, p_n = \psi(n).$$
     When applied to the even ground state $u$, parity eliminates all but two terms: $[M^2, Q] u = D_0 b - B_1 e$ with $B_1 = \sum n\psi(n) u_n$. Because the ground-state eigenvalue is orders of magnitude smaller than $D_0$ ($\lambda_{\min}(24) \sim 10^{-43} \ll |D_0| \sim 10^{-20}$), the quadratic spectral moment satisfies the forced linear balance:
     $$Q M^2 u \approx -D_0 b + B_1 e,$$
     proving that the entire quadratic moment $M^2 u$ is sourced by an amplitude proportional to $D_0$.
   - *Subexponentiality of the Difference $\Delta_N$:* While $-\log|D_0|$ grows from $23.24 \to 45.92$ and $-\log|D_1|$ grows from $12.60 \to 32.76$, their direct difference:
     $$\Delta_N = -\log|D_0| + \log|D_1| = \log\frac{|D_1|}{|D_0|}$$
     drifts only from $10.64$ to $13.16$ across $N \in \{8, \dots, 24\}$, while the consecutive effective decay rates narrow monotonically toward each other ($|\alpha_N(D_0) - \alpha_N(D_1)|$ drops $0.2088 \to 0.1499 \to 0.1408 \to 0.1308$). Because $s_N = (\kappa N)^2 (D_0 / D_1) \iff D_0 / D_1 = s_N / (\kappa N)^2$, the subexponential behavior of $s_N$ implies that $D_0$ and $D_1$ share the same underlying leading exponential barrier suppression factor.
   - *Structured Signed Cancellation:* Mode-by-mode decomposition confirms cancellations down to $10^{-20}$ ($D_0$) and $10^{-15}$ ($D_1$) between $\mathcal{O}(1)$ sub-sums $S^\pm$.
   - *Asymmetry in Bulk vs. Edge Decoupling:* $D_0$ is governed by bulk cancellation with an exponentially tiny edge correction ($-1.52 \times 10^{-8}$ at $N = 24$), whereas $D_1$ is governed by an exact bulk-edge balance ($\mathrm{bulk} + \mathrm{edge} \approx 0$).
   - *Sobolev Trace Non-Sharpness:* The Cauchy–Schwarz bound $|D_1| \le \sqrt{2}\kappa \|T'_v\|_{L^2} \sqrt{\sum m^2}$ grows as $\sim N^{3/2}$ while $D_1 \to 0$, giving a non-sharpness ratio of $7.58 \times 10^{-18}$ at $N = 24$. Because $\|T_N\|_{L^2} = 1$ while $T_N(0) \to 0$, endpoint suppression is encoded specifically in the variational ground-state eigenvector rather than generic Sobolev norm constraints.

### Conditional Vanishing of the Volterra Boundary Jump
If the conjectured infinite-order boundary flatness holds (Conjecture 6.3), then the corresponding Volterra convolution:

$$K_\infty(\omega) = 2 \int_0^\omega T_\infty(t) T_\infty(\omega - t) \, dt$$

vanishes smoothly at both $\omega = 0$ and $\omega = 1$:

$$\lim_{\omega \to 0} K_\infty(\omega) = 0, \qquad \lim_{\omega \to 1} K_\infty(\omega) = 0,$$

with no jump discontinuities of any finite order at $\omega = 1$. This would eliminate the boundary jump that historically produced the oscillatory factor $1 - \cos(rL)$ and the $A_0/r^2$ tail in the finite-rank Galerkin models.

### Proposition 6.8 (Tri-Partite Decomposition and Continuous-Quadrature Balance)
*Let $\mathcal{Q}(v) = \mathcal{Q}_{\mathrm{pole}}(v) + \mathcal{Q}_{\mathrm{prime}}(v) + \mathcal{Q}_{\mathrm{arch}}(v)$ be the Connes–van Suijlekom quadratic form on the Galerkin subspace of dimension $2N+1$. For every finite dimension $N$, the algebraic matrix sum matches the minimum eigenvalue identically by definition:*

$$\mathcal{Q}_{\mathrm{matrix}}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{\mathrm{matrix}}(v_N) \equiv \lambda_{\min}(N).$$

*Numerical Evidence for Continuum Equilibrium.*
What is significant is the **independent numerical quadrature of the continuous Archimedean component**: when $\mathcal{Q}_{\mathrm{arch}}(v_N)$ is evaluated independently by numerical quadrature of the continuous-variable integral $\frac{1}{\pi} \int_0^{80} h_+(r) \Phi_{v_N}(r)^2 \, dr$ (using the companion analysis script `cell46.py` and logged in `cell46.out` [10]), the algebraic pole and prime contributions and the independently quadrature-evaluated Archimedean contribution cancel to a residual of order $10^{-43}$:

1. **Stabilization of the Continuous Archimedean Integral:**
   Because $R_{v_{24}}(r)$ is extremely strongly suppressed over the computed range, the continuous Archimedean integral:

   $$A_{\mathrm{arch}}(R_{\max}) = \frac{1}{\pi} \int_0^{R_{\max}} h_+(r) \Phi_{v_{24}}(r)^2 \, dr$$

   freezes completely without truncation remainder as $R_{\max}$ increases:
   - $R_{\max} = 10$: $-1.480396530465$,
   - $R_{\max} = 20$: $-1.479797764647$ (tail increment $5.99 \times 10^{-4}$),
   - $R_{\max} = 40$: $-1.479797763974798$ (tail increment $2.68 \times 10^{-16}$),
   - $R_{\max} = 60$: $-1.479797763974798$ (tail increment $4.49 \times 10^{-29}$),
   - $R_{\max} = 80$: $-1.479797763974798326397825$ (tail increment $7.57 \times 10^{-40}$).

2. **Dimension-by-Dimension Continuous-Quadrature Sum:**
   Across all Galerkin dimensions $N \in \{4, 8, 12, 16, 20, 24\}$, the three independent components cancel from $\mathcal{O}(1)$ down to the order of $\lambda_{\min}(N)$:
   - $N = 4$: $\mathcal{Q}_{\mathrm{pole}} = +2.206186$, $\mathcal{Q}_{\mathrm{prime}} = -0.316153$, $\mathcal{Q}_{\mathrm{arch}} = -1.890032$, summing to $\mathcal{Q}_{\mathrm{total}} = 7.82 \times 10^{-15}$ (matrix eigenvalue $\lambda_{\min} = 8.83 \times 10^{-15}$),
   - $N = 8$: $\mathcal{Q}_{\mathrm{pole}} = +1.813949$, $\mathcal{Q}_{\mathrm{prime}} = -0.154916$, $\mathcal{Q}_{\mathrm{arch}} = -1.659033$, summing to $\mathcal{Q}_{\mathrm{total}} = 5.38 \times 10^{-23}$ (matrix eigenvalue $\lambda_{\min} = 6.71 \times 10^{-23}$),
   - $N = 12$: $\mathcal{Q}_{\mathrm{pole}} = +1.675166$, $\mathcal{Q}_{\mathrm{prime}} = -0.108101$, $\mathcal{Q}_{\mathrm{arch}} = -1.567065$, summing to $\mathcal{Q}_{\mathrm{total}} = 1.32 \times 10^{-29}$ (matrix eigenvalue $\lambda_{\min} = 1.78 \times 10^{-29}$),
   - $N = 16$: $\mathcal{Q}_{\mathrm{pole}} = +1.609630$, $\mathcal{Q}_{\mathrm{prime}} = -0.088194$, $\mathcal{Q}_{\mathrm{arch}} = -1.521436$, summing to $\mathcal{Q}_{\mathrm{total}} = 5.11 \times 10^{-35}$ (matrix eigenvalue $\lambda_{\min} = 7.12 \times 10^{-35}$),
   - $N = 20$: $\mathcal{Q}_{\mathrm{pole}} = +1.572288$, $\mathcal{Q}_{\mathrm{prime}} = -0.077529$, $\mathcal{Q}_{\mathrm{arch}} = -1.494759$, summing to $\mathcal{Q}_{\mathrm{total}} = 8.81 \times 10^{-40}$ (matrix eigenvalue $\lambda_{\min} = 1.32 \times 10^{-39}$),
   - $N = 24$: $\mathcal{Q}_{\mathrm{pole}} = +1.551652$, $\mathcal{Q}_{\mathrm{prime}} = -0.071854$, $\mathcal{Q}_{\mathrm{arch}} = -1.479798$, summing to $\mathcal{Q}_{\mathrm{total}} = 1.29 \times 10^{-43}$ (matrix eigenvalue $\lambda_{\min} = 2.53 \times 10^{-43}$).

   The residual discrepancy between the continuous-quadrature sum ($1.29 \times 10^{-43}$) and the matrix eigenvalue ($2.53 \times 10^{-43}$) is of order $10^{-43}$; while small in absolute terms, it represents a proportional factor of $\approx 1.96$ at the residual scale whose precise numerical origins (such as continuous quadrature tolerances versus matrix truncation parameters) remain to be isolated in dedicated verification benchmarks.

3. **Observed Continuum Equilibrium Candidates ($c = 13$):**
   The $N = 24$ values appear numerically stabilized and suggest the following candidate continuum values:
   $$\mathcal{Q}_{\mathrm{pole}}^{(24)} \approx +1.55165219571747,$$
   $$\mathcal{Q}_{\mathrm{prime}}^{(24)} \approx -0.07185443174267,$$
   $$\mathcal{Q}_{\mathrm{arch}}^{(24)} \approx -1.47979776397480,$$
   producing an apparent numerical balance at $N = 24$:
   $$\frac{\mathcal{Q}_{\mathrm{pole}}^{(24)}}{|\mathcal{Q}_{\mathrm{prime}}^{(24)}| + |\mathcal{Q}_{\mathrm{arch}}^{(24)}|} \approx 1.00000000000000,$$
   with values that appear stabilized over the computed dimensions.

4. **Prime-Power Decomposition of the Negative Barrier:**
   Direct point-evaluation of the Volterra convolution $K_{v_{24}}(\omega_q)$ at all prime powers $q \le 13$ matches the matrix-computed prime form to 52 decimal digits ($|\text{diff}| = 1.67 \times 10^{-52}$). The lowest prime $q = 2$ provides **$98.65\%$** of the entire prime energy ($-0.0708858$), $q = 3$ accounts for **$1.34\%$** ($-0.0009658$), while contributions above $q = 7$ decay exponentially below $10^{-13}$ ($q = 11$: $-9.52 \times 10^{-28}$). At the Volterra endpoint $\omega = 0$ ($q = 13$), $K_{v_{24}}(0) = 0$ identically. $\blacksquare$

### Numerical Observation 6.9 (Numerical Multi-$c$ Observations)
*Investigation across prime cutoffs $c \in \{5, 7, 11, 13, 17\}$ and dimensions $N \in \{4, 8, 12, 16, 20\}$ using the multi-parameter analysis suite (`cell47.py` and output log `cell47.out` [10]) reveals three striking empirical patterns governing the ground state:*

1. **Approximate Constancy of the Scaling Ratio $\kappa$ Across Cutoffs:**
   *Across prime cutoffs $c \ge 7$ at $N = 20$, the scaling ratio $\kappa_c(N) = \lambda_{\min}(N) / A_0(N)$ shows remarkable stability:*

   $$\kappa_7 = 0.0024026, \quad \kappa_{11} = 0.0023670, \quad \kappa_{13} = 0.0024145, \quad \kappa_{17} = 0.0023362.$$

   *While the ground-state eigenvalue $\lambda_{\min}(20)$ drops across 17 orders of magnitude (from $6.85 \times 10^{-27}$ at $c = 7$ to $1.15 \times 10^{-43}$ at $c = 17$), the values of $\kappa_c(20)$ remain within approximately $1.6\%$ of their mean:*

   $$\kappa \approx 0.00238 \pm 0.00004.$$

   *This suggests, but does not establish, approximate cutoff-independence of $\kappa$ across prime cutoffs $c \ge 7$.*

2. **Observed WKB Semiclassical Scaling at $N = 20$:**
   *The WKB barrier tunneling action $\mathcal{S}_{\mathrm{WKB}}(N, c) = \int_0^{t_{\mathrm{turn}}} \sqrt{T''/T} \, dt$ closely tracks the semiclassical relation:*

   $$\frac{\mathcal{S}_{\mathrm{WKB}}(N, c)}{L} \approx \frac{\pi N}{4}.$$

   *At $N = 20$, $\frac{\pi \times 20}{4} = 5\pi \approx 15.70796$. Numerical evaluations yield:*
   - $c = 11$: $\mathcal{S}_{\mathrm{WKB}} / L = 15.3258$,
   - $c = 13$: $\mathcal{S}_{\mathrm{WKB}} / L = 15.6681$ (*$99.75\%$ match to $5\pi$*),
   - $c = 17$: $\mathcal{S}_{\mathrm{WKB}} / L = 15.8090$ (*$99.36\%$ match to $5\pi$*).

   *The three tested cutoffs give values within approximately $2.5\%$ of $\pi N / 4$. This suggests a possible semiclassical scaling relation requiring further testing in both $N$ and $c$. The ratio $\text{Actual Suppression} / \mathcal{S}_{\mathrm{WKB}}$ decreases monotonically across the tested cutoffs, from $1.121$ to $1.054$ ($1.121 \to 1.084 \to 1.063 \to 1.059 \to 1.054$), moving toward unity. Across approximately 47 decimal orders of magnitude at $c = 17$, the WKB action agrees with the observed logarithmic boundary suppression to within $5.3\%$. The effective turning point appears to stabilise near $t_{\mathrm{turn}} / L \approx 0.41$ across the tested cutoffs.*

3. **Prime Energy Share Across Cutoffs:**
   *For every cutoff $c$, the matrix decomposition satisfies the algebraic balance $\mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}} = \lambda_{\min}(20) \sim 10^{-17}\text{ to }10^{-44}$. The fraction of negative dispersive energy shouldered by the discrete prime powers $f_{\mathrm{prime}}(c) = |\mathcal{Q}_{\mathrm{prime}}| / \mathcal{Q}_{\mathrm{pole}}$ is monotonically increasing over the tested cutoffs:*
   - $c = 5$: $2.79\%$ prime / $97.21\%$ arch,
   - $c = 7$: $3.42\%$ prime / $96.58\%$ arch,
   - $c = 11$: $4.47\%$ prime / $95.53\%$ arch,
   - $c = 13$: $4.93\%$ prime / $95.07\%$ arch,
   - $c = 17$: $5.76\%$ prime / $94.24\%$ arch.

   *As the scaling interval $[0, \log c]$ expands to encompass higher primes, the discrete prime-power sum shoulders an increasing fraction of the negative energy counterbalancing the geometric dilation pole.*

---

## 7. Conclusion and Analytical Roadmap toward Weil Positivity

The findings of this paper resolve the finite-$N$ algebraic reduction of the Archimedean kernel and substantially clarify its numerical tail behaviour in relation to the continuous Weil quadratic form:

1. **Exact Resolvent Formula and Pointwise Kernel Positivity (Theorem):**
   The finite-$N$ Archimedean Fourier kernel is an exact non-asymptotic square on the real axis:
   $$K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0 \qquad \forall r \in \mathbb{R}.$$
   This establishes algebraically, independently of numerical quadrature, that the Fourier-side Archimedean kernel is pointwise non-negative for every real coefficient vector. It does not by itself imply positivity of the integrated Archimedean quadratic form, since the Archimedean density $h_+(r)$ is not everywhere non-negative.

2. **The Galerkin Cutoff as a Confinement Barrier (Observation & Conjecture):**
   The finite-rank spectral gap $\lambda_{\min}(N) > 0$ is an artifact of band limitation. The finite rank $N$ prevents the trigonometric polynomial from satisfying the Dirichlet boundary condition $T(0) = 0$ identically. The boundary energy leaks out as $A_0(N) = \frac{2}{L} [T_{v_N}(0)]^2$, driving the observed numerical eigenvalue gap:
   $$\lambda_{\min}(N) \sim \kappa_c A_0(N) \longrightarrow 0 \quad \text{(Numerical Conjecture)}.$$

3. **Conjectured Emergence of Smooth Compact Support (Conjecture):**
   Numerical evidence indicates that as $N \to \infty$, the ground state $T_\infty(t)$ develops infinite-order flat boundary contact ($\widetilde{T}_\infty \in C_c^\infty(\mathbb{R})$ with $\operatorname{supp} \widetilde{T}_\infty = [0, L]$), modeled semiclassically by WKB quantum barrier penetration ($\mathcal{S}_{\mathrm{WKB}} \approx 44.36$). Conditional on this boundary flatness, all boundary jumps in the Volterra kernel vanish, and the continuous resolvent is conjectured to decay faster than every inverse power of $r$.

4. **Observed Tri-Partite Zero-Energy Balance (Observation):**
   On finite Galerkin subspaces, the matrix quadratic form decomposes into pole, prime, and Archimedean terms. When the Archimedean contribution is evaluated independently via continuous-variable quadrature, the independently computed $\mathcal{O}(1)$ terms cancel to a residual of order $10^{-43}$ at $N = 24$:

   $$\mathcal{Q}_{\mathrm{pole}}^{(24)} + \mathcal{Q}_{\mathrm{prime}}^{(24)} + \mathcal{Q}_{\mathrm{arch}}^{(24)} = 1.29 \times 10^{-43},$$

   which is of the same $10^{-43}$ scale as the matrix eigenvalue $2.53 \times 10^{-43}$, with the residual discrepancy discussed above. The positive geometric dilation energy from the zeta pole ($+1.55165$) is counterbalanced by the combined dispersive negative contributions of the prime powers ($-0.07185$) and Archimedean places ($-1.47980$).

### Reconnection with the Riemann Hypothesis (Weil Positivity)

In André Weil's 1952 explicit formula framework and Alain Connes' noncommutative geometry formulation:
- The Riemann Hypothesis is mathematically equivalent to the **positivity of the Weil quadratic form** $\Delta_{\mathrm{Weil}}(f, f) \ge 0$ on the space of test functions on the idele class group $\mathbb{A}_{\mathbb{Q}} / \mathbb{Q}^\times$.
- In the Connes–van Suijlekom truncation, the quadratic form is regularized on a finite scaling interval $[0, L] = [0, \log c]$ with cutoff $N$. For the finite dimensions $N$ and prime cutoffs $c$ examined numerically, the ground-state eigenvalue is strictly positive: $\lambda_{\min}(N) > 0$.
- As $N$ increases, numerical evidence shows that $\lambda_{\min}(N)$ decreases rapidly toward zero and is empirically proportional to the boundary leakage energy $A_0(N)$ (approaching zero strictly from above).
- **The Defensible Mathematical Core:** The finite-$N$ Archimedean contribution has an exact positive-square representation, and numerical evidence indicates that the smallest Galerkin eigenvalue is driven toward zero by increasingly strong suppression of a boundary leakage term.
- **The Analytical Roadmap toward Continuous Weil Positivity:** Bridging these finite-dimensional results to a formal proof of Weil positivity requires addressing three major analytical steps:
  1. *Operator Convergence and Spectral Pollution:* Proving an operator convergence theorem (such as strong resolvent convergence) connecting the sequence of finite-rank Galerkin operators $Q_{c, N}$ to a continuous self-adjoint operator $Q_c$ on $L^2([0, L])$. Positivity of finite Galerkin projections $\lambda_{\min}(N) > 0$ does not by itself rule out spectral pollution or negative spectrum in the infinite-dimensional limit.
  2. *Proof of Conjectures 5.1–6.6:* Establishing uniform bounds on the mode coefficients $|v_{N, m}|$ to justify term-by-term differentiation, thereby proving infinite-order boundary flatness $\widetilde{T}_\infty \in C_c^\infty(\mathbb{R})$, eliminating the Volterra boundary jump, and confirming that zero is an isolated edge mode.
  3. *The Double Scaling Limit:* Controlling the joint limit $(N, c) \to \infty$ to extend positivity from compact scaling intervals $[0, \log c]$ to the full, unbounded idele class group.

---

## 8. Computational Reproducibility and Software Availability

To ensure complete computational transparency and reproducibility, the entire mathematical software pipeline and all raw high-precision calculation transcripts supporting this study are permanently archived and openly accessible in the public repository [10]:

> **Software Repository:** <https://github.com/akivag613/connes-cvs->  
> **Mirror Repository:** <https://github.com/nrensen/connes-cvs->

The numerical calculations reported in this manuscript were performed using Python and the `mpmath` arbitrary-precision arithmetic library, with working precisions stated in the corresponding scripts and verification logs. The paper's empirical observations, asymptotic fits, and spectral decompositions are reproduced by standalone computational scripts (`cell*.py`), whose complete numerical output transcripts are preserved in matching log files (`cell*.out`). 

**Table 4: Mapping of Manuscript Results to Computational Scripts and Output Logs**

| Manuscript Section & Result | Mathematical / Numerical Focus | Python Script | Verification Log |
| :--- | :--- | :--- | :--- |
| Section 2.2 & 3 | Four-term Volterra reduction & exact resolvent identity | `cell32.py` | `cell32.out` |
| Section 5.1–5.3 (Table 1 & 2) | Ground state eigensystem & mode decay ($N=1\dots 24$) | `cell34.py`, `cell40.py`, `cell41.py` | `cell34.out`, `cell40.out`, `cell41.out` |
| Section 6.1–6.3 (Proposition 6.1, Observation 6.2, Conjecture 6.3) | Spatial wave profile & boundary jet derivatives $D_0\dots D_3$ | `cell42.py`, `cell43.py` | `cell42.out`, `cell43.out` |
| Section 6.4–6.5 (Observation 6.4, Proposition 6.5) | Effective WKB barrier potential & Bauer–Bessel Legendre multipoles | `cell44.py` | `cell44.out` |
| Section 6.6 (Conjecture 6.6) | High-frequency resolvent decay & Taylor jet extinction $A_0\dots A_4$ | `cell45.py` | `cell45.out` |
| Section 6.7 (Proposition 6.7, Cells 51–54) | Continuous Archimedean integral, Cauchy transform, dual-scale boundary layer & first-jet commutator anatomy | `cell46.py`, `cell51.py`, `cell52.py`, `cell53.py`, `cell54.py` | `cell46.out`, `cell51.out`, `cell52.out`, `cell53.out`, `cell54.out` |
| Section 6.8 (Observation 6.8) | Universal ratio $\kappa$, multi-cutoff WKB scaling & prime partition | `cell47.py` | `cell47.out` |

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
