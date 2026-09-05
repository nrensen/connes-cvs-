# An Exact Resolvent and Commutator Toolkit for the Truncated Connes–van Suijlekom Weil Quadratic Form
### Fourier–Resolvent Factorization, Reflected Autocorrelation, and Commutator Algebra

**Authors:** Research Record / Connes–CvS Investigation Series  
**Date:** September 2026  
**Software & Reproducibility Suite:** `https://github.com/akivag613/connes-cvs-` (mirror: `nrensen/connes-cvs-`)  
**Status:** Standalone Manuscript / The Rigorous Toolkit (Companion to Paper 4B: *The Research Programme*)

---

### Abstract

The truncated Weil quadratic form developed by Connes–van Suijlekom (2025) and Connes–Consani–Moscovici (2026) discretizes the explicit formula of prime number theory using a finite-rank Galerkin projection of band $N$ on a logarithmic scaling interval $[0, L] = [0, \log c]$. In finite-dimensional computational implementations of this model, the omitted Archimedean tail has previously been treated by oscillatory numerical quadrature or empirical asymptotic inverse-power fits.

In this paper, we establish the **exact algebraic solution** to the finite-$N$ Archimedean kernel and the associated commutator algebra, proving unconditionally and independently of numerical quadrature:

1. **Exact Rational Resolvent & Operator Identity (Theorem):** Starting from the four-term analytic reduction of the Archimedean Volterra integral, we prove algebraically that the reduced Fourier kernel $R_v(r) = K_{\mathrm{Fourier}}(v, r, L) / (1 - \cos(rL))$ is identically equal to the squared Cauchy resolvent:
   $$R_v(r) \equiv \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^{N} \frac{r v_m}{r^2 - a_m^2} \right]^2, \qquad a_m = \frac{2\pi m}{L},$$
   on the punctured complex plane $\mathbb{C} \setminus \{0, \pm a_1, \dots, \pm a_N\}$. Furthermore, introducing the Neumann Laplacian $\mathcal{L} = -d^2/dt^2$ on $[0, L]$ with $T'(0) = T'(L) = 0$, the rational generating function $D(z)$ is the boundary evaluation of the operator resolvent:
   $$D(z) \equiv \big[(I + z\mathcal{L})^{-1} T_v\big](0) = \int_0^\infty e^{-s} \big[ e^{-sz\mathcal{L}} T_v \big](0) \, ds \quad (\operatorname{Re} z > 0),$$
   whose Taylor expansion around $z = 0$ reproduces the entire endpoint-jet hierarchy $D_k = T_v^{(2k)}(0)$.
2. **Universal Fourier Factorization and Unconditional Pointwise Positivity (Theorem):** The entire Fourier-side amplitude $\Phi_v(r)$ factors directly in terms of the boundary Neumann resolvent evaluated at the inverted spectral variable $z = -1/r^2$:
   $$\Phi_v(r) \equiv \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D\left(-\frac{1}{r^2}\right), \qquad K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 = \frac{4}{L} \frac{\sin^2(rL/2)}{r^2} D\left(-\frac{1}{r^2}\right)^2 \ge 0,$$
   proving algebraically and unconditionally that the Fourier-side Archimedean kernel $K_{\mathrm{Fourier}}(v, r, L)$ is pointwise non-negative for all real $r$ and all real coefficient vectors $v \in \mathbb{R}^{N+1}$.
3. **Spectral Lattice Sampling Identity (Theorem):** At the lattice nodes $r = a_m$, the apparent poles cancel cleanly against the envelope zeros via removable singularities, yielding the exact sampling identity:
   $$K_{\mathrm{Fourier}}(v, 0, L) = L v_0^2 = L u_0^2, \qquad K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 = L u_m^2 \quad (m = 1, \dots, N),$$
   sampling the kernel directly into the squared Fourier coefficients.
4. **Exact Archimedean Cauchy Transform and Closed-Form Pole Decomposition (Theorems):** We evaluate the continuous Archimedean Cauchy transform $J(q) = \frac{1}{\pi}\int_0^\infty \frac{2q}{q^2 + r^2} K_{\mathrm{Fourier}}(v, r, L) \, dr$ in exact closed algebraic form via spatial Laplace transform of the reflected autocorrelation kernel (with an independent verification via contour integration in the complex frequency plane), isolating the origin residue ($2v_0^2/q$), the discrete lattice poles ($2qv_m^2/(q^2+a_m^2)$), and the boundary leakage term. Combined with the Weierstrass partial fraction expansion of the digamma function, this expresses the continuous Archimedean quadratic form $\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty [ \frac{\|v\|_2^2}{n+1} - J(q_n) ]$ (with $C_{\mathrm{arch}} = -\gamma - \log \pi$) as an unconditionally convergent algebraic series with fast $\mathcal{O}(n^{-2})$ summand convergence, yielding an explicit $\mathcal{O}(M^{-1})$ finite-$M$ remainder bound and eliminating the need for numerical quadrature.
5. **Finite-$T$ Archimedean Cutoff Defect and Endpoint-Jet Representation (Theorem):** In finite-rank Galerkin implementations where the continuous Archimedean density is evaluated with a finite integration cutoff $T > 0$, we prove the exact identity $v^T Q_{\mathrm{arch}}^{(T)} v \equiv \frac{1}{\pi}\int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$. Consequently, the discrepancy between the Galerkin eigenvalue $\lambda_N$ and the continuous tripartite functional $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N)$ is identically the Archimedean cutoff tail $\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) \equiv -\frac{1}{\pi}\int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr \equiv -\delta_T^{\mathrm{tail}}(v_N)$, resolving the exact relationship for the continuous Weil functional evaluated on the finite Galerkin test vector $v_N$. For $T > a_N = 2\pi N/L$, this tail admits an exact Laurent series $\sum_{k=0}^\infty A_k(N) \mathcal{J}_k(T, L)$ in terms of endpoint derivatives $D_j = T_v^{(2j)}(0)$ and universal moment integrals $\mathcal{J}_k(T, L)$, converging geometrically with step ratio $(a_N/T)^2$. This proves that the observed $10^{-43}$ numerical residual for $v_N$ is identically continuous Archimedean cutoff tail leakage $\delta_T^{\mathrm{tail}}(v_N)$, rigorously separating the finite cutoff artifact from the infinite-dimensional continuum minimizer projection problem.
6. **Exact Commutator Algebra and Parity Factorization (Theorem):** For the coordinate operator $M = \operatorname{diag}(n)$ and the Galerkin matrix $Q$, the commutator $[M^k, Q]$ has rank at most $2k$. Under the positive-gap condition $\lambda < \min \operatorname{spec}(Q_{\mathrm{odd}})$, parity reflection decouples the system into even and odd sectors, yielding the exact odd-sector resolvent identity:
   $$M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi, \qquad B_1 = -D_0 \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$
7. **Exact First-Jet Resolvent Identity and Small-Denominator Cancellation (Theorems):** We prove algebraically that the first-jet cancellation ratio is identically equal to the relative first correction of the large-$r$ Archimedean resolvent:
   $$\frac{D_1}{D_0} \equiv -\frac{1}{2} \frac{A_1}{A_0} \equiv \frac{D'(0)}{D(0)} \equiv -\kappa^2 \frac{F'(0)}{F(0)} \equiv -\frac{\sqrt{2}\displaystyle\sum_{m=1}^N a_m^2 v_m}{v_0 + \sqrt{2}\displaystyle\sum_{m=1}^N v_m}.$$
   In terms of the discrete signed atomic measure $\mu_N = v_0 \delta_0 + \sqrt{2} \sum_{m=1}^N v_m \delta_{a_m^2}$, this ratio equals $-\frac{\int_0^\infty x \, d\mu_N(x)}{\int_0^\infty d\mu_N(x)}$ as an algebraic quotient of finite sums. Furthermore, we establish the exact small-denominator cancellation theorem: in the spectral expansion of the resolvent $(Q_{\mathrm{even}} - \lambda I)^\dagger$ coupled to the odd arithmetic energy, the apparent denominator $(E_k - \lambda)$ for any even eigenstate $k$ cancels identically against the odd resolvent difference in the numerator via the first resolvent identity away from odd resonances.

---

## 1. Introduction

The explicit formula of Guinand and Weil relates the nontrivial zeros of the Riemann zeta function $\zeta(s)$ to arithmetic prime-power sums, pole contributions, and Archimedean gamma-factor terms. In André Weil's formulation (1952), the Riemann Hypothesis (RH) is equivalent to the non-negativity of the associated quadratic functional:

$$W(g) \ge 0$$

on all admissible test functions $g = f * f^*$ on the idele class group.

In Alain Connes’ non-commutative geometry program, Weil positivity is realized through an operator-theoretic spectral framework on prolate spheroidal wave spaces. Recent work by Connes and van Suijlekom (2025) and Connes, Consani, and Moscovici (2026) discretizes the continuous Weil form using a finite-rank Galerkin projection: for a logarithmic prime cutoff $L = \log c$ (with $c > 1$) and a finite frequency band $N \ge 1$, the continuous form is projected onto an explicit $(2N+1) \times (2N+1)$ matrix $Q_{c, N}$.

The total quadratic form decomposes into prime, pole, and Archimedean components:

$$\langle v, Q v \rangle = \langle v, Q_{\mathrm{prime}} v \rangle + \langle v, Q_{\mathrm{pole}} v \rangle + \langle v, Q_{\mathrm{arch}} v \rangle.$$

While the prime and pole contributions admit exact closed finite algebraic representations, the Archimedean term involves an integral over the positive real axis:

$$\langle v, Q_{\mathrm{arch}} v \rangle = \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr,$$

where

$$h_+(r) = \operatorname{Re} \psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi$$

is the smooth Archimedean density.

### The Role of the Toolkit and Relationship to Connes–van Suijlekom

In finite-dimensional numerical implementations of the Connes–van Suijlekom model, evaluating the Archimedean integral has previously been performed using oscillatory numerical quadrature or empirical inverse-power Laurent fits $\sum A_k / r^{2k+2}$. Furthermore, whether $K_{\mathrm{Fourier}}(v, r, L)$ possessed sign-oscillations that could induce negative eigenvalues at large frequencies remained an open question.

The Connes–van Suijlekom Galerkin discretization [2] provides the underlying divided-difference framework and basis projections. The primary contribution of this manuscript is to establish the explicit Fourier–resolvent factorization:

$$\boxed{\text{CvS finite Fourier amplitude} \implies \text{explicit squared amplitude} \implies \text{Neumann-resolvent formula} \implies \text{closed } J(q),}$$

establishing exact closed-form algebraic solutions for the Archimedean Cauchy transform, the Weierstrass pole decomposition, and the low-rank commutator projections. No novelty is claimed for the general autocorrelation/Fourier-square principle itself.

We emphasize that **every theorem and identity established in this paper is an exact, unconditional finite-$N$ mathematical statement**. We do not assert that finite-$N$ kernel non-negativity proves Weil positivity or the Riemann Hypothesis. All asymptotic questions concerning the infinite-dimensional limit $N \to \infty$ (such as solitary wave profiles, semiclassical WKB barrier tunneling, and formal Wiener–Hopf continuum scaling) are systematically investigated in the companion paper, *The Dirichlet Continuum Limit, Barrier Mechanics, and Asymptotic Weil Positivity in the Connes–van Suijlekom Galerkin Truncation* [Paper 4B].

The toolkit established herein proves:
- The exact squared rational resolvent identity for $K_{\mathrm{Fourier}}$, proving pointwise non-negativity globally on $\mathbb{R}$.
- The operator representation of the generating function as a Neumann resolvent.
- The closed-form evaluation of the Archimedean Cauchy transform $J(q)$ via reflected autocorrelation and the unconditionally convergent Weierstrass pole series for $\mathcal{Q}_{\mathrm{arch}}(v)$, completely eliminating numerical quadrature.
- The finite-$T$ Archimedean cutoff defect identity and geometric endpoint-jet Laurent series, resolving the exact relationship for the Weil functional evaluated on the finite Galerkin test vector $v_N$.
- The exact rank-$2k$ commutator algebra of the Galerkin matrix, yielding explicit closed-form representations for coordinate moments and the first-jet ratio.
- The exact algebraic cancellation of small bound-state denominators in the odd-even resolvent coupling away from resonances.

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

### 2.2 The Fourier-Side Representation and Spatial Autocorrelation

Conceptually, the Fourier-side Archimedean kernel represents the power spectral density of the compactly supported trigonometric wave $T_v(t) \mathbf{1}_{[0, L]}(t)$:

$$K_{\mathrm{Fourier}}(v, r, L) \equiv \frac{1}{L} \left| \int_0^L T_v(t) e^{-i r t} \, dt \right|^2.$$

Because $T_v(t) = v_0 + \sqrt{2}\sum_{m=1}^N v_m \cos(2\pi m t / L)$ satisfies the reflection symmetry $T_v(L - t) = T_v(t)$ on $[0, L]$, the Fourier power spectrum is identically the cosine transform of the spatial autocorrelation function of $T_v$. In terms of the causal Volterra convolution $K_v^{\mathrm{phys}}(x) = 2\int_0^x T_v(t) T_v(x - t) dt$, the autocorrelation corresponds to the reflected kernel $\frac{1}{L} K_v^{\mathrm{phys}}(L - y)$:

$$K_{\mathrm{Fourier}}(v, r, L) = \frac{1}{L} \int_0^L K_v^{\mathrm{phys}}(L - y) \cos(r y) \, dy.$$

By direct integration of the Volterra convolution against $\cos(r L \omega)$ on $[0, 1]$, the boundary terms at the endpoint $x = L$ factor out cleanly, isolating the common oscillatory factor $1 - \cos(rL)$:

$$K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r),$$

where $R_v(r)$ is a rational function of $r^2$. Defining the fundamental lattice scale:

$$\kappa = \frac{2\pi}{L}, \qquad a_m = \kappa m = \frac{2\pi m}{L},$$

and the elementary mode integrals:

$$S_{\mathrm{bar}}(m) = \frac{a_m}{a_m^2 - r^2}, \qquad C_{\mathrm{bar}}(m) = \frac{r^2 + a_m^2}{L (r^2 - a_m^2)^2},$$

the rational kernel was originally obtained by expanding the four-part interaction sum:

$$\begin{aligned}
R_v(r) &= \frac{2 v_0^2}{L r^2} + 2 \sum_{m=1}^{N} v_m^2 C_{\mathrm{bar}}(m) - \frac{2\sqrt{2} v_0}{\pi} \sum_{m=1}^{N} \frac{v_m S_{\mathrm{bar}}(m)}{m} - \frac{1}{\pi} \sum_{m=1}^{N} \frac{v_m^2 S_{\mathrm{bar}}(m)}{m} \\
&\quad + \frac{4}{\pi} \sum_{1 \le m < n \le N} v_m v_n \frac{m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n)}{n^2 - m^2}.
\end{aligned}$$

---

## 3. The Exact Rational Resolvent Identity

Using the analytic expression for the reduced kernel $R_v(r)$ derived from the Archimedean Volterra integral, we state the first central theorem of this paper. This result proves **algebraically, in closed form, and completely independently of numerical quadrature**, that the four-term interaction sum matches the square of a single Cauchy resolvent at every finite dimension $N$.

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

and verify equivalence with the four distinct interaction blocks of $R_v(r)$ term-by-term:

**1. The $v_0^2$ Block:**
The leading term is $\frac{2}{L} (\frac{v_0}{r})^2 = \frac{2 v_0^2}{L r^2}$, which matches the first term of $R_v(r)$ identically.

**2. The $v_0 v_m$ Cross-Term Block:**
In the square expansion, the cross-term between $v_0$ and $v_m$ is:
$$\frac{2}{L} \left[ 2 \frac{v_0}{r} \cdot \sqrt{2} \frac{r v_m}{r^2 - a_m^2} \right] = \frac{4\sqrt{2} v_0 v_m}{L (r^2 - a_m^2)}.$$
In $R_v(r)$, this term is given by $-\frac{2\sqrt{2}}{\pi} v_0 v_m \frac{S_{\mathrm{bar}}(m)}{m}$. Recalling $S_{\mathrm{bar}}(m) = \frac{a_m}{a_m^2 - r^2}$ and $a_m = \frac{2\pi m}{L}$:
$$-\frac{2\sqrt{2}}{\pi m} S_{\mathrm{bar}}(m) = -\frac{2\sqrt{2}}{\pi m} \frac{2\pi m / L}{a_m^2 - r^2} = \frac{4\sqrt{2}}{L (r^2 - a_m^2)}.$$
Multiplying by $v_0 v_m$ shows that this block matches identically.

**3. The $v_m^2$ Diagonal Block:**
In the square expansion, the diagonal term for mode $m$ is:
$$\frac{2}{L} \cdot 2 \cdot \left(\frac{r v_m}{r^2 - a_m^2}\right)^2 = \frac{4 v_m^2}{L} \frac{r^2}{(r^2 - a_m^2)^2}.$$
In $R_v(r)$, the coefficient of $v_m^2$ is:
$$2 C_{\mathrm{bar}}(m) - \frac{1}{\pi m} S_{\mathrm{bar}}(m) = \frac{2(r^2 + a_m^2)}{L (r^2 - a_m^2)^2} - \frac{1}{\pi m} \frac{a_m}{a_m^2 - r^2}.$$
Using $a_m / (\pi m) = 2/L$, the second term becomes $-\frac{2}{L(a_m^2 - r^2)} = \frac{2(r^2 - a_m^2)}{L (r^2 - a_m^2)^2}$. Adding:
$$\frac{2(r^2 + a_m^2) + 2(r^2 - a_m^2)}{L (r^2 - a_m^2)^2} = \frac{4 r^2}{L (r^2 - a_m^2)^2}.$$
This matches the diagonal term of the square identically.

**4. The $v_m v_n$ ($m < n$) Off-Diagonal Block:**
In the square expansion, the pairwise cross-term between mode $m$ and mode $n$ ($m \ne n$) is:
$$\frac{2}{L} \cdot 2 \cdot \sqrt{2} \frac{r v_m}{r^2 - a_m^2} \cdot \sqrt{2} \frac{r v_n}{r^2 - a_n^2} = \frac{8 v_m v_n}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)}.$$
In $R_v(r)$, the corresponding term is $\frac{4 v_m v_n}{\pi} \frac{m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n)}{n^2 - m^2}$. Partial fractions on the numerator using $S_{\mathrm{bar}}(m) = \frac{\kappa m}{\kappa^2 m^2 - r^2}$:
$$m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n) = \frac{\kappa m^2}{\kappa^2 m^2 - r^2} - \frac{\kappa n^2}{\kappa^2 n^2 - r^2} = \kappa \left[ \frac{n^2}{r^2 - \kappa^2 n^2} - \frac{m^2}{r^2 - \kappa^2 m^2} \right].$$
Combining over a common denominator:
$$n^2 (r^2 - \kappa^2 m^2) - m^2 (r^2 - \kappa^2 n^2) = r^2 (n^2 - m^2).$$
Therefore:
$$\frac{m S_{\mathrm{bar}}(m) - n S_{\mathrm{bar}}(n)}{n^2 - m^2} = \frac{\kappa r^2 (n^2 - m^2)}{(n^2 - m^2)(r^2 - a_m^2)(r^2 - a_n^2)} = \frac{2\pi}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)}.$$
Multiplying by $\frac{4 v_m v_n}{\pi}$:
$$\frac{4 v_m v_n}{\pi} \cdot \frac{2\pi}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)} = \frac{8 v_m v_n}{L} \frac{r^2}{(r^2 - a_m^2)(r^2 - a_n^2)}.$$
This matches the cross-terms of the square identically. Since all four blocks match identically, the algebraic identity $R_v(r) \equiv \frac{1}{r^2} A(1/r^2)$ is exact. Although the bracketed representation inside the square contains the odd factor $1/r$, every cross-term in its square expansion is an even function of $r$, so that $R_v(r)$ is identically an even rational function of $r$, hence a rational function of $r^2$. $\blacksquare$

### Theorem 3.2 (Neumann Resolvent Representation and Heat-Kernel Identity)
*Let $\mathcal{L} = -\frac{d^2}{dt^2}$ denote the Neumann Laplacian on the physical interval $[0, L]$ with domain $\mathcal{D}(\mathcal{L}) = \{f \in H^2(0, L) : f'(0) = f'(L) = 0\}$. The normalized cosine eigenbasis of $\mathcal{L}$ is:*

$$\phi_0(t) = 1, \qquad \phi_m(t) = \sqrt{2} \cos\left(\frac{2\pi m t}{L}\right) \quad (m \ge 1),$$

*with eigenvalues $\mathcal{L} \phi_m = a_m^2 \phi_m$, where $a_m = \frac{2\pi m}{L}$. For any canonical coefficient vector $v \in \mathbb{R}^{N+1}$, the spatial trigonometric wave $T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos(a_m t)$ expands in this basis as $T_v(t) = \sum_{m=0}^N v_m \phi_m(t)$.*

*Then:*
1. **Boundary Resolvent Identity:**
   $$D(z) \equiv \big[(I + z\mathcal{L})^{-1} T_v\big](0) \qquad \forall z \in \mathbb{C} \setminus \{-a_1^{-2}, \dots, -a_N^{-2}\}.$$
2. **Heat-Resolvent Integral Representation:**
   *For $\operatorname{Re}(z) > 0$:*
   $$D(z) = \int_0^\infty e^{-s} \big[ e^{-s z \mathcal{L}} T_v \big](0) \, ds.$$
3. **Boundary Resolvent Expansion:**
   *The Taylor coefficients of $D(z)$ around $z = 0$ reproduce the endpoint-jet hierarchy:*
   $$D(z) = \sum_{k=0}^\infty (-1)^k \big[\mathcal{L}^k T_v\big](0) z^k = \sum_{k=0}^\infty T_v^{(2k)}(0) z^k = \sum_{k=0}^\infty D_k z^k.$$

*Proof.* (1) By spectral decomposition, $(I + z\mathcal{L})^{-1} \phi_m = \frac{1}{1 + a_m^2 z} \phi_m$. Evaluating at $t = 0$ gives $[(I + z\mathcal{L})^{-1} T_v](0) = v_0 \phi_0(0) + \sum_{m=1}^N \frac{v_m}{1 + a_m^2 z} \phi_m(0) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 z} = D(z)$.

(2) For $\operatorname{Re}(z) > 0$, using $(1 + a_m^2 z)^{-1} = \int_0^\infty e^{-s(1 + a_m^2 z)} ds = \int_0^\infty e^{-s} e^{-s z a_m^2} ds$, summing over eigenmodes gives $D(z) = \int_0^\infty e^{-s} [ e^{-s z \mathcal{L}} T_v ](0) ds$.

(3) Expanding $(I + z\mathcal{L})^{-1} = \sum_{k=0}^\infty (-1)^k z^k \mathcal{L}^k$ around $z = 0$ and noting $[\mathcal{L}^k T_v](0) = (-1)^k T_v^{(2k)}(0) = (-1)^k D_k$ yields the Taylor expansion. $\blacksquare$

---

## 4. Unconditional Finite-$N$ Positivity and the Spectral Lattice Identity

### Theorem 4.1 (Unconditional Finite-$N$ Pointwise Kernel Positivity and Entire Amplitude)
*The Fourier-side Archimedean kernel $K_{\mathrm{Fourier}}(v, r, L)$ is unconditionally pointwise non-negative on the real line for all $v \in \mathbb{R}^{N+1}$:*

$$K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0 \qquad \forall r \in \mathbb{R},$$

*where $\Phi_v(r)$ is an entire function of exponential type at most $L/2$ given by:*

$$\Phi_v(r) = \frac{2}{\sqrt{L}} \left[ v_0 \frac{\sin(rL/2)}{r} + \sqrt{2} \sum_{m=1}^{N} v_m \frac{r \sin(rL/2)}{r^2 - a_m^2} \right].$$

### Proof of Theorem 4.1
Using $1 - \cos(rL) = 2 \sin^2(rL/2)$, we factor the full kernel:

$$K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r) = 2 \sin^2(rL/2) \cdot \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} \right]^2.$$

Bringing $\sqrt{2}\sin(rL/2)$ inside the square yields $\Phi_v(r)^2$. 

At $r = 0$, $\lim_{r\to 0} \frac{\sin(rL/2)}{r} = \frac{L}{2}$ and $\lim_{r\to 0} \frac{r\sin(rL/2)}{r^2 - a_m^2} = 0$, giving $\Phi_v(0) = \sqrt{L} v_0$.

At $r = \pm a_m$, Taylor expansion of $\sin(rL/2)$ around $r = a_m$ gives $\sin(rL/2) = (-1)^m \frac{L}{2}(r - a_m) + \mathcal{O}((r - a_m)^3)$. Because $r^2 - a_m^2 = (r - a_m)(r + a_m)$, both apparent poles at $r = \pm a_m$ are removable:

$$\lim_{r\to \pm a_m} \frac{r \sin(rL/2)}{r^2 - a_m^2} = (-1)^m \frac{L}{4}.$$

Because all singularities are removable, $\Phi_v(r)$ is entire on $\mathbb{C}$ of exponential type at most $L/2$. Because $v$ is real, $\Phi_v(r) \in \mathbb{R}$ for all $r \in \mathbb{R}$, which forces $\Phi_v(r)^2 \ge 0$ unconditionally. $\blacksquare$

*Remark 4.1.1 (Distinction Between Kernel Positivity and Functional Positivity).*
The unconditional pointwise non-negativity $K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0$ on the real line guarantees that the Fourier spectral density is non-negative everywhere. However, this does *not* imply that the integrated Archimedean functional $\mathcal{Q}_{\mathrm{arch}}(v) = \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$ is non-negative, because the Archimedean weight $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$ is not positive everywhere (it is negative for $r < r_0 \approx 9.77$ and positive for $r > r_0$). Consequently, kernel non-negativity $K_{\mathrm{Fourier}} \ge 0$ does not imply positivity of the full Weil quadratic form $\mathcal{Q}(v) \ge 0$, nor does it imply the Riemann Hypothesis. It provides an unconditionally positive spectral density entering the explicit formula.

*Remark 4.1.2 (Three-Level Analytic Hierarchy: Meromorphic $R_v$, Entire $K_{\mathrm{Fourier}}$, and Entire $\Phi_v$).*
The analytic structure of the Fourier Archimedean kernel forms a precise three-level hierarchy:
1. **The Reduced Resolvent $R_v(r)$:** A meromorphic rational function on $\mathbb{C}$ possessing second-order poles at $r = 0$ and $r = \pm a_m$ ($m = 1, \dots, N$).
2. **The Archimedean Kernel $K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r)$:** The second-order zeros of the trigonometric envelope $1 - \cos(rL) = 2\sin^2(rL/2)$ at $r = 0$ and at each discrete lattice node $r = \pm a_m = \pm 2\pi m/L$ exactly cancel each second-order pole of $R_v(r)$ via removable singularities, extending $K_{\mathrm{Fourier}}$ to an entire function on the full complex plane $\mathbb{C}$.
3. **The Entire Amplitude $\Phi_v(r)$:** Factoring $K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2$, the physical Fourier amplitude $\Phi_v(r) = \frac{1}{\sqrt{L}} \int_0^L T_v(t) e^{-irt} dt$ is itself an entire function of exponential type at most $L/2$ by the Paley–Wiener theorem, ensuring that $K_{\mathrm{Fourier}} \ge 0$ unconditionally on $\mathbb{R}$.

### Corollary 4.2 (Universal Fourier Factorization)
*The entire Fourier amplitude $\Phi_v(r)$ factors directly into the product of the universal sinc envelope and the boundary Neumann resolvent evaluated at $z = -1/r^2$:*

$$\Phi_v(r) \equiv \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D\left(-\frac{1}{r^2}\right), \qquad R_v(r) \equiv \frac{2}{L r^2} \left[ D\left(-\frac{1}{r^2}\right) \right]^2.$$

*Proof.* Factoring $1/r$ inside $\Phi_v(r)$:
$$\frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} = \frac{1}{r} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - a_m^2 / r^2} \right] = \frac{1}{r} D\left(-\frac{1}{r^2}\right).$$
Substituting into $\Phi_v(r)$ and squaring yields the result. $\blacksquare$

### Theorem 4.3 (Spectral Lattice Sampling Identity)
*At the discrete Fourier frequencies $a_m = 2\pi m / L$, the Archimedean kernel samples the squared Fourier coefficients orthogonally:*

$$K_{\mathrm{Fourier}}(v, 0, L) = L v_0^2 = L u_0^2, \qquad K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 = L u_m^2 \quad (m = 1, \dots, N).$$

*Proof.* At $r = 0$, $\Phi_v(0) = \sqrt{L} v_0$, so $K_{\mathrm{Fourier}}(v, 0, L) = L v_0^2$. For $m \in \{1, \dots, N\}$, $\sin(a_m L/2) = \sin(\pi m) = 0$, so all terms vanish except the removable diagonal limit $\lim_{r\to a_m} \sqrt{2} v_m \frac{r \sin(rL/2)}{r^2 - a_m^2} = (-1)^m \sqrt{2} v_m \frac{L}{4}$. Thus $\Phi_v(a_m) = \frac{2}{\sqrt{L}} (-1)^m \sqrt{2} v_m \frac{L}{4} = (-1)^m \sqrt{\frac{L}{2}} v_m$. Squaring gives $K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 = L u_m^2$. $\blacksquare$

---

## 5. Exact Archimedean Cauchy Transform and Closed-Form Pole Decomposition

Evaluating the integrated Archimedean quadratic form $\mathcal{Q}_{\mathrm{arch}}(v) = \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$ against the digamma weight $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$ has historically required numerical quadrature with truncation cutoffs. The rational resolvent factorization $K_{\mathrm{Fourier}}(r) = (1 - \cos(rL)) R_v(r)$ permits an exact closed-form evaluation of its Cauchy transform and an unconditionally convergent algebraic pole decomposition.

### Theorem 5.1 (Exact Archimedean Cauchy Transform Identity)
*Let $c > 1$, $L = \log c$, $N \ge 1$, and $a_m = \frac{2\pi m}{L}$. For any canonical coefficient vector $v \in \mathbb{R}^{N+1}$ and any pole parameter $q > 0$, the Cauchy transform of the Fourier Archimedean kernel evaluates in exact closed algebraic form:*

$$J(q) \equiv \frac{1}{\pi} \int_0^\infty \frac{2 q}{q^2 + r^2} K_{\mathrm{Fourier}}(v, r, L) \, dr = \frac{2 v_0^2}{q} + \sum_{m=1}^N \frac{2 q v_m^2}{q^2 + a_m^2} - \frac{2(1 - e^{-q L})}{L q^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2} \right]^2.$$

### Primary Proof via Reflected Spatial Autocorrelation and Poisson Kernel
By definition of the Fourier amplitude $\Phi_v(r) = \frac{1}{\sqrt{L}} \int_0^L T_v(t) e^{-i r t} \, dt$, the Fourier Archimedean kernel is the power spectrum:
$$K_{\mathrm{Fourier}}(v, r, L) = \frac{1}{L} \left| \int_0^L T_v(t) e^{-i r t} \, dt \right|^2 = \frac{1}{L} \int_0^L \int_0^L T_v(t_1) T_v(t_2) \cos(r(t_1 - t_2)) \, dt_1 dt_2.$$
Substituting this double integral into the Cauchy transform $J(q) = \frac{1}{\pi} \int_0^\infty \frac{2q}{q^2 + r^2} K_{\mathrm{Fourier}}(v, r, L) \, dr$ and interchanging integration via Fubini's theorem (justified since $T_v$ is smooth and bounded on $[0, L]$):
$$J(q) = \frac{1}{L} \int_0^L \int_0^L T_v(t_1) T_v(t_2) \left[ \frac{2q}{\pi} \int_0^\infty \frac{\cos(r(t_1 - t_2))}{q^2 + r^2} \, dr \right] dt_1 dt_2.$$
The inner integral is the classical Poisson kernel:
$$\frac{2q}{\pi} \int_0^\infty \frac{\cos(r y)}{q^2 + r^2} \, dr = e^{-q |y|} \qquad (q > 0),$$
which yields the exact spatial representation:
$$J(q) = \frac{1}{L} \int_0^L \int_0^L T_v(t_1) T_v(t_2) e^{-q |t_1 - t_2|} \, dt_1 dt_2.$$

Expanding $T_v(t) = \sum_{m=0}^N c_m \cos(a_m t)$ where $c_0 = v_0$, $c_m = \sqrt{2} v_m$ ($m \ge 1$), and $a_m = \frac{2\pi m}{L}$ (with $a_0 = 0$):
$$J(q) = \frac{1}{L} \sum_{m,n=0}^N c_m c_n I_{mn}(q), \qquad I_{mn}(q) \equiv \int_0^L \int_0^L \cos(a_m t_1) \cos(a_n t_2) e^{-q |t_1 - t_2|} \, dt_1 dt_2.$$
By reflection symmetry in $t_1, t_2$, the integral decomposes over the triangle $t_2 \le t_1$:
$$I_{mn}(q) = \int_0^L dt_1 \cos(a_m t_1) \int_0^{t_1} dt_2 \cos(a_n t_2) e^{-q(t_1 - t_2)} + \int_0^L dt_2 \cos(a_n t_2) \int_0^{t_2} dt_1 \cos(a_m t_1) e^{-q(t_2 - t_1)}.$$
Evaluating the inner elementary integral:
$$\int_0^{t_1} \cos(a_n t_2) e^{q t_2} \, dt_2 = \frac{q e^{q t_1} \cos(a_n t_1) + a_n e^{q t_1} \sin(a_n t_1) - q}{q^2 + a_n^2}.$$
Multiplying by $e^{-q t_1} \cos(a_m t_1)$ and integrating $t_1 \in [0, L]$:
1. **Orthogonality Term:**
   $$\frac{q}{q^2 + a_n^2} \int_0^L \cos(a_m t_1) \cos(a_n t_1) \, dt_1 = \frac{q}{q^2 + a_n^2} \frac{L}{2} (1 + \delta_{m0}) \delta_{mn}.$$
2. **Sine Cross-Term:**
   $$\frac{a_n}{q^2 + a_n^2} \int_0^L \cos(a_m t_1) \sin(a_n t_1) \, dt_1 = 0 \qquad \forall m, n \ge 0,$$
   since $a_k = 2\pi k / L$ implies $\sin((a_n \pm a_m)t_1)$ has integer periods over $[0, L]$.
3. **Boundary Exponential Term:**
   $$-\frac{q}{q^2 + a_n^2} \int_0^L e^{-q t_1} \cos(a_m t_1) \, dt_1 = -\frac{q}{q^2 + a_n^2} \operatorname{Re}\left( \frac{1 - e^{-(q - i a_m)L}}{q - i a_m} \right) = -\frac{q^2 (1 - e^{-q L})}{(q^2 + a_n^2)(q^2 + a_m^2)},$$
   using $e^{i a_m L} = e^{i 2\pi m} = 1$.

Adding the symmetric transpose ($m \leftrightarrow n$), the diagonal term doubles and the boundary term produces:
$$I_{mn}(q) = \frac{q L (1 + \delta_{m0})}{q^2 + a_m^2} \delta_{mn} - \frac{2 q^2 (1 - e^{-q L})}{(q^2 + a_m^2)(q^2 + a_n^2)}.$$
Summing against $\frac{1}{L} c_m c_n$:
- For the diagonal sum:
  $$\frac{1}{L} \sum_{m=0}^N c_m^2 \frac{q L (1 + \delta_{m0})}{q^2 + a_m^2} = \frac{2 v_0^2}{q} + \sum_{m=1}^N \frac{2 q v_m^2}{q^2 + a_m^2}.$$
- For the boundary sum:
  $$-\frac{2(1 - e^{-q L})}{L} \sum_{m,n=0}^N c_m c_n \frac{q^2}{(q^2 + a_m^2)(q^2 + a_n^2)} = -\frac{2(1 - e^{-q L})}{L q^2} \left[ \sum_{m=0}^N c_m \frac{q^2}{q^2 + a_m^2} \right]^2 = -\frac{2(1 - e^{-q L})}{L q^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2} \right]^2.$$
This completes the direct, elementary derivation. $\blacksquare$

### Alternative Proof via Complex Frequency Contour Integration
Alternatively, Theorem 5.1 can be derived by contour integration in the complex frequency plane. Decompose the kernel into $K_{\mathrm{Fourier}}(v, r, L) = f_1(r) + f_2(r)$ where $f_1(z) = H(z)(1 - e^{i z L})$, $f_2(z) = H(z)(1 - e^{-i z L})$, and $H(z) = \frac{1}{2} R_v(z) = \frac{1}{L z^2} [D(-1/z^2)]^2$. Extending over $(-\infty, \infty)$:

$$J(q) = \frac{q}{\pi} \lim_{\epsilon \to 0^+} \left[ \int_{C_\epsilon^+} \frac{f_1(z)}{z^2 + q^2} \, dz + \int_{C_\epsilon^+} \frac{f_2(z)}{z^2 + q^2} \, dz \right],$$

where $C_\epsilon^+$ is the real axis indented into the upper half-plane $\mathbb{C}^+$ around $z = 0$ and $z = \pm a_m$.

1. **Upper Contour for $f_1(z)$:** In $\mathbb{C}^+$, $|e^{izL}| = e^{-L \operatorname{Im}(z)} < 1$, so $f_1(z) = \mathcal{O}(|z|^{-2})$. Closing with a large semicircle in $\mathbb{C}^+$ encloses only $z = iq$. By Cauchy's residue theorem:
   $$\int_{C_\epsilon^+} \frac{f_1(z)}{z^2 + q^2} \, dz = 2\pi i \operatorname{Res}_{z = i q} \left[ \frac{f_1(z)}{z^2 + q^2} \right] = \frac{\pi}{q} H(i q) (1 - e^{-q L}).$$
   Evaluating $H(iq) = -\frac{1}{L q^2} [ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2} ]^2$.

2. **Lower Contour for $f_2(z)$:** In $\mathbb{C}^-$, $|e^{-izL}| = e^{L \operatorname{Im}(z)} < 1$. Closing with a clockwise semicircle in $\mathbb{C}^-$ encloses $z = -iq$, the origin $z = 0$, and the discrete lattice nodes $z = \pm a_m$.
   - At $z = -iq$: $-2\pi i \operatorname{Res}_{z = -iq} = \frac{\pi}{q} H(iq) (1 - e^{-qL})$, matching $f_1$. Summing with $f_1$ gives the third term of $J(q)$:
     $$\frac{q}{\pi} \left[ 2 \frac{\pi}{q} H(iq) (1 - e^{-qL}) \right] = -\frac{2(1 - e^{-q L})}{L q^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2} \right]^2.$$
   - At $z = 0$: Since $H(z) = \frac{v_0^2}{L z^2} + \mathcal{O}(1)$ and $1 - e^{-izL} = izL + \mathcal{O}(z^2)$, $f_2(z) = \frac{iv_0^2}{z} + \mathcal{O}(1)$. Residue at $z = 0$ gives $-2\pi i \frac{iv_0^2}{q^2} = \frac{2\pi v_0^2}{q^2}$. Multiplying by $\frac{q}{\pi}$ gives the leading term:
     $$\frac{2 v_0^2}{q}.$$
   - At $z = \pm a_m$: The double pole of $H(z)$ is regularized by the zero of $1 - e^{-izL}$, yielding residue $\operatorname{Res}_{z = \pm a_m} f_2(z) = \frac{i}{2} v_m^2$. Evaluating in the integrand gives $-2\pi i \frac{i v_m^2}{q^2 + a_m^2} = \frac{2\pi v_m^2}{q^2 + a_m^2}$. Multiplying by $\frac{q}{\pi}$ gives:
     $$\frac{2 q v_m^2}{q^2 + a_m^2}.$$

Summing all contributions matches the result of the primary proof identically. $\blacksquare$

### Theorem 5.2 (Spatial Laplace Duality via the Reflected Autocorrelation Kernel)
*The Archimedean Cauchy transform $J(q)$ is identically the spatial Laplace transform of the reflected physical kernel:*

$$J(q) \equiv \frac{1}{L} \int_0^L K_v^{\mathrm{phys}}(L - y) e^{-q y} \, dy \qquad \forall q > 0.$$

*Proof.* In the spatial representation of $J(q)$ derived in the primary proof of Theorem 5.1:
$$J(q) = \frac{1}{L} \int_0^L \int_0^L T_v(t_1) T_v(t_2) e^{-q|t_1 - t_2|} \, dt_1 dt_2.$$
Splitting into $t_1 \ge t_2$ and setting $y = t_1 - t_2 \in [0, L]$ (with $t_2 = t \in [0, L-y]$), symmetry gives:
$$J(q) = \frac{2}{L} \int_0^L e^{-qy} \left( \int_0^{L-y} T_v(t) T_v(t + y) \, dt \right) dy.$$
By the reflection symmetry $T_v(L - s) = T_v(s)$ on $[0, L]$, setting $s = t + y$ yields $T_v(t + y) = T_v(L - y - t)$. Recalling the physical Volterra convolution $K_v^{\mathrm{phys}}(x) = 2 \int_0^x T_v(t) T_v(x - t) dt$ (Section 2.1), evaluating at $x = L - y$ gives:
$$K_v^{\mathrm{phys}}(L - y) = 2 \int_0^{L-y} T_v(t) T_v(L - y - t) \, dt = 2 \int_0^{L-y} T_v(t) T_v(t + y) \, dt.$$
Substituting this into the double integral yields:
$$J(q) = \frac{1}{L} \int_0^L K_v^{\mathrm{phys}}(L - y) e^{-q y} \, dy.$$

*Verification on Constant Wave:* Taking $v = (1, 0, \dots, 0)$, we have $T_v(t) = 1$, whence $K_v^{\mathrm{phys}}(x) = 2 \int_0^x 1 \cdot 1 \, dt = 2x$, and $K_v^{\mathrm{phys}}(L - y) = 2(L - y)$. Evaluating the Laplace transform:
$$\frac{1}{L} \int_0^L 2(L - y) e^{-qy} \, dy = \frac{2}{q} - \frac{2(1 - e^{-qL})}{L q^2},$$
which matches the closed-form evaluation from Theorem 5.1 identically. $\blacksquare$

### Corollary 5.3 (Exact Closed-Form Pole Decomposition of $\mathcal{Q}_{\mathrm{arch}}(v)$)
*Using the Weierstrass partial fraction expansion for the Digamma function:*

$$h_+(r) = \operatorname{Re}\psi\left(\frac{1}{4} + \frac{ir}{2}\right) - \log \pi = C_{\mathrm{arch}} + \sum_{n=0}^\infty \left[ \frac{1}{n+1} - \frac{2 q_n}{q_n^2 + r^2} \right], \qquad q_n = 2n + \frac{1}{2},$$

*where $C_{\mathrm{arch}} = -\gamma - \log \pi \approx -1.7219455$, the continuous Archimedean quadratic form decomposes into an exact, unconditionally convergent algebraic series without numerical quadrature:*

$$\mathcal{Q}_{\mathrm{arch}}(v) \equiv C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty \left[ \frac{\|v\|_2^2}{n+1} - J(q_n) \right],$$

*where each term $J(q_n)$ is evaluated in closed algebraic form via Theorem 5.1.*

*Proof.* 
1. **Term-by-Term Finiteness:**
   For each $n \ge 0$, define the rational kernel:
   $$w_n(r) \equiv \frac{1}{n+1} - \frac{2 q_n}{q_n^2 + r^2} = \frac{r^2 - 4n^2 + 1/4}{(n+1)((2n + 1/2)^2 + r^2)}.$$
   Since $K_{\mathrm{Fourier}}(v, r, L) \ge 0$ is continuous on $[0, \infty)$ and decays as $\mathcal{O}(r^{-2})$ as $r \to \infty$ (Theorem 3.1), and $w_n(r)$ is bounded on $[0, \infty)$, each individual integral:
   $$I_n \equiv \frac{1}{\pi} \int_0^\infty w_n(r) K_{\mathrm{Fourier}}(v, r, L) \, dr = \frac{\|v\|_2^2}{n+1} - J(q_n)$$
   is absolutely convergent and finite.

2. **Rigorous Finite-$M$ Remainder Bound and Summation-Integration Interchange:**
   To establish the identity without relying on heuristic dominated-convergence envelopes, consider the partial sum:
   $$h_M(r) \equiv C_{\mathrm{arch}} + \sum_{n=0}^M w_n(r), \qquad w_n(r) = \frac{1}{n+1} - \frac{2q_n}{q_n^2 + r^2}.$$
   The remainder function $R_M(r) \equiv h_+(r) - h_M(r) = \sum_{n=M+1}^\infty w_n(r)$ decomposes algebraically as:
   $$w_n(r) = \left( \frac{1}{n+1} - \frac{2}{q_n} \right) + \frac{2r^2}{q_n(q_n^2 + r^2)} = -\frac{3/4}{(n+1)(n + 1/4)} + \frac{2r^2}{q_n(q_n^2 + r^2)}.$$
   Summing over $n \ge M+1$, the negative part is an $r$-independent convergent numerical tail bounded by:
   $$\sum_{n=M+1}^\infty \frac{3/4}{(n+1)(n + 1/4)} \le \frac{3}{4} \sum_{n=M+1}^\infty \frac{1}{n^2} \le \frac{3}{4M}.$$
   For the positive part, since $K_{\mathrm{Fourier}}(v, r, L) \le \frac{C_K}{1 + r^2}$ on $[0, \infty)$ (Theorem 4.1), integrating each term against $K_{\mathrm{Fourier}}(r)$ yields:
   $$\int_0^\infty \frac{2r^2}{q_n(q_n^2 + r^2)} K_{\mathrm{Fourier}}(v, r, L) \, dr \le C_K \int_0^\infty \frac{2r^2}{(1 + r^2) q_n (q_n^2 + r^2)} \, dr = \frac{\pi C_K}{q_n(q_n + 1)} \le \frac{\pi C_K}{4n^2}.$$
   Summing over $n \ge M+1$:
   $$\sum_{n=M+1}^\infty \frac{\pi C_K}{4n^2} \le \frac{\pi C_K}{4M}.$$
   Consequently, the integrated remainder satisfies the uniform finite-$M$ bound:
   $$\left| \mathcal{Q}_{\mathrm{arch}}(v) - \frac{1}{\pi} \int_0^\infty h_M(r) K_{\mathrm{Fourier}}(v, r, L) \, dr \right| \le \frac{1}{\pi} \int_0^\infty |R_M(r)| K_{\mathrm{Fourier}}(v, r, L) \, dr \le \frac{C_{\mathrm{rem}}}{M},$$
   where $C_{\mathrm{rem}} = \frac{3}{4\pi} \|K_{\mathrm{Fourier}}\|_{L^1} + \frac{C_K}{4} < \infty$.
   Taking $M \to \infty$ proves unconditionally that:
   $$\mathcal{Q}_{\mathrm{arch}}(v) = \lim_{M\to\infty} \frac{1}{\pi} \int_0^\infty h_M(r) K_{\mathrm{Fourier}}(v, r, L) \, dr = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty \left[ \frac{\|v\|_2^2}{n+1} - J(q_n) \right].$$

3. **Large-$q$ Remainder and Absolute Convergence:**
   To establish the convergence rate, from Theorem 5.1 we expand $J(q)$ for large $q > 0$:
   $$J(q) = \frac{2v_0^2}{q} + \sum_{m=1}^N \frac{2qv_m^2}{q^2 + a_m^2} - \frac{2(1 - e^{-qL})}{L q^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2} \right]^2.$$
   Using the exact algebraic identity $\frac{2qv_m^2}{q^2 + a_m^2} = \frac{2v_m^2}{q} - \frac{2 a_m^2 v_m^2}{q(q^2 + a_m^2)}$ and $\|v\|_2^2 = v_0^2 + \sum_{m=1}^N v_m^2$, we have:
   $$\frac{2\|v\|_2^2}{q} - J(q) = \sum_{m=1}^N \frac{2 a_m^2 v_m^2}{q(q^2 + a_m^2)} + \frac{2(1 - e^{-qL})}{L q^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2} \right]^2.$$
   For each $m$, $0 \le \frac{2 a_m^2 v_m^2}{q(q^2 + a_m^2)} \le \frac{2 a_m^2 v_m^2}{q^3}$. For the boundary term, $|v_0 + \sqrt{2}\sum_{m=1}^N \frac{q^2 v_m}{q^2 + a_m^2}| \le \sqrt{2N+1}\|v\|_2$. Thus:
   $$\left| J(q) - \frac{2\|v\|_2^2}{q} \right| \le \frac{C_J}{q^2}, \qquad C_J \equiv \frac{2(2N+1)\|v\|_2^2}{L} + 2 \sum_{m=1}^N a_m^2 v_m^2.$$
   Evaluating at $q_n = 2n + 1/2 = 2(n + 1/4)$, we obtain $\frac{2\|v\|_2^2}{q_n} = \frac{\|v\|_2^2}{n + 1/4}$, so the series summand satisfies:
   $$\left| \frac{\|v\|_2^2}{n+1} - J(q_n) \right| \le \|v\|_2^2 \left( \frac{1}{n + 1/4} - \frac{1}{n+1} \right) + \frac{C_J}{(2n + 1/2)^2} = \frac{\frac{3}{4} \|v\|_2^2}{(n+1)(n + 1/4)} + \frac{C_J}{(2n + 1/2)^2} \le \frac{C_{\mathrm{tail}}}{n^2}.$$
   Since $\sum_{n=1}^\infty n^{-2} = \frac{\pi^2}{6} < \infty$, the series converges absolutely and unconditionally with fast $\mathcal{O}(n^{-2})$ convergence. $\blacksquare$

### Corollary 5.4 (Exact Closed-Form Digamma Identity for $\mathcal{Q}_{\mathrm{arch}}(v)$)
*By summing the harmonic difference and the discrete lattice pole series in Corollary 5.3 analytically in terms of the digamma function $\psi(z)$, the Archimedean quadratic form decomposes into the exact discrete lattice sum plus an exponentially suppressed boundary leakage term:*

$$\mathcal{Q}_{\mathrm{arch}}(v) = \mathcal{Q}_{\mathrm{arch}}^{\mathrm{lattice}}(v) + \mathcal{Q}_{\mathrm{arch}}^{\mathrm{leakage}}(v),$$

*where the discrete lattice contribution is:*
$$\mathcal{Q}_{\mathrm{arch}}^{\mathrm{lattice}}(v) \equiv h_+(0) v_0^2 + \sum_{m=1}^N v_m^2 h_+(a_m),$$
*and the boundary leakage is:*
$$\mathcal{Q}_{\mathrm{arch}}^{\mathrm{leakage}}(v) \equiv \sum_{n=0}^\infty \frac{2(1 - e^{-q_n L})}{L q_n^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q_n^2 v_m}{q_n^2 + a_m^2} \right]^2 = \sum_{n=0}^\infty \frac{2(1 - e^{-q_n L})}{L q_n^2} \left[ D\left(\frac{1}{q_n^2}\right) \right]^2 \ge 0,$$
*with $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$ and $q_n = 2n + 1/2$.*

*Proof.* In the summand of Corollary 5.3, write $\frac{\|v\|_2^2}{n+1} - J(q_n) = (\frac{\|v\|_2^2}{n+1} - \frac{2\|v\|_2^2}{q_n}) + (\frac{2\|v\|_2^2}{q_n} - J(q_n))$. Summing the first part gives $\|v\|_2^2 \sum_{n=0}^\infty (\frac{1}{n+1} - \frac{1}{n+1/4}) = \|v\|_2^2 (\psi(1/4) + \gamma)$, which combined with $C_{\mathrm{arch}}\|v\|_2^2 = (-\gamma - \log\pi)\|v\|_2^2$ yields $(\psi(1/4) - \log\pi)\|v\|_2^2 = h_+(0)\|v\|_2^2$. For the second part, using the identity $\frac{2\|v\|_2^2}{q_n} - J(q_n) = \sum_{m=1}^N \frac{2 a_m^2 v_m^2}{q_n(q_n^2 + a_m^2)} + \text{boundary term}$, summing over $n$ using $\sum_{n=0}^\infty \frac{2 a_m^2}{q_n(q_n^2 + a_m^2)} = \operatorname{Re}\psi(1/4 + i a_m/2) - \psi(1/4) = h_+(a_m) - h_+(0)$ yields $\sum_{m=1}^N v_m^2 [h_+(a_m) - h_+(0)]$. Adding to $h_+(0)\|v\|_2^2$ cancels the $-h_+(0)\sum v_m^2$ term, leaving $h_+(0) v_0^2 + \sum_{m=1}^N v_m^2 h_+(a_m)$ plus the boundary sum. $\blacksquare$

### Corollary 5.4.1 (Exact Lattice Lower Bound)
*Because $L > 0$, $q_n > 0$, and the boundary summand involves the square $[D(1/q_n^2)]^2$, every term in the boundary leakage series $\mathcal{Q}_{\mathrm{arch}}^{\mathrm{leakage}}(v)$ is non-negative:*
$$\frac{2(1 - e^{-q_n L})}{L q_n^2} \left[ D\left(\frac{1}{q_n^2}\right) \right]^2 \ge 0 \qquad \forall n \ge 0.$$
*Consequently, the discrete lattice sum provides an unconditional lower bound on the Archimedean quadratic form for every canonical vector $v$ and all finite $N$:*

$$\mathcal{Q}_{\mathrm{arch}}(v) \ge h_+(0) v_0^2 + \sum_{m=1}^N v_m^2 h_+(a_m).$$

### Theorem 5.5 (Finite-$T$ Archimedean Cutoff Defect and Endpoint-Jet Representation)
*Let $Q_{\mathrm{arch}}^{(T)}$ be the finite-cutoff Galerkin Archimedean matrix of band $N$ defined via divided differences of the cutoff test function $\psi_{\mathrm{arch}}^{(T)}(x) = \frac{1}{2\pi^2}\int_{-T}^T h_+(r) \operatorname{Re}\hat{S}_x(r) \, dr$ for $T > 0$.*

1. **Exact Divided-Difference Integral Identity:**
   *For any canonical coefficient vector $v \in \mathbb{R}^{N+1}$ corresponding to the even full vector $u \in \mathbb{R}^{2N+1}$, the matrix quadratic form is identically equal to the continuous Fourier integral truncated at $T$:*
   $$v^T Q_{\mathrm{arch}}^{(T)} v \equiv u^T Q_{\mathrm{arch}}^{(T)} u \equiv \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr.$$

2. **Exact Archimedean Cutoff Defect:**
   *Let $Q^{(T)} = Q_{\mathrm{arch}}^{(T)} + Q_{\mathrm{pole}} + Q_{\mathrm{prime}}$ be the full Galerkin matrix with cutoff $T$, and let $u_N$ be an even ground-state eigenvector with eigenvalue $\lambda_N = u_N^T Q^{(T)} u_N$. Let $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{(\infty)}(v_N)$ be the exact continuous tripartite functional evaluated at $T = \infty$ (via Corollary 5.4). Then the discrepancy between the Galerkin eigenvalue and the continuous functional is identically the Archimedean cutoff tail:*
   $$\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) \equiv -\frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr \equiv -\delta_T^{\mathrm{tail}}(v_N).$$
   *In particular, this identity proves that for the fixed finite-dimensional Galerkin vector $v_N$, the entire discrepancy between the finite-cutoff Galerkin eigenvalue $\lambda_N = u_N^T Q^{(T)} u_N$ and the continuous infinite-cutoff functional $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N)$ is identically the continuous Archimedean cutoff tail $\delta_T^{\mathrm{tail}}(v_N)$. (This identity applies strictly to the functional evaluated on $v_N$, and does not assert the absence of projection errors between $v_N$ and any infinite-dimensional continuum minimizer.)*

3. **Taylor Endpoint-Jet Laurent Decomposition:**
   *For any cutoff $T > a_N = \frac{2\pi N}{L}$, the rational resolvent $R_v(r) = \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} \right]^2$ converges uniformly on $[T, \infty)$ as the Laurent series:*
   $$R_v(r) = \sum_{k=0}^\infty \frac{A_k(N)}{r^{2k+2}}, \qquad A_k(N) = \frac{2}{L} (-1)^k \sum_{j=0}^k D_j D_{k-j},$$
   *where $D_0 = T_v(0) = v_0 + \sqrt{2}\sum_{m=1}^N v_m$ and $D_j = T_v^{(2j)}(0) = (-1)^j \sqrt{2}\sum_{m=1}^N a_m^{2j} v_m$ ($j \ge 1$) are the physical endpoint derivatives (Theorem 3.2). Consequently, the tail defect admits the exact series representation:*
   $$\delta_T^{\mathrm{tail}}(v_N) = \sum_{k=0}^\infty A_k(N) \mathcal{J}_k(T, L),$$
   *where the universal tail moments:*
   $$\mathcal{J}_k(T, L) \equiv \frac{1}{\pi} \int_T^\infty \frac{h_+(r) (1 - \cos(rL))}{r^{2k+2}} \, dr$$
   *are independent of the eigenvector $v_N$.*

4. **Explicit Geometric Remainder Bound:**
   *For any integration cutoff $T > a_N = \frac{2\pi N}{L}$ and any truncation order $K \ge 0$, the tail truncation remainder is unconditionally bounded by:*
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{C_K(v, N, L, T)}{T^{2K+3}} \frac{1}{1 - (a_N/T)^2},$$
   *where the order-dependent constant is:*
   $$C_K(v, N, L, T) \equiv \frac{8 (\log T + 1)}{\pi L} \|v\|_1^2 a_N^{2K+2}.$$
   *Equivalently, factoring the powers of $T$ into dimensionless geometric form:*
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{C_0(v, L, T)}{T} \left( \frac{a_N}{T} \right)^{2K+2} \frac{1}{1 - (a_N/T)^2},$$
   *where $C_0(v, L, T) \equiv \frac{8 (\log T + 1)}{\pi L} \|v\|_1^2$ is strictly independent of the truncation order $K$ and bandwidth $a_N$. Consequently, the Laurent series converges unconditionally and geometrically with base ratio $(a_N/T)^2 < 1$. Away from non-generic algebraic cancellations among the endpoint derivatives, successive non-vanishing terms scale with this characteristic geometric step ratio.*

*Proof.*
1. **Divided-Difference Integral Identity and Constant Normalization:**
   In the Connes–van Suijlekom Galerkin discretization [2, 3], the Archimedean coordinate function $\psi_{\mathrm{arch}}^{(T)}(x)$ is defined by integrating the smooth test-function transform $S_x(r) = \int_{-L/2}^{L/2} e^{i r y} \frac{\cos(2\pi x y / L)}{L/2} dy$ over $[-T, T]$ against the Archimedean weight $h_+(r)$:
   $$\psi_{\mathrm{arch}}^{(T)}(x) = \frac{1}{2\pi^2} \int_{-T}^T h_+(r) \operatorname{Re} S_x(r) \, dr.$$
   The matrix entries $(Q_{\mathrm{arch}}^{(T)})_{mn}$ are defined by divided differences for $m \ne n$ and derivatives for $m = n$. Contracting with the full symmetric Fourier coefficient vector $u \in \mathbb{R}^{2N+1}$ ($u_{-m} = u_m$) and interchanging the finite sum with the compact Riemann integral on $[-T, T]$:
   $$u^T Q_{\mathrm{arch}}^{(T)} u = \sum_{m,n=-N}^N u_m u_n (Q_{\mathrm{arch}}^{(T)})_{mn} = \frac{1}{2\pi^2} \int_{-T}^T h_+(r) \left[ \sum_{m,n=-N}^N u_m u_n \frac{\operatorname{Re} S_m(r) - \operatorname{Re} S_n(r)}{m - n} \right] dr.$$
   For any symmetric vector $u$, the discrete divided-difference bilinear form inside the brackets evaluates identically to $\pi K_{\mathrm{Fourier}}(v, r, L)$ (Theorem 4.1 & Theorem 4.3). Because both $h_+(r)$ and $K_{\mathrm{Fourier}}(v, r, L)$ are even functions of $r$, the integral over $[-T, T]$ is exactly $2 \int_0^T$:
   $$\frac{1}{2\pi^2} \int_{-T}^T h_+(r) \Big[ \pi K_{\mathrm{Fourier}}(v, r, L) \Big] dr = \frac{1}{2\pi^2} \cdot 2 \pi \int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr = \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr.$$
   The normalization constants cancel identically: $\frac{1}{2\pi^2} \times \pi \times 2 = \frac{1}{\pi}$, establishing identity (1) unconditionally for every $T > 0$.

2. **Exact Archimedean Cutoff Defect:**
   The full discrete Galerkin matrix is $Q^{(T)} = Q_{\mathrm{arch}}^{(T)} + Q_{\mathrm{pole}} + Q_{\mathrm{prime}}$. Since the pole and prime matrices have no frequency cutoff and match their continuous functional counterparts identically ($u_N^T Q_{\mathrm{pole}} u_N = \mathcal{Q}_{\mathrm{pole}}(v_N)$ and $u_N^T Q_{\mathrm{prime}} u_N = \mathcal{Q}_{\mathrm{prime}}(v_N)$), subtracting the exact continuous functional $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{(\infty)}(v_N)$ gives:
   $$\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) = u_N^T Q_{\mathrm{arch}}^{(T)} u_N - \mathcal{Q}_{\mathrm{arch}}^{(\infty)}(v_N) = \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr - \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr.$$
   Combining the integrals yields $-\frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr = -\delta_T^{\mathrm{tail}}(v_N)$.

3. **Taylor Endpoint-Jet Laurent Decomposition:**
   For $r \ge T > a_N$, we have $a_m / r \le a_N / T < 1$ for all $m \in \{1, \dots, N\}$. Expanding each term in the rational resolvent (Theorem 3.1):
   $$\frac{r v_m}{r^2 - a_m^2} = \frac{v_m}{r} \frac{1}{1 - a_m^2 / r^2} = \sum_{j=0}^\infty \frac{a_m^{2j} v_m}{r^{2j+1}}.$$
   Summing over $m$ and adding $v_0/r$:
   $$\frac{v_0}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} = \frac{1}{r} \left[ v_0 + \sqrt{2}\sum_{m=1}^N v_m + \sum_{j=1}^\infty \frac{\sqrt{2}\sum_{m=1}^N a_m^{2j} v_m}{r^{2j}} \right] = \sum_{j=0}^\infty \frac{(-1)^j D_j}{r^{2j+1}},$$
   where the endpoint derivatives (Theorem 3.2) are explicitly given by:
   $$D_0 \equiv T_v(0) = v_0 + \sqrt{2}\sum_{m=1}^N v_m, \qquad D_j \equiv T_v^{(2j)}(0) = (-1)^j \sqrt{2}\sum_{m=1}^N a_m^{2j} v_m \quad (j \ge 1).$$
   (Note that the constant mode $v_0$ has vanishing derivatives of positive order, so $v_0$ contributes strictly to $D_0$ and does not appear in $D_j$ for $j \ge 1$.)
   Squaring this series and taking the Cauchy product yields:
   $$R_v(r) = \frac{2}{L} \left[ \sum_{j=0}^\infty \frac{(-1)^j D_j}{r^{2j+1}} \right]^2 = \sum_{k=0}^\infty \frac{A_k(N)}{r^{2k+2}}, \qquad A_k(N) = \frac{2}{L} (-1)^k \sum_{j=0}^k D_j D_{k-j}.$$
   Since $K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r)$ (Theorem 4.1), substituting into the tail integral gives $\delta_T^{\mathrm{tail}}(v_N) = \sum_{k=0}^\infty A_k(N) \mathcal{J}_k(T, L)$. Term-by-term integration is justified by the dominated convergence theorem, as the series converges uniformly on $[T, \infty)$ and $h_+(r)(1 - \cos(rL)) = \mathcal{O}(\log r)$.

4. **Explicit Geometric Remainder Bound:**
   To establish the explicit remainder bound with full parameter accounting, let $r \ge T > a_N$. Since $a_m \le a_N$ for all $m \in \{1, \dots, N\}$, the endpoint derivatives satisfy the uniform upper bounds:
   $$|D_0| = \left| v_0 + \sqrt{2}\sum_{m=1}^N v_m \right| \le |v_0| + \sqrt{2}\sum_{m=1}^N |v_m| \le \sqrt{2} \sum_{m=0}^N |v_m| = \sqrt{2} \|v\|_1,$$
   and for each $j \ge 1$:
   $$|D_j| = \sqrt{2} \left| \sum_{m=1}^N a_m^{2j} v_m \right| \le \sqrt{2} a_N^{2j} \sum_{m=1}^N |v_m| \le \sqrt{2} \|v\|_1 a_N^{2j}.$$
   Thus the bound $|D_j| \le \sqrt{2} \|v\|_1 a_N^{2j}$ holds uniformly for all $j \ge 0$ (with $a_N^0 \equiv 1$).
   The Laurent coefficients $A_k(N) = \frac{2}{L}(-1)^k \sum_{j=0}^k D_j D_{k-j}$ are the Cauchy convolution of the endpoint jets. By the triangle inequality:
   $$|A_k(N)| \le \frac{2}{L} \sum_{j=0}^k |D_j D_{k-j}| \le \frac{2}{L} \sum_{j=0}^k \Big( \sqrt{2} \|v\|_1 a_N^{2j} \Big) \Big( \sqrt{2} \|v\|_1 a_N^{2(k-j)} \Big) = \frac{4 (k + 1)}{L} \|v\|_1^2 a_N^{2k}.$$
   Subtracting the partial sum of order $K$ from the rational resolvent $R_v(r) = \sum_{k=0}^\infty \frac{A_k(N)}{r^{2k+2}}$:
   $$\left| R_v(r) - \sum_{k=0}^K \frac{A_k(N)}{r^{2k+2}} \right| \le \sum_{k=K+1}^\infty \frac{|A_k(N)|}{r^{2k+2}} \le \frac{4 \|v\|_1^2}{L r^2} \sum_{k=K+1}^\infty (k+1) \left( \frac{a_N}{r} \right)^{2k}.$$
   Summing the series for $x = (a_N/r)^2 \le (a_N/T)^2 < 1$:
   $$\sum_{k=K+1}^\infty (k+1) x^k = x^{K+1} \sum_{l=0}^\infty (K + 2 + l) x^l = x^{K+1} \frac{K + 2 - (K+1)x}{(1 - x)^2} \le \frac{(K+2) x^{K+1}}{1 - x}.$$
   Therefore:
   $$\left| R_v(r) - \sum_{k=0}^K \frac{A_k(N)}{r^{2k+2}} \right| \le \frac{4 (K+2) \|v\|_1^2}{L} \frac{a_N^{2K+2}}{r^{2K+4}} \frac{1}{1 - (a_N/T)^2}.$$
   Substituting into the tail integral $\delta_T^{\mathrm{tail}}(v_N) = \frac{1}{\pi}\int_T^\infty h_+(r)(1 - \cos(rL)) R_v(r) dr$, and using $0 \le 1 - \cos(rL) \le 2$ and $h_+(r) \le \log(r/2) + \frac{1}{24 r^2} \le \log r$ for $r \ge T \ge 2\pi$:
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{8 (K+2) \|v\|_1^2}{\pi L} \frac{a_N^{2K+2}}{1 - (a_N/T)^2} \int_T^\infty \frac{\log r}{r^{2K+4}} \, dr.$$
   Integrating by parts:
   $$\int_T^\infty \frac{\log r}{r^{2K+4}} \, dr = \frac{\log T}{(2K+3) T^{2K+3}} + \frac{1}{(2K+3)^2 T^{2K+3}} \le \frac{\log T + 1}{(2K+3) T^{2K+3}}.$$
   Since $\frac{K+2}{2K+3} \le \frac{2}{3} < 1$, combining factors yields:
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{C_K(v, N, L, T)}{T^{2K+3}} \frac{1}{1 - (a_N/T)^2} = \frac{C_0(v, L, T)}{T} \left( \frac{a_N}{T} \right)^{2K+2} \frac{1}{1 - (a_N/T)^2},$$
   where $C_K(v, N, L, T) \equiv \frac{8(\log T + 1)}{\pi L} \|v\|_1^2 a_N^{2K+2}$ and $C_0(v, L, T) \equiv \frac{8(\log T + 1)}{\pi L} \|v\|_1^2$. This completes the rigorous proof. $\blacksquare$

### Remark 5.6 (Leading Non-Oscillatory Asymptotics and Alternating Jet Compensation)
At leading order $k = 0$, evaluating the non-oscillatory part of $\mathcal{J}_0(T, L) = \frac{1}{\pi} \int_T^\infty \frac{h_+(r)}{r^2} dr + \text{oscillatory}$ using $h_+(r) = \log(r/(2\pi)) + \mathcal{O}(r^{-2})$ yields the leading asymptotic formula:
$$\mathcal{J}_0(T, L) \sim \frac{1}{\pi T} \left( \log \frac{T}{2\pi} + 1 \right) \implies A_0 \mathcal{J}_0(T, L) \sim \frac{2 D_0^2}{\pi L T} \left( \log \frac{T}{2\pi} + 1 \right).$$
For $N = 24$, $c = 13$, $L = \log 13 \approx 2.565$, and $T = 400$, this leading term gives $A_0 \mathcal{J}_0 \approx 4.08 \times 10^{-43}$, which exceeds the observed defect $\delta_T^{\mathrm{tail}} \approx 1.67 \times 10^{-43}$ by a factor of $\approx 2.45$.

This difference is resolved by the sub-leading jets: while the rigorous remainder bound (4) guarantees overall geometric convergence governed by the spectral radius $a_N$, individual coefficients $A_k = \frac{2}{L}(-1)^k \sum_{j=0}^k D_j D_{k-j}$ can experience internal cancellations between the bulk and boundary modes. Indeed, because the first-jet ratio $D_1/D_0 \approx 5.20 \times 10^5$ is large (Theorem 7.1), the first sub-leading Laurent coefficient $A_1 = -\frac{4}{L} D_0 D_1 \approx -1.05 \times 10^{-34}$ is negative and substantially larger in absolute magnitude, contributing a substantial negative correction $A_1 \mathcal{J}_1 \approx -7.78 \times 10^{-43}$. The second jet $A_2$ contributes $+1.09 \times 10^{-42}$, and successive alternating terms converge geometrically with the characteristic step ratio $(a_{24}/T)^2 = (58.788 / 400)^2 \approx 0.0216 \approx 1/46.3$. At truncation order $K = 8$, the jet sum matches the target defect to $2.23 \times 10^{-45}$, confirming that the tiny residual is governed by an exquisitely structured cancellation among the endpoint jets.

---

## 6. Exact Commutator Algebra, Parity Reduction, and Odd-Sector Resolvent Identity

We now analyze the algebraic structure of the Connes–van Suijlekom Galerkin matrix $Q$ directly via its commutators with the coordinate operator $M = \operatorname{diag}(-N, \dots, N)$.

### Theorem 6.1 (General Commutator Rank and Parity Factorization)
*Let $M = \operatorname{diag}(-N, \dots, N)$ on $\mathbb{C}^{2N+1}$, and let $Q$ be the Connes–van Suijlekom Galerkin matrix in full Fourier coordinates:*

$$Q_{mn} = \begin{cases} \dfrac{\psi(m) - \psi(n)}{m - n}, & m \ne n, \\[8pt] \psi'(n), & m = n, \end{cases}$$

*where $\psi(-n) = -\psi(n)$ and $\psi'(-n) = \psi'(n)$. For every integer $k \ge 1$, the commutator $[M^k, Q]$ has rank at most $2k$ and admits the exact algebraic representation:*

$$[M^k, Q] = \sum_{j=0}^{k-1} \Big( (M^j p)(M^{k-1-j} e)^T - (M^j e)(M^{k-1-j} p)^T \Big),$$

*where $e = (1, \dots, 1)^T$ and $p = (\psi(-N), \dots, \psi(N))^T$.*

*Proof.* For $m \ne n$, divided differences give $\frac{m^k - n^k}{m - n} = \sum_{j=0}^{k-1} m^j n^{k-1-j}$. Multiplying by $\psi(m) - \psi(n)$:
$$[M^k, Q]_{mn} = (m^k - n^k) \frac{\psi(m) - \psi(n)}{m - n} = \sum_{j=0}^{k-1} \left( m^j \psi(m) \cdot n^{k-1-j} - m^j \cdot n^{k-1-j} \psi(n) \right).$$
On the diagonal $m = n$, both sides vanish identically, proving the outer-product representation. $\blacksquare$

### Corollary 6.1.1 (Strict Parity Decoupling)
*Under the parity reflection operator $\mathcal{P} x_n = x_{-n}$, $e$ is even ($\mathcal{P} e = e$), $p$ is odd ($\mathcal{P} p = -p$), and $M$ is odd ($\mathcal{P} M \mathcal{P} = -M$). Consequently, for any even vector $u$ ($\mathcal{P} u = u$):*

$$(M^r e)^T u = \begin{cases} 0, & r \text{ odd}, \\[6pt] (-1)^{r/2} \dfrac{D_{r/2}}{\kappa^r}, & r \text{ even}, \end{cases} \qquad (M^r p)^T u = \begin{cases} B_r \equiv \displaystyle\sum_{n=-N}^N n^r \psi(n) u_n, & r \text{ odd}, \\[8pt] 0, & r \text{ even}. \end{cases}$$

*Every inner product appearing in $[M^k, Q] u$ reduces strictly to either an endpoint jet $D_m$ or an arithmetic moment $B_r$.*

### Theorem 6.2 (Exact Odd-Sector Resolvent Identity for $B_1$)
*Let $u$ be the even ground-state eigenvector satisfying $Q u = \lambda u$ with $\lambda = \lambda_{\min}(N) > 0$. Assume the odd-sector positive-gap condition:*

$$\lambda < \min \operatorname{spec}(Q_{\mathrm{odd}}),$$

*so that $Q_{\mathrm{odd}} - \lambda I$ is strictly positive definite on $\mathcal{H}_{\mathrm{odd}}$ and boundedly invertible.*

1. **Odd-Sector Resolvent Identity:**
   *The $k = 1$ commutator applied to $u$ yields:*
   $$M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi.$$
2. **Arithmetic Energy Proportionality:**
   *Consequently, the arithmetic moment scalar $B_1$ is strictly proportional to $D_0$:*
   $$B_1 = -D_0 \mathcal{E}_{\mathrm{arith}}, \qquad \mathcal{E}_{\mathrm{arith}} \equiv \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle > 0.$$

*Proof.* Setting $k = 1$ in Theorem 6.1 gives $[M, Q] = p e^T - e p^T$. Applying to $u$ and using $p^T u = 0$ yields $[M, Q] u = D_0 \psi$. Since $[M, Q] u = \lambda M u - Q M u$, we obtain $(Q - \lambda I)(M u) = -D_0 \psi$. Because $M u$ and $\psi$ are purely odd, and because $\lambda < \min \operatorname{spec}(Q_{\mathrm{odd}})$ ensures that $(Q_{\mathrm{odd}} - \lambda I)$ is strictly positive definite and boundedly invertible on $\mathcal{H}_{\mathrm{odd}}$, inverting on $\mathcal{H}_{\mathrm{odd}}$ establishes $M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi$. Taking the inner product with $\psi$ yields $B_1 = \psi^T M u = -D_0 \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle = -D_0 \mathcal{E}_{\mathrm{arith}}$. Since $Q_{\mathrm{odd}} - \lambda I > 0$ and $\psi \ne 0$, $\mathcal{E}_{\mathrm{arith}} > 0$ is strictly positive. $\blacksquare$

### Proposition 6.3 (Exact $K$-Commutator and First-Jet Resolvent Norm Identity)
*Let $K \equiv M = \operatorname{diag}(-N, \dots, N)$ denote the Fourier-index operator, $d \equiv e = (1, \dots, 1)^T$, $\boldsymbol\psi \equiv p = (\psi(-N), \dots, \psi(N))^T$, and let $c \equiv u$ denote the normalized even ground state ($Q c = \lambda c$).*

1. **Exact Rank-2 Commutator:**
   *The commutator $[Q, K]$ satisfies the exact rank-2 identity on $\mathbb{C}^{2N+1}$:*
   $$[Q, K] = -\boldsymbol\psi d^T + d \boldsymbol\psi^T.$$
2. **Ground-State Action:**
   *Acting on the even ground state $c$, because $\boldsymbol\psi^T c = 0$ and $d^T c = D_0$:*
   $$(Q - \lambda I) K c = -D_0 \boldsymbol\psi.$$
3. **Odd-Sector Resolvent Expansion & Norm Identity:**
   *Expanding along the orthonormal odd-sector eigenbasis $\{e_j\}_{j \ge 1}$ of $Q_{\mathrm{odd}}$ ($Q e_j = \lambda_j e_j$):*
   $$K c = -D_0 \sum_{j \ge 1} \frac{\langle \boldsymbol\psi, e_j \rangle}{\lambda_j - \lambda} e_j.$$
   *Consequently, taking Euclidean norms yields the exact closed-form identity:*
   $$\|K c\|^2 = D_0^2 \sum_{j \ge 1} \frac{|\langle \boldsymbol\psi, e_j \rangle|^2}{(\lambda_j - \lambda)^2} = D_0^2 \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-2} \boldsymbol\psi \rangle \equiv D_0^2 M_2,$$
   *where $M_2 \equiv \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-2} \boldsymbol\psi \rangle$.*

*Proof.* For $m \ne n$, $[Q, K]_{mn} = Q_{mn}(n - m) = \frac{\psi(m) - \psi(n)}{m - n}(n - m) = -(\psi(m) - \psi(n))$. The rank-2 outer product $(-\boldsymbol\psi d^T + d \boldsymbol\psi^T)_{mn} = -\psi(m) \cdot 1 + 1 \cdot \psi(n) = -(\psi(m) - \psi(n))$. On the diagonal $m = n$, both sides vanish identically. Applying to $c$ and noting $\boldsymbol\psi^T c = 0$ gives $(Q - \lambda I) K c = [Q, K] c = -\boldsymbol\psi (d^T c) = -D_0 \boldsymbol\psi$. Since $K c$ and $\boldsymbol\psi$ are purely odd, inverting $Q - \lambda I$ on $\mathcal{H}_{\mathrm{odd}}$ and taking the squared norm yields $\|K c\|^2 = D_0^2 \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-2} \boldsymbol\psi \rangle$. $\blacksquare$

### Corollary 6.3.1 (Exact Mode-by-Mode Commutator Projections and Parseval Decomposition)
*For every orthonormal odd eigenmode $e_j \in \mathcal{H}_{\mathrm{odd}}$ with eigenvalue $\mu_j$ and odd spectral gap $\Delta_j \equiv \mu_j - \lambda$, the projection of the coordinate action $Kc$ onto $e_j$ satisfies the exact identity:*

$$\langle e_j, K c \rangle \equiv - D_0 \frac{\langle \boldsymbol\psi, e_j \rangle}{\mu_j - \lambda} \equiv - D_0 \frac{a_j}{\Delta_j}, \qquad a_j \equiv \langle \boldsymbol\psi, e_j \rangle.$$

*Consequently, the Euclidean norm of $Kc$ satisfies the exact Parseval decomposition:*

$$\|K c\|^2 = \sum_{j=1}^N |\langle e_j, K c \rangle|^2 = D_0^2 \sum_{j=1}^N \frac{a_j^2}{\Delta_j^2} \equiv D_0^2 M_2.$$

*In particular, for the lowest odd mode $e_1$, the coordinate projection factors into the product of the two dimensionless tunneling invariants:*

$$|\langle e_1, K c \rangle| = \left(\frac{D_0}{\sqrt{\Delta_1}}\right) \left(\frac{|a_1|}{\sqrt{\Delta_1}}\right) \equiv R_D \cdot C_N,$$

*where $R_D \equiv D_0 / \sqrt{\Delta_1}$ and $C_N \equiv |a_1| / \sqrt{\Delta_1}$.*

*Proof.* Taking the inner product of $(Q - \lambda I) K c = -D_0 \boldsymbol\psi$ (Proposition 6.3) with $e_j$ yields:
$$\langle e_j, (Q - \lambda I) K c \rangle = \langle (Q - \lambda I) e_j, K c \rangle = (\mu_j - \lambda) \langle e_j, K c \rangle = -D_0 \langle e_j, \boldsymbol\psi \rangle = -D_0 a_j.$$
Dividing by $\Delta_j = \mu_j - \lambda \ne 0$ establishes the projection identity. Squaring and summing over the complete orthonormal basis $\{e_j\}_{j=1}^N$ of $\mathcal{H}_{\mathrm{odd}}$ yields the Parseval formula. Factoring $\Delta_1 = \sqrt{\Delta_1}\sqrt{\Delta_1}$ in the $j=1$ term gives $|\langle e_1, K c \rangle| = R_D \cdot C_N$. $\blacksquare$

*Remark 6.3.1 (Finite-Rank Commutator Hierarchy).*
The appearance of low-rank commutators is an inherent consequence of the divided-difference structure: for any power $K^p$, $\frac{m^p - n^p}{m - n}$ is a polynomial of degree $p-1$ in $m$ and $n$, ensuring that $[Q, K^p]$ has rank at most $2p$.

---

## 7. The First-Jet Ratio $D_1/D_0$ and Exact Resolvent Small-Denominator Cancellation

We now establish the exact algebraic structure of the first-jet cancellation ratio $D_1 / D_0$ and the non-singular cancellation of small denominators in the resolvent spectral expansion.

### Theorem 7.1 (Exact Archimedean Resolvent Jet Representation)
*The first-jet cancellation ratio $D_1 / D_0$ is identically equal to the relative first correction of the large-$r$ Archimedean resolvent:*

$$\frac{D_1}{D_0} \equiv -\frac{1}{2} \frac{A_1}{A_0} \equiv \frac{D'(0)}{D(0)} \equiv -\kappa^2 \frac{F'(0)}{F(0)} \equiv -\frac{\sqrt{2}\displaystyle\sum_{m=1}^N a_m^2 v_m}{v_0 + \sqrt{2}\displaystyle\sum_{m=1}^N v_m},$$

*where $D(w) \equiv \big[(I + w\mathcal{L})^{-1} T_v\big](0) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 w}$ is the operator resolvent generating function, and $F(z) = e^T (I - z M^2)^{-1} v = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - z m^2}$ is the discrete mode generating function. The two generating functions are related identically by:*

$$D(w) \equiv F(-\kappa^2 w), \qquad \kappa = \frac{2\pi}{L}.$$

*Proof.* From Theorem 3.1, the large-$r$ expansion of the reduced Fourier kernel is:

$$R_v(r) = \frac{A_0}{r^2} + \frac{A_1}{r^4} + \mathcal{O}(r^{-6}),$$

where $A_0 = \frac{2}{L} D_0^2$ and $A_1 = -\frac{4}{L} D_0 D_1$. Taking the ratio gives $-A_1 / (2 A_0) = -(-\frac{4}{L} D_0 D_1) / (\frac{4}{L} D_0^2) = D_1 / D_0$. 

In terms of the operator resolvent generating function $D(w) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 w}$, Taylor expansion around $w = 0$ gives $D(w) = D_0 + D_1 w + \mathcal{O}(w^2)$, whence $D(0) = D_0 = v_0 + \sqrt{2}\sum_{m=1}^N v_m$ and $D'(0) = D_1 = -\sqrt{2}\sum_{m=1}^N a_m^2 v_m$. Thus:

$$\frac{D'(0)}{D(0)} = \frac{D_1}{D_0} = -\frac{\sqrt{2}\sum_{m=1}^N a_m^2 v_m}{v_0 + \sqrt{2}\sum_{m=1}^N v_m}.$$

In terms of the unscaled coordinate generating function $F(z) = e^T (I - z M^2)^{-1} v = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - z m^2}$, expanding around $z = 0$ gives $F(z) = D_0 + (\sqrt{2}\sum_{m=1}^N m^2 v_m) z + \mathcal{O}(z^2) = D_0 - \frac{D_1}{\kappa^2} z + \mathcal{O}(z^2)$. Differentiating at $z = 0$ yields $F'(0) = -D_1 / \kappa^2$, whence $F'(0) / F(0) = -D_1 / (\kappa^2 D_0)$, which gives $\frac{D_1}{D_0} = -\kappa^2 \frac{F'(0)}{F(0)}$. 

Equivalently, substituting $a_m = \kappa m$ into $D(w)$ gives $D(w) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - (-\kappa^2 w) m^2} \equiv F(-\kappa^2 w)$. Differentiating via the chain rule gives $D'(0) = -\kappa^2 F'(0)$, immediately establishing $\frac{D'(0)}{D(0)} = -\kappa^2 \frac{F'(0)}{F(0)} = \frac{D_1}{D_0}$. $\blacksquare$

*Remark 7.1.1 (Algebraic Ratio vs. Measure Expectation).*
In terms of the signed atomic measure $\mu_N = v_0 \delta_0 + \sqrt{2}\sum_{m=1}^N v_m \delta_{a_m^2}$ supported on the discrete Neumann Laplacian spectrum $\{0, a_1^2, \dots, a_N^2\}$, this ratio evaluates formally as $-\frac{\int_0^\infty x \, d\mu_N(x)}{\int_0^\infty d\mu_N(x)}$. Because the Fourier coefficients $v_m$ alternate and take negative values, $\mu_N$ is a signed measure rather than a positive probability distribution. The identity is a strictly algebraic relation between finite sums, not a probabilistic expectation.

---

### Theorem 7.2 (Exact Denominator Cancellation in the Odd-Even Spectral Coupling Away from Resonances)
*Let $Q_{\mathrm{even}}$ and $Q_{\mathrm{odd}}$ denote the restrictions of the Connes–van Suijlekom Galerkin matrix $Q$ to the even and odd subspaces of $\mathbb{C}^{2N+1}$. For any even eigenvector $u^{(k)}$ of $Q_{\mathrm{even}}$ with eigenvalue $E_k$ satisfying the non-resonance condition $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$ and $\lambda \notin \operatorname{spec}(Q_{\mathrm{odd}})$, let $D_0^{(k)} = e^T u^{(k)}$ and define the arithmetic energy vector $\mathbf{w}(\lambda) = b + \mathcal{E}_{\mathrm{arith}}(\lambda) e$, where $b_n = n \psi(n)$ and $\mathcal{E}_{\mathrm{arith}}(\lambda) = \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1}\psi \rangle$.*

*Then in the spectral expansion of $(Q_{\mathrm{even}} - \lambda I)^\dagger \mathbf{w}(\lambda)$, the numerator for mode $k$ factors identically as:*

$$u^{(k)T} \mathbf{w}(\lambda) = -(E_k - \lambda) D_0^{(k)} \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$

*Consequently, for every $E_k \ne \lambda$ with $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$, the factor $(E_k - \lambda)$ cancels identically, yielding the exact algebraic summand with the $(E_k - \lambda)$ denominator eliminated:*

$$\frac{(e^T u^{(k)}) (u^{(k)T} \mathbf{w}(\lambda))}{E_k - \lambda} \equiv - [D_0^{(k)}]^2 \cdot \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$

*Proof.* Expanding $(Q_{\mathrm{even}} - \lambda I)^\dagger$ along the orthonormal eigenbasis $\{u^{(k)}\}_{k \ge 1, \text{ even}}$ gives numerators $(e^T u^{(k)}) (u^{(k)T} \mathbf{w}(\lambda)) = D_0^{(k)} [B_1^{(k)} + \mathcal{E}_{\mathrm{arith}}(\lambda) D_0^{(k)}]$, where $B_1^{(k)} = \psi^T M u^{(k)}$. Since $[M, Q] u^{(k)} = D_0^{(k)} \psi$ by Theorem 6.1 and $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$, the operator $(Q_{\mathrm{odd}} - E_k I)$ is boundedly invertible on $\mathcal{H}_{\mathrm{odd}}$, yielding $M u^{(k)} = -D_0^{(k)} (Q_{\mathrm{odd}} - E_k I)^{-1} \psi$, whence $B_1^{(k)} = \psi^T M u^{(k)} = -D_0^{(k)} \mathcal{E}_{\mathrm{arith}}(E_k)$. The bracket is therefore:

$$B_1^{(k)} + \mathcal{E}_{\mathrm{arith}}(\lambda) D_0^{(k)} = -D_0^{(k)} \big[ \mathcal{E}_{\mathrm{arith}}(E_k) - \mathcal{E}_{\mathrm{arith}}(\lambda) \big].$$

By the first resolvent identity for the odd operator $Q_{\mathrm{odd}}$ (valid since $E_k, \lambda \notin \operatorname{spec}(Q_{\mathrm{odd}})$):

$$(Q_{\mathrm{odd}} - E_k I)^{-1} - (Q_{\mathrm{odd}} - \lambda I)^{-1} = (E_k - \lambda) (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1}.$$

Taking the inner product with $\psi$:

$$\mathcal{E}_{\mathrm{arith}}(E_k) - \mathcal{E}_{\mathrm{arith}}(\lambda) = (E_k - \lambda) \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$

Substituting this into the numerator yields $u^{(k)T} \mathbf{w}(\lambda) = -(E_k - \lambda) D_0^{(k)} \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle$. Dividing by $(E_k - \lambda)$ completes the algebraic cancellation. $\blacksquare$

*Remark 7.2.1 (Non-Resonance and Spectral Disjointness).*
The hypothesis $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$ ensures that the odd resolvent $(Q_{\mathrm{odd}} - E_k I)^{-1}$ is well-defined. Theorem 7.2 eliminates the explicit $(E_k - \lambda)$ denominator algebraically, but the odd-sector resolvent $(Q_{\mathrm{odd}} - E_k I)^{-1}$ possesses genuine poles if $E_k \in \operatorname{spec}(Q_{\mathrm{odd}})$. In all numerical Galerkin calculations across all tested dimensions and cutoffs ($N \in \{1, \dots, 24\}$, $c \in \{5, \dots, 17\}$), the even and odd spectra are strictly disjoint ($\operatorname{spec}(Q_{\mathrm{even}}) \cap \operatorname{spec}(Q_{\mathrm{odd}}) = \emptyset$), so no spectral degeneracies or resonances occur between the sectors.

---

### Theorem 7.3 (Exact $K^2$-Commutator Resolvent Representation of $D_1/D_0$)
*Let $K = \operatorname{diag}(-N, \dots, N)$, $d = (1, \dots, 1)^T$, $k = K d$, and let $c$ be the normalized even ground state satisfying $Q c = \lambda c$. Assume that the ground-state eigenvalue $\lambda = \lambda_{\min}(N)$ is simple in the even sector, and let $R_{\mathrm{even}} \equiv (Q_{\mathrm{even}} - \lambda I)_{c^\perp}^{-1}$ denote the inverse of $Q_{\mathrm{even}} - \lambda I$ restricted to the orthogonal complement $c^\perp$.*

1. **Rank-4 Commutator Identity:**
   *The commutator $[Q, K^2]$ admits the exact rank-4 representation:*
   $$[Q, K^2] = -(K\boldsymbol\psi) d^T - \boldsymbol\psi k^T + k \boldsymbol\psi^T + d (K\boldsymbol\psi)^T.$$
2. **Ground-State Action and Automatic Orthogonality:**
   *Applying $[Q, K^2]$ to $c$ yields:*
   $$(Q - \lambda I) K^2 c = -D_0 (K\boldsymbol\psi + M_1 d) \equiv -D_0 s_2,$$
   *where $M_1 \equiv \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi \rangle$. The source vector $s_2$ is identically orthogonal to the ground state:*
   $$\langle c, s_2 \rangle = \langle c, K\boldsymbol\psi \rangle + M_1 \langle c, d \rangle = -D_0 M_1 + M_1 D_0 \equiv 0.$$
3. **Exact Resolvent Representation of $D_1/D_0$:**
   *Consequently, $K^2 c = D_0^2 M_2 c - D_0 R_{\mathrm{even}} s_2$, and contracting with $d^T$ yields the exact closed form:*
   $$\boxed{\frac{D_1}{D_0} = \kappa^2 \left[ \langle d, R_{\mathrm{even}}(K\boldsymbol\psi + M_1 d) \rangle - D_0^2 M_2 \right],} \qquad \kappa = \frac{2\pi}{L}.$$

*Proof.* For $m \ne n$, $[Q, K^2]_{mn} = Q_{mn}(n^2 - m^2) = -(\psi(m) - \psi(n))(m + n)$. The $(m, n)$ entry of the right-hand side is $-(m \psi(m)) \cdot 1 - \psi(m) \cdot n + m \cdot \psi(n) + 1 \cdot (n \psi(n)) = -(m + n)\psi(m) + (m + n)\psi(n) = -(\psi(m) - \psi(n))(m + n)$, and on the diagonal both sides vanish, proving statement 1.

Applying to $c$, parity gives $k^T c = 0$ and $\boldsymbol\psi^T c = 0$. Meanwhile $d^T c = D_0$, and by Proposition 6.3, $(K\boldsymbol\psi)^T c = \boldsymbol\psi^T K c = -D_0 \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \boldsymbol\psi \rangle = -D_0 M_1$. Thus $[Q, K^2] c = -D_0 K\boldsymbol\psi - D_0 M_1 d = -D_0 s_2$, with $\langle c, s_2 \rangle = -D_0 M_1 + M_1 D_0 = 0$, proving statement 2.

Because $s_2 \in c^\perp$ and $\lambda$ is simple in the even sector, $(Q_{\mathrm{even}} - \lambda I)$ is invertible on $c^\perp$ with resolvent $R_{\mathrm{even}}$. Writing $K^2 c = \alpha c - D_0 R_{\mathrm{even}} s_2$, where $\alpha = c^T K^2 c = \|K c\|^2 = D_0^2 M_2$ by Proposition 6.3, taking the inner product with $d$ yields $d^T K^2 c = \alpha (d^T c) - D_0 \langle d, R_{\mathrm{even}} s_2 \rangle = D_0^3 M_2 - D_0 \langle d, R_{\mathrm{even}} s_2 \rangle$. Recalling from Theorem 7.1 that $d^T K^2 c = \sum m^2 c_m = \sqrt{2}\sum_{m=1}^N m^2 v_m = -D_1 / \kappa^2$, we obtain $-D_1 / \kappa^2 = D_0^3 M_2 - D_0 \langle d, R_{\mathrm{even}} s_2 \rangle$. Multiplying by $-\kappa^2$ and dividing by $D_0$ proves statement 3. $\blacksquare$

### Corollary 7.3.1 (Exact Even-Sector Mode Projections and Small-Denominator Cancellation)
*For every excited even eigenmode $u^{(k)} \in \mathcal{H}_{\mathrm{even}} \cap c^\perp$ ($k \ge 1$) with eigenvalue $E_k$ and gap $\Delta_{\mathrm{even}, k} \equiv E_k - \lambda$, the projection of $K^2 c$ satisfies the exact identity:*

$$\langle u^{(k)}, K^2 c \rangle \equiv - D_0 \frac{\langle u^{(k)}, s_2 \rangle}{E_k - \lambda} \equiv - D_0 \frac{b_k}{\Delta_{\mathrm{even}, k}}, \qquad b_k \equiv \langle u^{(k)}, s_2 \rangle.$$

*Consequently, the excited-state contribution to the resolvent inner product in Theorem 7.3 admits the exact representation:*

$$\tau_k \equiv \frac{d_k b_k}{\Delta_{\mathrm{even}, k}} \equiv - \frac{d_k \langle u^{(k)}, K^2 c \rangle}{D_0}, \qquad d_k \equiv \langle u^{(k)}, d \rangle.$$

*The small denominator $\Delta_{\mathrm{even}, k}$ is algebraically extinguished by $b_k \propto \Delta_{\mathrm{even}, k}$, leaving each spectral summand $\tau_k$ finite and directly determined by the projection $\langle u^{(k)}, K^2 c \rangle$.*

*Proof.* Taking the inner product of $(Q - \lambda I) K^2 c = -D_0 s_2$ from Theorem 7.3 with $u^{(k)}$ yields $(E_k - \lambda) \langle u^{(k)}, K^2 c \rangle = -D_0 \langle u^{(k)}, s_2 \rangle = -D_0 b_k$. Dividing by $\Delta_{\mathrm{even}, k}$ gives the projection identity. Substituting $b_k / \Delta_{\mathrm{even}, k} = -\langle u^{(k)}, K^2 c \rangle / D_0$ into $\tau_k = d_k b_k / \Delta_{\mathrm{even}, k}$ immediately yields the formula. $\blacksquare$

*Remark 7.3.1 (Elimination of the Exponential Ground-State Singularity).*
In Theorem 7.3, the exponentially tiny tunneling factor $D_0$ factors out algebraically. The catastrophic ground-state singularity $1/\lambda$ is completely eliminated from the jet hierarchy: the remaining resolvents are reduced resolvents on the odd sector and on the orthogonal complement $c^\perp$ of the even ground state. Their quantitative behaviour with $N$ is a separate spectral question investigated numerically below.

---

## 8. Conclusion and Toolkit Index

This manuscript establishes the exact, rigorous operator-theoretic foundation for the finite-rank Connes–van Suijlekom truncated Weil quadratic form. Every result in this paper is an exact algebraic identity or an unconditionally proven finite-dimensional theorem:

**Table: Index of Exact Mathematical Tools**

| Tool Name | Formal Statement | Mathematical Status | Primary Significance |
| :--- | :--- | :--- | :--- |
| **Rational Resolvent Identity** | $R_v(r) \equiv \frac{2}{L} [ \frac{v_0}{r} + \sqrt{2}\sum \frac{rv_m}{r^2 - a_m^2} ]^2$ | Theorem 3.1 | Reduces 4-term Volterra integral to single squared rational resolvent. |
| **Neumann Resolvent Identity** | $D(z) \equiv [(I + z\mathcal{L})^{-1} T_v](0)$ | Theorem 3.2 | Identifies generating function as boundary operator resolvent; heat representation. |
| **Pointwise Kernel Positivity** | $K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0$ | Theorem 4.1 | Proves Fourier Archimedean kernel is globally non-negative for all real $r$. |
| **Universal Factorization** | $\Phi_v(r) = \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D(-1/r^2)$ | Corollary 4.2 | Factors entire amplitude into universal sinc and operator resolvent. |
| **Lattice Sampling Identity** | $K_{\mathrm{Fourier}}(a_m) = \frac{L}{2} v_m^2$ | Theorem 4.3 | Samples squared Fourier coefficients orthogonally at discrete lattice frequencies. |
| **Cauchy Transform Identity** | Closed algebraic form for $J(q)$ | Theorem 5.1 | Evaluates Archimedean Cauchy transform in closed form via reflected autocorrelation and contour integration. |
| **Spatial Laplace Duality** | $J(q) = \frac{1}{L}\int_0^L K_v^{\mathrm{phys}}(L-y) e^{-qy} dy$ | Theorem 5.2 | Duality between complex frequency Cauchy transform and reflected spatial autocorrelation. |
| **Weierstrass Pole Series** | $\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum [\frac{\|v\|_2^2}{n+1} - J(q_n)]$ | Corollary 5.3 | Unconditionally convergent algebraic series with fast $\mathcal{O}(n^{-2})$ convergence. |
| **Digamma Closed Form** | $\mathcal{Q}_{\mathrm{arch}}(v) = \mathcal{Q}_{\mathrm{arch}}^{\mathrm{lattice}}(v) + \mathcal{Q}_{\mathrm{arch}}^{\mathrm{leakage}}(v)$ | Corollary 5.4 | Closed-form evaluation summing harmonic and lattice modes analytically. |
| **Lattice Lower Bound** | $\mathcal{Q}_{\mathrm{arch}}(v) \ge h_+(0) v_0^2 + \sum v_m^2 h_+(a_m)$ | Corollary 5.4.1 | Establishes discrete diagonal lattice sum as an unconditional lower bound on $\mathcal{Q}_{\mathrm{arch}}(v)$. |
| **Finite-$T$ Cutoff Defect** | $\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) \equiv -\delta_T^{\mathrm{tail}}(v_N)$ | Theorem 5.5 | Proves residual on $v_N$ is identically Archimedean cutoff tail leakage; geometric jet expansion. |
| **Rank-$2k$ Commutator** | $[M^k, Q] = \sum (M^j p)(M^{k-1-j} e)^T - (M^j e)(M^{k-1-j} p)^T$ | Theorem 6.1 | Determines exact commutator structure of Galerkin matrix. |
| **Odd-Sector Resolvent** | $Mu = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi$ | Theorem 6.2 | Solves first spectral moment and arithmetic energy uniquely in odd sector. |
| **First-Jet Norm Identity** | $\|Kc\|^2 = D_0^2 \langle \boldsymbol\psi, (Q_{\mathrm{odd}} - \lambda I)^{-2}\boldsymbol\psi \rangle$ | Proposition 6.3 | Gives the exact ratio $\|Kc\| / |D_0|$ as an odd-sector resolvent norm. |
| **Odd-Sector Mode Projections** | $\langle e_j, K c \rangle = -D_0 \frac{a_j}{\Delta_j}$, $\|Kc\|^2 = D_0^2 M_2$ | Corollary 6.3.1 | Exact coordinate projections, Parseval decomposition, and tunneling product $|\langle e_1, Kc \rangle| = R_D C_N$. |
| **Exact First-Jet Identity** | $D_1/D_0 \equiv -\frac{1}{2} A_1/A_0 \equiv D'/D \equiv -\kappa^2 F'/F$ | Theorem 7.1 | Relates endpoint derivative ratio directly to large-$r$ resolvent tail and generating functions ($D(w) \equiv F(-\kappa^2 w)$). |
| **Denominator Cancellation** | $(E_k - \lambda)$ cancels identically in resolvent coupling | Theorem 7.2 | Proves algebraic cancellation of the $(E_k - \lambda)$ denominator away from odd resonances. |
| **Exact $K^2$ Resolvent Representation** | $D_1/D_0 = \kappa^2 [\langle d, R_{\mathrm{even}} s_2 \rangle - D_0^2 M_2]$ | Theorem 7.3 | Eliminates the ground-state factor $1/\lambda$ and reduces $D_1/D_0$ to excited-sector resolvent matrix elements. |
| **Even-Sector Mode Projections** | $\langle u_k, K^2 c \rangle = -D_0 \frac{b_k}{\Delta_k}$, $\tau_k = -\frac{d_k \langle u_k, K^2 c \rangle}{D_0}$ | Corollary 7.3.1 | Exact second-jet coordinate projections and resolvent summand representation. |

---

## 9. Computational Reproducibility and Software Availability

To ensure complete computational transparency and reproducibility, the entire mathematical software pipeline and all raw high-precision calculation transcripts supporting this study are permanently archived in the public repository [10]:

> **Software Repository:** <https://github.com/akivag613/connes-cvs->  
> **Mirror Repository:** <https://github.com/nrensen/connes-cvs->

The exact algebraic identities and theorems established in this manuscript were implemented and validated using Python and the `mpmath` arbitrary-precision arithmetic library:

**Table: Mapping of Toolkit Results to Verification Scripts and Logs**

| Toolkit Theorem / Identity | Mathematical Focus | Python Script | Verification Log |
| :--- | :--- | :--- | :--- |
| Theorem 3.1 & Theorem 4.1 | 4-term Volterra reduction & squared resolvent validation | `cell32.py` | `cell32.out` |
| Theorem 3.2 | Operator resolvent & Taylor jet hierarchy $D_k = T_v^{(2k)}(0)$ | `cell38.py`, `cell43.py` | `cell38.out`, `cell43.out` |
| Theorem 4.3 | Lattice sampling identity $K_{\mathrm{Fourier}}(a_m) = \frac{L}{2} v_m^2$ | `cell32.py` | `cell32.out` |
| Theorem 5.1 & Theorem 5.2 | Exact Cauchy transform $J(q)$ & reflected Laplace duality | `cell56.py` | `cell56.out` |
| Corollary 5.3, 5.4 & 5.4.1 | Weierstrass pole series, digamma closed form & lattice lower bound | `cell56.py` | `cell56.out` |
| Theorem 5.5 | Finite-$T$ Archimedean cutoff defect & endpoint-jet reconstruction | `cell57.py` | `cell57.out` |
| Theorem 6.1 & Corollary 6.1.1 | Rank-$2k$ commutator algebra & strict parity decoupling | `cell54.py`, `cell55.py` | `cell54.out`, `cell55.out` |
| Theorem 6.2 | Odd-sector resolvent identity $Mu = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1}\psi$ | `cell54.py`, `cell55.py` | `cell54.out`, `cell55.out` |
| Proposition 6.3 & Corollary 6.3.1 | Exact $K$ commutator, Parseval norm identity & mode projections | `cell59.py`, `cell61.py` | `cell59.out`, `cell61.out` |
| Theorem 7.1 | Exact Archimedean resolvent first-jet identity | `cell54.py`, `cell55.py` | `cell54.out`, `cell55.out` |
| Theorem 7.2 | Denominator cancellation $(E_k - \lambda)$ in odd-even coupling | `cell55.py` | `cell55.out` |
| Theorem 7.3 & Corollary 7.3.1 | Exact $K^2$ commutator, even mode projections & $D_1/D_0$ representation | `cell59.py`, `cell61.py` | `cell59.out`, `cell61.out` |

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
11. Research Record / Connes–CvS Series, *The Dirichlet Continuum Limit, Barrier Mechanics, and Asymptotic Weil Positivity in the Connes–van Suijlekom Galerkin Truncation*, Companion Paper (Paper 4B), GitHub: `nrensen/connes-cvs-` (2026).
