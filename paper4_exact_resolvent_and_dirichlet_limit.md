# An Exact Resolvent and Commutator Toolkit for the Truncated Connes–van Suijlekom Weil Quadratic Form

**Authors:** Research Record / Connes–CvS Investigation Series  
**Date:** September 2026  
**Software & Reproducibility Suite:** `https://github.com/akivag613/connes-cvs-` (mirror: `nrensen/connes-cvs-`)  
**Status:** Standalone Manuscript / The Rigorous Toolkit (Companion to Paper 4B: *The Research Programme*)

---

### Abstract

The truncated Weil quadratic form developed by Connes–van Suijlekom (2025) and Connes–Consani–Moscovici (2026) discretizes the explicit formula of prime number theory using a finite-rank Galerkin projection of band $N$ on a logarithmic scaling interval $[0, L] = [0, \log c]$. The omitted Archimedean tail of this truncation has historically been treated as a difficult oscillatory numerical integration problem or as an empirical asymptotic inverse-power expansion.

In this paper, we establish the **exact algebraic solution** to the finite-$N$ Archimedean kernel and the associated commutator algebra, proving unconditionally and independently of numerical quadrature:

1. **Exact Rational Resolvent & Operator Identity (Theorem):** Starting from the four-term analytic reduction of the Archimedean Volterra integral, we prove algebraically that the reduced Fourier kernel $R_v(r) = K_{\mathrm{Fourier}}(v, r, L) / (1 - \cos(rL))$ is identically equal to the squared Cauchy resolvent:
   $$R_v(r) \equiv \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2} \sum_{m=1}^{N} \frac{r v_m}{r^2 - a_m^2} \right]^2, \qquad a_m = \frac{2\pi m}{L},$$
   on the punctured complex plane $\mathbb{C} \setminus \{0, \pm a_1, \dots, \pm a_N\}$. Furthermore, introducing the Neumann Laplacian $\mathcal{L} = -d^2/dt^2$ on $[0, L]$ with $T'(0) = T'(L) = 0$, the rational generating function $D(z)$ is the boundary evaluation of the operator resolvent:
   $$D(z) \equiv \big[(I + z\mathcal{L})^{-1} T_v\big](0) = \int_0^\infty e^{-s} \big[ e^{-sz\mathcal{L}} T_v \big](0) \, ds \quad (\operatorname{Re} z > 0),$$
   whose Taylor expansion around $z = 0$ reproduces the entire endpoint-jet hierarchy $D_k = T_v^{(2k)}(0)$.
2. **Universal Fourier Factorization and Unconditional Pointwise Positivity (Theorem):** The entire Fourier-side amplitude $\Phi_v(r)$ factors directly in terms of the boundary Neumann resolvent evaluated at the inverted spectral variable $z = -1/r^2$:
   $$\Phi_v(r) \equiv \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D\left(-\frac{1}{r^2}\right), \qquad K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 = \frac{4}{L} \frac{\sin^2(rL/2)}{r^2} D\left(-\frac{1}{r^2}\right)^2 \ge 0,$$
   proving algebraically and unconditionally that the Fourier-side Archimedean kernel $K_{\mathrm{Fourier}}(v, r, L)$ is pointwise non-negative for all real $r$ and all real coefficient vectors $v \in \mathbb{R}^{N+1}$.
3. **Spectral Lattice Orthogonality (Theorem):** At the lattice nodes $r = a_m$, the apparent poles cancel cleanly against the envelope zeros via removable singularities, yielding the exact sampling identity:
   $$K_{\mathrm{Fourier}}(v, 0, L) = L v_0^2 = L u_0^2, \qquad K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 = L u_m^2 \quad (m = 1, \dots, N),$$
   uncoupling the kernel into the squared Fourier coefficients.
4. **Exact Archimedean Cauchy Transform and Closed-Form Pole Decomposition (Theorems):** We evaluate the continuous Archimedean Cauchy transform $J(q) = \frac{1}{\pi}\int_0^\infty \frac{2q}{q^2 + r^2} K_{\mathrm{Fourier}}(v, r, L) \, dr$ in exact closed algebraic form via contour integration in the complex frequency plane, isolating the origin residue ($2v_0^2/q$), the discrete lattice poles ($2qv_m^2/(q^2+a_m^2)$), and the imaginary pole at $z = iq$. We establish the spatial Laplace duality $J(q) \equiv \int_0^L K_v^{\mathrm{phys}}(y) e^{-qy} dy$. Combined with the Weierstrass partial fraction expansion of the digamma function, this expresses the continuous Archimedean quadratic form $\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty [ \frac{\|v\|_2^2}{n+1} - J(q_n) ]$ (with $C_{\mathrm{arch}} = -\gamma - \log \pi$) as an unconditionally convergent algebraic series with fast $\mathcal{O}(n^{-2})$ absolute convergence, eliminating the need for numerical quadrature.
5. **Finite-$T$ Archimedean Cutoff Defect and Endpoint-Jet Representation (Theorem):** In finite-rank Galerkin implementations where the continuous Archimedean density is evaluated with a finite integration cutoff $T > 0$, we prove the exact identity $v^T Q_{\mathrm{arch}}^{(T)} v \equiv \frac{1}{\pi}\int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$. Consequently, the discrepancy between the Galerkin eigenvalue $\lambda_N$ and the continuous tripartite functional $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N)$ is identically the Archimedean cutoff tail $\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) \equiv -\frac{1}{\pi}\int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr$. For $T > a_N = 2\pi N/L$, this tail admits an exact Laurent series $\sum_{k=0}^\infty A_k(N) \mathcal{J}_k(T, L)$ in terms of endpoint derivatives $D_j = T_v^{(2j)}(0)$ and universal moment integrals $\mathcal{J}_k(T, L)$, converging geometrically with step ratio $(a_N/T)^2$. This proves that the observed $10^{-43}$ numerical residual at $N = 24$ is 100% continuous cutoff leakage rather than a finite-rank subspace projection defect.
6. **Exact Commutator Algebra and Parity Factorization (Theorem):** For the coordinate operator $M = \operatorname{diag}(n)$ and the Galerkin matrix $Q$, the commutator $[M^k, Q]$ has rank at most $2k$. Parity reflection decouples the system into even and odd sectors, yielding the exact odd-sector resolvent identity:
   $$M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi, \qquad B_1 = -D_0 \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$
7. **Exact First-Jet Resolvent Identity and Small-Denominator Cancellation (Theorems):** We prove algebraically that the first-jet cancellation ratio is identically the relative first correction of the large-$r$ Archimedean resolvent:
   $$\frac{D_1}{D_0} \equiv -\frac{1}{2} \frac{A_1}{A_0} = -\frac{D'(0)}{D(0)} = -\kappa^2 \frac{F'(0)}{F(0)} = -\frac{\sqrt{2}\displaystyle\sum_{m=1}^N a_m^2 v_m}{v_0 + \sqrt{2}\displaystyle\sum_{m=1}^N v_m}.$$
   In terms of the discrete signed atomic measure $\mu_N = v_0 \delta_0 + \sqrt{2} \sum_{m=1}^N v_m \delta_{a_m^2}$, this ratio equals $-\frac{\int_0^\infty x \, d\mu_N(x)}{\int_0^\infty d\mu_N(x)}$ as an algebraic quotient of finite sums. Furthermore, we establish the exact small-denominator cancellation theorem: in the spectral expansion of the resolvent $(Q_{\mathrm{even}} - \lambda I)^\dagger$ coupled to the odd arithmetic energy, the apparent denominator $(E_k - \lambda)$ for any even eigenstate $k$ cancels identically against the odd resolvent difference in the numerator via the first resolvent identity, proving algebraically that the resolvent coupling is non-singular.

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

### The Problem of the Archimedean Tail and the Role of the Toolkit

In finite implementations, numerical quadrature of the Archimedean integral suffered from severe oscillatory cancellation errors, while asymptotic inverse-power expansions $\sum A_k / r^{2k+2}$ presented non-trivial remainder bounds. Furthermore, whether $K_{\mathrm{Fourier}}(v, r, L)$ possessed sign-oscillations that could induce negative eigenvalues at large frequencies remained an open question.

The objective of this manuscript is to provide a **100% rigorous, pure-mathematics toolkit** that resolves all finite-$N$ algebraic and operator questions unconditionally. We establish:
- The exact squared rational resolvent identity for $K_{\mathrm{Fourier}}$, proving pointwise non-negativity globally.
- The operator representation of the generating function as a Neumann resolvent.
- The closed-form contour evaluation of the Archimedean Cauchy transform $J(q)$ and the unconditionally convergent Weierstrass pole series for $\mathcal{Q}_{\mathrm{arch}}(v)$, completely eliminating numerical quadrature.
- The finite-$T$ Archimedean cutoff defect identity and geometric endpoint-jet Laurent series, resolving the exact relationship between finite Galerkin truncations and continuous Weil functionals.
- The exact rank-$2k$ commutator algebra of the Galerkin matrix, proving that barrier suppression factors cancel identically in the first-jet ratio.
- The exact cancellation of small bound-state denominators and universal semigroup squeezing bounds.

All exploratory and asymptotic questions concerning the infinite-dimensional limit $N \to \infty$ (such as solitary wave profiles, semiclassical WKB barrier tunneling, and formal Wiener–Hopf continuum scaling) are systematically investigated in the companion paper, *The Dirichlet Continuum Limit, Barrier Mechanics, and Asymptotic Weil Positivity in the Connes–van Suijlekom Galerkin Truncation* [Paper 4B].

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

By direct integration of the Volterra convolution against $\cos(r L \omega)$, the boundary terms at the endpoint $\omega = 1$ ($x = L$) factor out cleanly, isolating the common oscillatory factor $1 - \cos(rL)$:

$$K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r),$$

where $R_v(r)$ is a rational function of $r^2$. Defining the fundamental lattice scale:

$$\kappa = \frac{2\pi}{L}, \qquad a_m = \kappa m = \frac{2\pi m}{L},$$

and the elementary mode integrals:

$$S_{\mathrm{bar}}(m) = \frac{a_m}{a_m^2 - r^2}, \qquad C_{\mathrm{bar}}(m) = \frac{r^2 + a_m^2}{L (r^2 - a_m^2)^2},$$

the rational kernel was originally defined by the four-part interaction sum:

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
This matches the cross-terms of the square identically. Since all four blocks match identically, the algebraic identity $R_v(r) \equiv \frac{1}{r^2} A(1/r^2)$ is exact. $\blacksquare$

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

### Proof of Theorem 5.1
Decompose the kernel into $K_{\mathrm{Fourier}}(v, r, L) = f_1(r) + f_2(r)$ where $f_1(z) = H(z)(1 - e^{i z L})$, $f_2(z) = H(z)(1 - e^{-i z L})$, and $H(z) = \frac{1}{2} R_v(z) = \frac{1}{L z^2} [D(-1/z^2)]^2$. Extending over $(-\infty, \infty)$:

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

Summing all contributions proves the identity. $\blacksquare$

### Theorem 5.2 (Spatial Laplace Duality)
*The Archimedean Cauchy transform $J(q)$ is identically the spatial Laplace transform of the physical Volterra kernel:*

$$J(q) \equiv \int_0^L K_v^{\mathrm{phys}}(y) e^{-q y} \, dy \qquad \forall q > 0.$$

*Proof.* By the spatial cosine representation of the Volterra convolution $K_{\mathrm{Fourier}}(v, r, L) = \int_0^L K_v^{\mathrm{phys}}(y) \cos(ry) dy$. Substituting into the Cauchy transform and interchanging integrals via Fubini's theorem (valid since $K_v^{\mathrm{phys}}$ is smooth on $[0, L]$):
$$J(q) = \frac{1}{\pi} \int_0^\infty \frac{2q}{q^2 + r^2} \left[ \int_0^L K_v^{\mathrm{phys}}(y) \cos(ry) \, dy \right] dr = \int_0^L K_v^{\mathrm{phys}}(y) \left[ \frac{2q}{\pi} \int_0^\infty \frac{\cos(ry)}{q^2 + r^2} \, dr \right] dy.$$
The inner integral is the standard Laplace contour integral $\frac{2q}{\pi} \frac{\pi}{2q} e^{-qy} = e^{-qy}$, establishing $J(q) = \int_0^L K_v^{\mathrm{phys}}(y) e^{-qy} dy$. $\blacksquare$

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

2. **Dominated Convergence and Interchange of Limits:**
   To justify interchanging summation and integration $\sum_{n=0}^\infty \int_0^\infty = \int_0^\infty \sum_{n=0}^\infty$, consider the partial sums $S_M(r) = C_{\mathrm{arch}} + \sum_{n=0}^M w_n(r)$.
   Notice that for each $n \ge 0$ and all $r \ge 0$:
   $$w_n(r) = \frac{1}{n+1} - \frac{2q_n}{q_n^2 + r^2} \ge \frac{1}{n+1} - \frac{2}{q_n} = \frac{1}{n+1} - \frac{1}{n + 1/4} = -\frac{3/4}{(n+1)(n + 1/4)}.$$
   Thus the negative part of each term is uniformly bounded by a convergent numerical series independent of $r$:
   $$\sum_{n=0}^\infty \max(0, -w_n(r)) \le \sum_{n=0}^\infty \frac{3/4}{(n+1)(n + 1/4)} = \psi(1) - \psi(1/4) = \frac{\pi}{2} + 3\log 2 < \infty.$$
   Because the partial sums $S_M(r)$ converge pointwise to $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$, which satisfies the logarithmic bound $|h_+(r)| \le C_1 + \frac{1}{2} \log(1 + r^2)$ for all $r \ge 0$, the partial sums satisfy the uniform $M$-independent upper envelope:
   $$|S_M(r)| \le |h_+(r)| + 2 \sum_{n=0}^\infty \max(0, -w_n(r)) \le C_0 + \log(1 + r) \equiv g(r).$$
   Since $K_{\mathrm{Fourier}}(v, r, L) \le \frac{C_K}{1 + r^2}$, the dominant function satisfies $\int_0^\infty g(r) K_{\mathrm{Fourier}}(r) dr < \infty$. By Lebesgue's Dominated Convergence Theorem, term-by-term integration is unconditionally valid:
   $$\mathcal{Q}_{\mathrm{arch}}(v) = \frac{1}{\pi} \int_0^\infty \lim_{M\to\infty} S_M(r) K_{\mathrm{Fourier}}(r) \, dr = \lim_{M\to\infty} \sum_{n=0}^M I_n = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty \left[ \frac{\|v\|_2^2}{n+1} - J(q_n) \right].$$

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
*By summing the harmonic difference and the discrete lattice pole series in Corollary 5.3 analytically in terms of the digamma function $\psi(z)$, the Archimedean quadratic form reduces to the exact discrete lattice sum plus an exponentially suppressed boundary leakage term:*

$$\mathcal{Q}_{\mathrm{arch}}(v) = h_+(0) v_0^2 + \sum_{m=1}^N v_m^2 h_+(a_m) + \sum_{n=0}^\infty \frac{2(1 - e^{-q_n L})}{L q_n^2} \left[ v_0 + \sqrt{2}\sum_{m=1}^N \frac{q_n^2 v_m}{q_n^2 + a_m^2} \right]^2,$$

*where $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$.*

*Proof.* In the summand of Corollary 5.3, write $\frac{\|v\|_2^2}{n+1} - J(q_n) = (\frac{\|v\|_2^2}{n+1} - \frac{2\|v\|_2^2}{q_n}) + (\frac{2\|v\|_2^2}{q_n} - J(q_n))$. Summing the first part gives $\|v\|_2^2 \sum_{n=0}^\infty (\frac{1}{n+1} - \frac{1}{n+1/4}) = $\|v\|_2^2 (\psi(1/4) + \gamma)$, which combined with $C_{\mathrm{arch}}\|v\|_2^2 = (-\gamma - \log\pi)\|v\|_2^2$ yields $(\psi(1/4) - \log\pi)\|v\|_2^2 = h_+(0)\|v\|_2^2$. For the second part, using the identity $\frac{2\|v\|_2^2}{q_n} - J(q_n) = \sum_{m=1}^N \frac{2 a_m^2 v_m^2}{q_n(q_n^2 + a_m^2)} + \text{boundary term}$, summing over $n$ using $\sum_{n=0}^\infty \frac{2 a_m^2}{q_n(q_n^2 + a_m^2)} = \operatorname{Re}\psi(1/4 + i a_m/2) - \psi(1/4) = h_+(a_m) - h_+(0)$ yields $\sum_{m=1}^N v_m^2 [h_+(a_m) - h_+(0)]$. Adding to $h_+(0)\|v\|_2^2$ cancels the $-h_+(0)\sum v_m^2$ term, leaving $h_+(0) v_0^2 + \sum_{m=1}^N v_m^2 h_+(a_m)$ plus the boundary sum. $\blacksquare$

### Theorem 5.5 (Finite-$T$ Archimedean Cutoff Defect and Endpoint-Jet Representation)
*Let $Q_{\mathrm{arch}}^{(T)}$ be the finite-cutoff Galerkin Archimedean matrix of band $N$ defined via divided differences of the cutoff test function $\psi_{\mathrm{arch}}^{(T)}(x) = \frac{1}{\pi^2}\int_0^T h_+(r) \operatorname{Re}\hat{S}_x(r) \, dr$ for $T > 0$.*

1. **Exact Divided-Difference Integral Identity:**
   *For any canonical coefficient vector $v \in \mathbb{R}^{N+1}$ corresponding to the even full vector $u \in \mathbb{R}^{2N+1}$, the matrix quadratic form is identically equal to the continuous Fourier integral truncated at $T$:*
   $$v^T Q_{\mathrm{arch}}^{(T)} v \equiv u^T Q_{\mathrm{arch}}^{(T)} u \equiv \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr.$$

2. **Exact Archimedean Cutoff Defect:**
   *Let $Q^{(T)} = Q_{\mathrm{arch}}^{(T)} + Q_{\mathrm{pole}} + Q_{\mathrm{prime}}$ be the full Galerkin matrix with cutoff $T$, and let $u_N$ be an even ground-state eigenvector with eigenvalue $\lambda_N = u_N^T Q^{(T)} u_N$. Let $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{(\infty)}(v_N)$ be the exact continuous tripartite functional evaluated at $T = \infty$ (via Corollary 5.4). Then the discrepancy between the Galerkin eigenvalue and the continuous functional is identically the Archimedean cutoff tail:*
   $$\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) \equiv -\frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr \equiv -\delta_T^{\mathrm{tail}}(v_N).$$
   *In particular, the residual is 100% Archimedean cutoff tail leakage and does not originate from a finite-rank subspace projection defect.*

3. **Taylor Endpoint-Jet Laurent Decomposition:**
   *For any cutoff $T > a_N = \frac{2\pi N}{L}$, the rational resolvent $R_v(r) = \frac{2}{L} \left[ \frac{v_0}{r} + \sqrt{2}\sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} \right]^2$ converges uniformly on $[T, \infty)$ as the Laurent series:*
   $$R_v(r) = \sum_{k=0}^\infty \frac{A_k(N)}{r^{2k+2}}, \qquad A_k(N) = \frac{2}{L} (-1)^k \sum_{j=0}^k D_j D_{k-j},$$
   *where $D_j = T_v^{(2j)}(0)$ are the physical endpoint derivatives (Theorem 3.2). Consequently, the tail defect admits the exact series representation:*
   $$\delta_T^{\mathrm{tail}}(v_N) = \sum_{k=0}^\infty A_k(N) \mathcal{J}_k(T, L),$$
   *where the universal tail moments:*
   $$\mathcal{J}_k(T, L) \equiv \frac{1}{\pi} \int_T^\infty \frac{h_+(r) (1 - \cos(rL))}{r^{2k+2}} \, dr$$
   *are independent of the eigenvector $v_N$.*

4. **Explicit Geometric Remainder Bound:**
   *For any integration cutoff $T > a_N = \frac{2\pi N}{L}$ and any truncation order $K \ge 0$, the tail truncation remainder is unconditionally bounded by:*
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{C(v, N, L)}{T^{2K+3}} \frac{1}{1 - (a_N/T)^2},$$
   *where the constant depends explicitly on the physical wave norm $\|v\|_1 = \sum_{m=0}^N |v_m|$, frequency bandwidth $a_N$, and interval length $L$:*
   $$C(v, N, L) \equiv \frac{8 (\log T + 1)}{\pi L} \|v\|_1^2 a_N^{2K+2}.$$
   *Equivalently, expressed in dimensionless scaling form:*
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{C_0(v, L, T)}{T} \left( \frac{a_N}{T} \right)^{2K+2} \frac{1}{1 - (a_N/T)^2},$$
   *with $C_0(v, L, T) = \frac{8 (\log T + 1)}{\pi L} \|v\|_1^2$. Consequently, the Laurent series converges unconditionally and geometrically with base ratio $(a_N/T)^2 < 1$. Away from non-generic algebraic cancellations among the endpoint derivatives, successive non-vanishing terms scale with this characteristic geometric step ratio.*

*Proof.*
1. **Divided-Difference Integral Identity:**
   In the Connes–van Suijlekom Galerkin discretization [2, 3], the Archimedean coordinate function $\psi_{\mathrm{arch}}^{(T)}(x)$ is defined by integrating the smooth test-function transform $S_x(r) = \int_{-L/2}^{L/2} e^{i r y} \frac{\cos(2\pi x y / L)}{L/2} dy$ over $[-T, T]$ against $h_+(r)$:
   $$\psi_{\mathrm{arch}}^{(T)}(x) = \frac{1}{2\pi^2} \int_{-T}^T h_+(r) S_x(r) \, dr.$$
   The matrix entries $(Q_{\mathrm{arch}}^{(T)})_{mn}$ are defined by divided differences for $m \ne n$ and derivatives for $m = n$. Interchanging the finite sum over $m, n \in \{-N, \dots, N\}$ with the compact Riemann integral on $[-T, T]$:
   $$u^T Q_{\mathrm{arch}}^{(T)} u = \sum_{m,n=-N}^N u_m u_n (Q_{\mathrm{arch}}^{(T)})_{mn} = \frac{1}{2\pi^2} \int_{-T}^T h_+(r) \left[ \sum_{m,n=-N}^N u_m u_n \frac{S_m(r) - S_n(r)}{m - n} \right] dr.$$
   For any even coefficient vector $u$ ($u_{-m} = u_m$), the discrete bilinear form inside the integral evaluates via the Fourier-side kernel identity (Theorem 4.1 & Theorem 4.3) identically to $2\pi K_{\mathrm{Fourier}}(v, r, L)$. Since $h_+(r)$ and $K_{\mathrm{Fourier}}(v, r, L)$ are even functions of $r$:
   $$\frac{1}{2\pi^2} \int_{-T}^T h_+(r) \big[ 2\pi K_{\mathrm{Fourier}}(v, r, L) \big] dr = \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr,$$
   establishing identity (1) identically for every $T > 0$.

2. **Exact Archimedean Cutoff Defect:**
   The full discrete Galerkin matrix is $Q^{(T)} = Q_{\mathrm{arch}}^{(T)} + Q_{\mathrm{pole}} + Q_{\mathrm{prime}}$. Since the pole and prime matrices have no frequency cutoff and match their continuous functional counterparts identically ($u_N^T Q_{\mathrm{pole}} u_N = \mathcal{Q}_{\mathrm{pole}}(v_N)$ and $u_N^T Q_{\mathrm{prime}} u_N = \mathcal{Q}_{\mathrm{prime}}(v_N)$), subtracting the exact continuous functional $\mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) = \mathcal{Q}_{\mathrm{pole}}(v_N) + \mathcal{Q}_{\mathrm{prime}}(v_N) + \mathcal{Q}_{\mathrm{arch}}^{(\infty)}(v_N)$ gives:
   $$\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v_N) = u_N^T Q_{\mathrm{arch}}^{(T)} u_N - \mathcal{Q}_{\mathrm{arch}}^{(\infty)}(v_N) = \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr - \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr.$$
   Combining the integrals yields $-\frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v_N, r, L) \, dr = -\delta_T^{\mathrm{tail}}(v_N)$.

3. **Taylor Endpoint-Jet Laurent Decomposition:**
   For $r \ge T > a_N$, we have $a_m / r \le a_N / T < 1$ for all $m \in \{1, \dots, N\}$. Expanding each term in the rational resolvent (Theorem 3.1):
   $$\frac{r v_m}{r^2 - a_m^2} = \frac{v_m}{r} \frac{1}{1 - a_m^2 / r^2} = \sum_{j=0}^\infty \frac{a_m^{2j} v_m}{r^{2j+1}}.$$
   Summing over $m$ and adding $v_0/r$:
   $$\frac{v_0}{r} + \sqrt{2} \sum_{m=1}^N \frac{r v_m}{r^2 - a_m^2} = \frac{1}{r} \left[ v_0 + \sqrt{2}\sum_{m=1}^N v_m + \sum_{j=1}^\infty \frac{\sqrt{2}\sum_{m=1}^N a_m^{2j} v_m}{r^{2j}} \right] = \sum_{j=0}^\infty \frac{(-1)^j D_j}{r^{2j+1}},$$
   where $D_0 = v_0 + \sqrt{2}\sum_{m=1}^N v_m$ and $D_j = (-1)^j \sqrt{2}\sum_{m=1}^N a_m^{2j} v_m = T_v^{(2j)}(0)$ are the physical endpoint derivatives (Theorem 3.2).
   Squaring this series and taking the Cauchy product yields:
   $$R_v(r) = \frac{2}{L} \left[ \sum_{j=0}^\infty \frac{(-1)^j D_j}{r^{2j+1}} \right]^2 = \sum_{k=0}^\infty \frac{A_k(N)}{r^{2k+2}}, \qquad A_k(N) = \frac{2}{L} (-1)^k \sum_{j=0}^k D_j D_{k-j}.$$
   Since $K_{\mathrm{Fourier}}(v, r, L) = (1 - \cos(rL)) R_v(r)$ (Theorem 4.1), substituting into the tail integral gives $\delta_T^{\mathrm{tail}}(v_N) = \sum_{k=0}^\infty A_k(N) \mathcal{J}_k(T, L)$. Term-by-term integration is justified by the dominated convergence theorem, as the series converges uniformly on $[T, \infty)$ and $h_+(r)(1 - \cos(rL)) = \mathcal{O}(\log r)$.

4. **Explicit Geometric Remainder Bound:**
   To establish the explicit remainder bound with full parameter accounting, let $r \ge T > a_N$. Since $a_m \le a_N$ for all $m$, the endpoint derivatives satisfy the uniform upper bound:
   $$|D_j| \le \sqrt{2} \sum_{m=1}^N a_m^{2j} |v_m| \le \sqrt{2} a_N^{2j} \sum_{m=0}^N |v_m| = \sqrt{2} \|v\|_1 a_N^{2j} \qquad (\forall j \ge 0).$$
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
   $$\left| \delta_T^{\mathrm{tail}}(v_N) - \sum_{k=0}^K A_k(N) \mathcal{J}_k(T, L) \right| \le \frac{C(v, N, L)}{T^{2K+3}} \frac{1}{1 - (a_N/T)^2} = \frac{C_0(v, L, T)}{T} \left( \frac{a_N}{T} \right)^{2K+2} \frac{1}{1 - (a_N/T)^2},$$
   where $C(v, N, L) = \frac{8(\log T + 1)}{\pi L} \|v\|_1^2 a_N^{2K+2}$ and $C_0(v, L, T) = \frac{8(\log T + 1)}{\pi L} \|v\|_1^2$. This completes the rigorous proof. $\blacksquare$

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
*Let $u$ be the even ground-state eigenvector satisfying $Q u = \lambda u$ with $\lambda = \lambda_{\min}(N) > 0$.*
1. **Odd-Sector Resolvent Identity:**
   *The $k = 1$ commutator applied to $u$ yields:*
   $$M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi.$$
2. **Arithmetic Energy Proportionality:**
   *Consequently, the arithmetic moment scalar $B_1$ is strictly proportional to $D_0$:*
   $$B_1 = -D_0 \mathcal{E}_{\mathrm{arith}}, \qquad \mathcal{E}_{\mathrm{arith}} \equiv \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle > 0.$$

*Proof.* Setting $k = 1$ in Theorem 6.1 gives $[M, Q] = p e^T - e p^T$. Applying to $u$ and using $p^T u = 0$ yields $[M, Q] u = D_0 \psi$. Since $[M, Q] u = \lambda M u - Q M u$, we obtain $(Q - \lambda I)(M u) = -D_0 \psi$. Because $M u$ and $\psi$ are odd, inverting on $\mathcal{H}_{\mathrm{odd}}$ establishes $M u = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi$. Taking the inner product with $\psi$ yields $B_1 = \psi^T M u = -D_0 \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle = -D_0 \mathcal{E}_{\mathrm{arith}}$. $\blacksquare$

---

## 7. The First-Jet Ratio $D_1/D_0$ and Exact Resolvent Small-Denominator Cancellation

We now establish the exact algebraic structure of the first-jet cancellation ratio $D_1 / D_0$ and the non-singular cancellation of small denominators in the resolvent spectral expansion.

### Theorem 7.1 (Exact Archimedean Resolvent Jet Representation)
*The first-jet cancellation ratio $D_1 / D_0$ is identically equal to the relative first correction of the large-$r$ Archimedean resolvent:*

$$\frac{D_1}{D_0} \equiv -\frac{1}{2} \frac{A_1}{A_0} = -\frac{D'(0)}{D(0)} = -\kappa^2 \frac{F'(0)}{F(0)} = -\frac{\sqrt{2}\displaystyle\sum_{m=1}^N a_m^2 v_m}{v_0 + \sqrt{2}\displaystyle\sum_{m=1}^N v_m},$$

*where $D(w) \equiv \big[(I + w\mathcal{L})^{-1} T_v\big](0) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 w}$ is the operator resolvent generating function, and $F(z) = e^T (I - z M^2)^{-1} v = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - z m^2}$ is the discrete mode generating function.*

*Proof.* From Theorem 3.1, the large-$r$ expansion of the reduced Fourier kernel is:

$$R_v(r) = \frac{A_0}{r^2} + \frac{A_1}{r^4} + \mathcal{O}(r^{-6}),$$

where $A_0 = \frac{2}{L} D_0^2$ and $A_1 = -\frac{4}{L} D_0 D_1$. Taking the ratio gives $-A_1 / (2 A_0) = -(-\frac{4}{L} D_0 D_1) / (\frac{4}{L} D_0^2) = D_1 / D_0$. 

In terms of the operator resolvent generating function $D(w) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + a_m^2 w}$, Taylor expansion around $w = 0$ gives $D(w) = D_0 - D_1 w + \mathcal{O}(w^2)$, whence $D(0) = D_0 = v_0 + \sqrt{2}\sum_{m=1}^N v_m$ and $D'(0) = -D_1 = -\sqrt{2}\sum_{m=1}^N a_m^2 v_m$. Thus:

$$-\frac{D'(0)}{D(0)} = \frac{D_1}{D_0} = -\frac{\sqrt{2}\sum_{m=1}^N a_m^2 v_m}{v_0 + \sqrt{2}\sum_{m=1}^N v_m}.$$

In terms of the unscaled coordinate generating function $F(z) = e^T (I - z M^2)^{-1} v = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 - z m^2}$, expanding around $z = 0$ gives $F(z) = D_0 + (\sqrt{2}\sum_{m=1}^N m^2 v_m) z + \mathcal{O}(z^2) = D_0 - \frac{D_1}{\kappa^2} z + \mathcal{O}(z^2)$. Differentiating at $z = 0$ yields $F'(0) = -D_1 / \kappa^2$, whence $F'(0) / F(0) = -D_1 / (\kappa^2 D_0)$, which gives $\frac{D_1}{D_0} = -\kappa^2 \frac{F'(0)}{F(0)}$. $\blacksquare$

*Remark 7.1.1 (Algebraic Ratio vs. Measure Expectation).*
In terms of the signed atomic measure $\mu_N = v_0 \delta_0 + \sqrt{2}\sum_{m=1}^N v_m \delta_{a_m^2}$ supported on the discrete Neumann Laplacian spectrum $\{0, a_1^2, \dots, a_N^2\}$, this ratio evaluates formally as $-\frac{\int_0^\infty x \, d\mu_N(x)}{\int_0^\infty d\mu_N(x)}$. Because the Fourier coefficients $v_m$ alternate and take negative values, $\mu_N$ is a signed measure rather than a positive probability distribution. The identity is a strictly algebraic relation between finite sums, not a probabilistic expectation.

---

### Theorem 7.2 (Exact Small-Denominator Cancellation in the Odd-Even Spectral Coupling)
*Let $Q_{\mathrm{even}}$ and $Q_{\mathrm{odd}}$ denote the restrictions of the Connes–van Suijlekom Galerkin matrix $Q$ to the even and odd subspaces of $\mathbb{C}^{2N+1}$. For any even eigenvector $u^{(k)}$ of $Q_{\mathrm{even}}$ with eigenvalue $E_k$ satisfying the non-resonance condition $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$ and $\lambda \notin \operatorname{spec}(Q_{\mathrm{odd}})$, let $D_0^{(k)} = e^T u^{(k)}$ and define the arithmetic energy vector $\mathbf{w}(\lambda) = b + \mathcal{E}_{\mathrm{arith}}(\lambda) e$, where $b_n = n \psi(n)$ and $\mathcal{E}_{\mathrm{arith}}(\lambda) = \langle \psi, (Q_{\mathrm{odd}} - \lambda I)^{-1}\psi \rangle$.*

*Then in the spectral expansion of $(Q_{\mathrm{even}} - \lambda I)^\dagger \mathbf{w}(\lambda)$, the numerator for mode $k$ factors identically as:*

$$u^{(k)T} \mathbf{w}(\lambda) = -(E_k - \lambda) D_0^{(k)} \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$

*Consequently, for every $E_k \ne \lambda$ with $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$, the factor $(E_k - \lambda)$ cancels identically, yielding the exact non-singular algebraic summand:*

$$\frac{(e^T u^{(k)}) (u^{(k)T} \mathbf{w}(\lambda))}{E_k - \lambda} \equiv - [D_0^{(k)}]^2 \cdot \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$

*Proof.* Expanding $(Q_{\mathrm{even}} - \lambda I)^\dagger$ along the orthonormal eigenbasis $\{u^{(k)}\}_{k \ge 1, \text{ even}}$ gives numerators $(e^T u^{(k)}) (u^{(k)T} \mathbf{w}(\lambda)) = D_0^{(k)} [B_1^{(k)} + \mathcal{E}_{\mathrm{arith}}(\lambda) D_0^{(k)}]$, where $B_1^{(k)} = \psi^T M u^{(k)}$. Since $[M, Q] u^{(k)} = D_0^{(k)} \psi$ by Theorem 6.1 and $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$, the operator $(Q_{\mathrm{odd}} - E_k I)$ is boundedly invertible on $\mathcal{H}_{\mathrm{odd}}$, yielding $M u^{(k)} = -D_0^{(k)} (Q_{\mathrm{odd}} - E_k I)^{-1} \psi$, whence $B_1^{(k)} = \psi^T M u^{(k)} = -D_0^{(k)} \mathcal{E}_{\mathrm{arith}}(E_k)$. The bracket is therefore:

$$B_1^{(k)} + \mathcal{E}_{\mathrm{arith}}(\lambda) D_0^{(k)} = -D_0^{(k)} \big[ \mathcal{E}_{\mathrm{arith}}(E_k) - \mathcal{E}_{\mathrm{arith}}(\lambda) \big].$$

By the first resolvent identity for the odd operator $Q_{\mathrm{odd}}$ (valid since $E_k, \lambda \notin \operatorname{spec}(Q_{\mathrm{odd}})$):

$$(Q_{\mathrm{odd}} - E_k I)^{-1} - (Q_{\mathrm{odd}} - \lambda I)^{-1} = (E_k - \lambda) (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1}.$$

Taking the inner product with $\psi$:

$$\mathcal{E}_{\mathrm{arith}}(E_k) - \mathcal{E}_{\mathrm{arith}}(\lambda) = (E_k - \lambda) \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle.$$

Substituting this into the numerator yields $u^{(k)T} \mathbf{w}(\lambda) = -(E_k - \lambda) D_0^{(k)} \langle \psi, (Q_{\mathrm{odd}} - E_k I)^{-1} (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi \rangle$. Dividing by $(E_k - \lambda)$ completes the algebraic cancellation. $\blacksquare$

*Remark 7.2.1 (Non-Resonance and Spectral Disjointness).*
The hypothesis $E_k \notin \operatorname{spec}(Q_{\mathrm{odd}})$ ensures that the odd resolvent $(Q_{\mathrm{odd}} - E_k I)^{-1}$ is well-defined. In all numerical Galerkin calculations across all tested dimensions and cutoffs ($N \in \{1, \dots, 24\}$, $c \in \{5, \dots, 17\}$), the even and odd spectra are strictly disjoint ($\operatorname{spec}(Q_{\mathrm{even}}) \cap \operatorname{spec}(Q_{\mathrm{odd}}) = \emptyset$). In particular, the lowest odd eigenvalue satisfies $E_{\mathrm{odd}, 1} \approx 0.05 \gg \lambda \sim 10^{-43}$, so no spectral degeneracies or resonances occur between the sectors.

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
| **Lattice Orthogonality** | $K_{\mathrm{Fourier}}(a_m) = \frac{L}{2} v_m^2$ | Theorem 4.3 | Samples squared Fourier coefficients orthogonally at lattice nodes. |
| **Cauchy Transform Identity** | Closed algebraic form for $J(q)$ | Theorem 5.1 | Evaluates Archimedean Cauchy transform in closed form via contour integration. |
| **Spatial Laplace Duality** | $J(q) = \int_0^L K_v^{\mathrm{phys}}(y) e^{-qy} dy$ | Theorem 5.2 | Duality between complex frequency Cauchy transform and physical Laplace transform. |
| **Weierstrass Pole Series** | $\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum [\frac{\|v\|_2^2}{n+1} - J(q_n)]$ | Corollary 5.3 | Unconditionally convergent algebraic series with fast $\mathcal{O}(n^{-2})$ convergence. |
| **Digamma Closed Form** | $\mathcal{Q}_{\mathrm{arch}}(v) = h_+(0) v_0^2 + \sum v_m^2 h_+(a_m) + \text{boundary}$ | Corollary 5.4 | Closed-form evaluation summing harmonic and lattice modes analytically. |
| **Finite-$T$ Cutoff Defect** | $\lambda_N - \mathcal{Q}_{\mathrm{total}}^{(\infty)} \equiv -\delta_T^{\mathrm{tail}} \equiv -\sum A_k \mathcal{J}_k(T, L)$ | Theorem 5.5 | Proves residual is 100% Archimedean cutoff leakage; geometric jet expansion. |
| **Rank-$2k$ Commutator** | $[M^k, Q] = \sum (M^j p)(M^{k-1-j} e)^T - (M^j e)(M^{k-1-j} p)^T$ | Theorem 6.1 | Determines exact commutator structure of Galerkin matrix. |
| **Odd-Sector Resolvent** | $Mu = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1} \psi$ | Theorem 6.2 | Solves first spectral moment and arithmetic energy uniquely in odd sector. |
| **Exact First-Jet Identity** | $D_1/D_0 \equiv -\frac{1}{2} A_1/A_0 \equiv -\kappa^2 F'/F$ | Theorem 7.1 | Relates endpoint derivative ratio directly to large-$r$ resolvent tail. |
| **Small-Denominator Cancellation** | $(E_k - \lambda)$ cancels identically in resolvent coupling | Theorem 7.2 | Proves algebraically that odd-even resolvent coupling is non-singular. |

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
| Theorem 4.3 | Lattice orthogonality $K_{\mathrm{Fourier}}(a_m) = \frac{L}{2} v_m^2$ | `cell32.py` | `cell32.out` |
| Theorem 5.1 & Theorem 5.2 | Exact Cauchy transform $J(q)$ & spatial Laplace duality | `cell56.py` | `cell56.out` |
| Corollary 5.3 & Corollary 5.4 | Weierstrass pole decomposition & exact digamma closed form | `cell56.py` | `cell56.out` |
| Theorem 5.5 | Finite-$T$ Archimedean cutoff defect & endpoint-jet reconstruction | `cell57.py` | `cell57.out` |
| Theorem 6.1 & Corollary 6.1.1 | Rank-$2k$ commutator algebra & strict parity decoupling | `cell54.py`, `cell55.py` | `cell54.out`, `cell55.out` |
| Theorem 6.2 | Odd-sector resolvent identity $Mu = -D_0 (Q_{\mathrm{odd}} - \lambda I)^{-1}\psi$ | `cell54.py`, `cell55.py` | `cell54.out`, `cell55.out` |
| Theorem 7.1 | Exact Archimedean resolvent first-jet identity | `cell54.py`, `cell55.py` | `cell54.out`, `cell55.out` |
| Theorem 7.2 | Small-denominator cancellation $(E_k - \lambda)$ in odd-even coupling | `cell55.py` | `cell55.out` |

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
