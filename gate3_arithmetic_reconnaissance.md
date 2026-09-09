# Gate 3 Arithmetic Reconnaissance: The Coordinate-Space Decomposition of $Q_{c, N}$ and the Weil Bridge

**Document Type:** Strategic Mathematical Working Note (Gate 3 Normalization & Structural Audit)  
**Date:** September 2026  
**Status:** Working Analytical Foundational Note (Finite-Subspace Identity Proved Exactly on $\mathcal{H}_{c, N}$; Continuum Extension Open)  
**Associated Manuscripts & Records:** [ROADMAP.md](file:///c:/data/github/connes-cvs-/ROADMAP.md) (Gate 3) | [paper4.md](file:///c:/data/github/connes-cvs-/paper4_exact_resolvent_and_dirichlet_limit.md) | [cell_history_map.md](file:///c:/data/github/connes-cvs-/cell_history_map.md)  

---

## Executive Summary & Calibrated Gate 3 Verdict

This working note conducts a rigorous structural and normalization audit of Gate 3 in the Connes–van Suijlekom (CvS) Galerkin framework:
$$\boxed{\textbf{What is the exact arithmetic content of the finite Galerkin matrix } Q_{c, N}\textbf{?}}$$

Specifically, it resolves whether the finite-rank Galerkin projection on the Fourier lattice $\Lambda_{\mathrm{Fourier}} = \{2\pi m / L\}$ introduces an irreversible aliasing distortion of André Weil's prime arithmetic support $\Lambda_{\mathrm{arith}} = \{\log p^k\}$, or whether the matrix represents a genuine projection of Weil's explicit quadratic functional.

### The Calibrated Verdict
$$\boxed{\textbf{GATE 3 FINITE-SUBSPACE IDENTITY: PROVED EXACTLY ON } \mathcal{H}_{c, N}}$$
$$\boxed{\textbf{GATE 2/3 DENSITY AND CONTINUUM EXTENSION: OPEN}}$$

Following complete derivations of all three components (prime, pole, and Archimedean), the structural findings are classified into:

### 1. Established Structural Theorems
- **Theorem A (Exact Prime-Power Autocorrelation Sampling):**  
  For every finite $N \ge 1$ and every test vector $v \in \mathbb{R}^{N+1}$, the prime block $Q_{\mathrm{prime}}^{(c)}$ evaluates the spatial autocorrelation of the wavepacket $f_v$ **pointwise and exactly** at each prime-power shift $y = \log p^k$:
  $$\langle v, Q_{\mathrm{prime}}^{(c)} v \rangle = -\sum_{p^k \le c} \frac{\log p}{p^{k/2}} K_v\left(1 - \frac{\log p^k}{L}\right) = -\frac{1}{\pi} \sum_{p^k \le c} \frac{\log p}{p^{k/2}} \widehat{g}_v\left(\frac{\log p^k}{2\pi}\right).$$
- **Theorem B (Zero Arithmetic Aliasing):**  
  The Fourier modes $2\pi m / L$ serve strictly as the coordinate basis for the test function $f_v$; arithmetic evaluation occurs on $\Lambda_{\mathrm{arith}}$ without spatial smearing or periodic aliasing. Fourier discretization does not compete with arithmetic localization.
- **Theorem C (Exact Archimedean Normalization & Local Regularization):**  
  The Weierstrass partial fraction series of the digamma Cauchy transforms collapses in coordinate space to the classical regularized Weil Archimedean distribution with the exact local constant $h_+(0) = \psi(1/4) - \log \pi$:
  $$\langle v, Q_{\mathrm{arch}, \infty} v \rangle = h_+(0) \|v\|_2^2 + \int_0^L \frac{K_v(1) - K_v(1 - y/L)}{2\sinh y} \, dy + K_v(1) \operatorname{artanh}(e^{-L}),$$
  matching $\frac{1}{2\pi}\int_{-\infty}^\infty h_+(r) g_v(r) dr$ identically with no undetermined constants.
- **Theorem D (Exact Pole Residue Identity):**  
  Tracing the normalizations from first principles eliminates the apparent factor-of-2 discrepancy, establishing:
  $$\langle v, Q_{\mathrm{pole}}^{(N)} v \rangle \equiv 2 g_v(i/2) \equiv g_v(i/2) + g_v(-i/2) \qquad \forall N \ge 1, \; \forall c > 1.$$
- **Theorem E (Subspace Identification):**  
  On the finite Galerkin subspace $\mathcal{H}_{c, N}$, the cutoff-free matrix equals the **prime-truncated Weil functional**:
  $$Q_{c, N}^{(\infty)} \equiv \mathcal{W}_c\big|_{\mathcal{H}_{c, N}}.$$
  Because test functions $g_v \in \mathcal{H}_{c, N}$ have Fourier support contained in $[-\frac{\log c}{2\pi}, \frac{\log c}{2\pi}]$, all prime powers $q > c$ vanish identically on $g_v$, yielding $\mathcal{W}_c[g_v] = \mathcal{W}[g_v] = \sum_\rho g_v(\gamma_\rho)$.

### 2. Issues Reconciled & Open Boundary Questions
- **Reconciled (Archimedean Regularization):** The exact local subtraction counter-term $K_v(1) - K_v(1 - y/L) \sim \mathcal{O}(y)$ cancels the $1/y$ singularity of $\frac{1}{2\sinh y}$ at $y = 0$, proving that the discrete series $\sum [\frac{\|v\|^2}{n+1} - J(q_n)]$ matches Weil's regularized distribution with the exact constant $h_+(0)$.
- **Reconciled (Pole Normalization):** The factor-of-four typo in preliminary drafts is resolved; $\langle v, Q_{\mathrm{pole}} v \rangle = 2 g_v(i/2)$ matches the Guinand–Weil residue $g(i/2) + g(-i/2)$ identically.
- **Open (Finite-$c$ Prime Tail):** For arbitrary test functions $g$ outside $\mathcal{H}_{c, N}$, the omitted prime sum $\mathcal{R}_{\mathrm{prime}}(c; g) = -\frac{1}{\pi}\sum_{q > c} \frac{\Lambda(q)}{\sqrt{q}}\widehat{g}(\frac{\log q}{2\pi})$ does not vanish; recovery of the full Weil functional $\mathcal{W}$ strictly requires $c \to \infty$.
- **Open (Gate 2/3 Boundary — Admissible Topology for Density):** Proving that $\bigcup_{c, N} \mathcal{H}_{c, N}$ is dense in Weil's admissible test class requires specifying the topological vector space (e.g. inductive limits of Paley–Wiener spaces with horizontal strip decay) and proving joint $(N, c) \to \infty$ approximation. $L^2$ density of trigonometric polynomials alone does not establish this.

---

## 1. The Finite-$N$ Prime Term: Pointwise Arithmetic Sampling

### 1.1 First-Principles Loewner Contraction
In the Connes–van Suijlekom formulation (Proposition 4.1), the Galerkin matrix $Q$ is assembled via divided differences of a source function $\psi(x)$:
$$(Q_\psi)_{mn} = \begin{cases} \dfrac{\psi(m) - \psi(n)}{m - n}, & m \ne n, \\[1.2ex] \psi'(m), & m = n, \end{cases} \qquad m, n \in \{-N, \dots, N\}.$$

The prime-piece source function $\psi_p^{(c)}(x)$ is an explicit finite sum over prime powers $q = p^k \le c$:
$$\psi_p^{(c)}(x) = -\frac{1}{\pi} \sum_{q = p^k \le c} \frac{\Lambda(q)}{\sqrt{q}} \sin(2\pi \omega_q x), \qquad \omega_q \equiv 1 - \frac{\log q}{L} \in (0, 1],$$
where $\Lambda(q) = \log p$ is the von Mangoldt function and $L = \log c$.

By linearity of the divided difference operator, the prime matrix decomposes as:
$$Q_{\mathrm{prime}}^{(c)} = \sum_{q = p^k \le c} Q_{q, N},$$
where each single prime-power matrix $Q_{q, N}$ is generated by the elementary single-frequency source:
$$\psi_q(x) = \frac{\alpha_q}{\pi} \sin(2\pi \omega_q x), \qquad \alpha_q \equiv -\frac{\Lambda(q)}{\sqrt{q}} = -\frac{\log p}{p^{k/2}}.$$

### 1.2 Coordinate-Space Contraction Lemma
Let $u = (u_{-N}, \dots, u_N)^\top \in \mathbb{R}^{2N+1}$ denote the Fourier coefficient vector with even symmetry $u_{-m} = u_m$, corresponding to the canonical even vector $v \in \mathbb{R}^{N+1}$ via $u_0 = v_0$, $u_m = v_m / \sqrt{2}$ ($m \ge 1$).

The vector $u$ generates the real-even trigonometric wave on $[0, 1]$:
$$\tau_v(s) = \sum_{m=-N}^N u_m e^{2\pi i m s} = v_0 + \sqrt{2}\sum_{m=1}^N v_m \cos(2\pi m s) \qquad (0 \le s \le 1).$$

The spatial Volterra auto-convolution on $[0, 1]$ is:
$$K_v(\omega) \equiv 2 \int_0^\omega \tau_v(s) \tau_v(\omega - s) \, ds.$$

**Theorem 1.1 (Exact Finite Source Calculus):**  
*For any frequency $\omega \in [0, 1]$ and amplitude $\alpha \in \mathbb{R}$, the divided-difference matrix $Q_{\alpha, \omega}$ attached to $\psi(x) = \frac{\alpha}{\pi}\sin(2\pi \omega x)$ satisfies:*
$$\langle u, Q_{\alpha, \omega} u \rangle = \alpha K_v(\omega) \qquad \forall N \ge 1.$$

*Proof:*  
Expanding the Volterra convolution in the exponential basis:
$$K_v(\omega) = 2 \sum_{m, n = -N}^N u_m u_n e^{2\pi i n \omega} \int_0^\omega e^{2\pi i (m - n) s} \, ds.$$
For off-diagonal entries $m \ne n$:
$$\int_0^\omega e^{2\pi i (m - n) s} \, ds = \frac{e^{2\pi i (m - n) \omega} - 1}{2\pi i (m - n)}.$$
Multiplying by $2 e^{2\pi i n \omega}$ yields:
$$2 e^{2\pi i n \omega} \frac{e^{2\pi i (m - n) \omega} - 1}{2\pi i (m - n)} = \frac{e^{2\pi i m \omega} - e^{2\pi i n \omega}}{\pi i (m - n)}.$$
Taking the real part (the imaginary part vanishes under the reflection symmetry $u_{-k} = u_k$):
$$\operatorname{Re}\left[ \frac{e^{2\pi i m \omega} - e^{2\pi i n \omega}}{\pi i (m - n)} \right] = \frac{\sin(2\pi \omega m) - \sin(2\pi \omega n)}{\pi (m - n)} = (Q_{1, \omega})_{mn}.$$
For diagonal entries $m = n$:
$$2 e^{2\pi i m \omega} \int_0^\omega 1 \, ds = 2\omega e^{2\pi i m \omega} \quad \implies \quad \operatorname{Re}[2\omega e^{2\pi i m \omega}] = 2\omega \cos(2\pi \omega m) = (Q_{1, \omega})_{mm}.$$
Multiplying by $\alpha$ proves the identity unconditionally for all $N$. $\blacksquare$

### 1.3 Precise Meaning of "Zero Arithmetic Aliasing"
Applying Theorem 1.1 to the prime source with $\omega_q = 1 - \frac{\log q}{L}$:
$$\boxed{\langle v, Q_{\mathrm{prime}}^{(c)} v \rangle = -\sum_{p^k \le c} \frac{\log p}{p^{k/2}} K_v\left(1 - \frac{\log p^k}{L}\right).}$$

In physical coordinates $t = L s \in [0, L]$, the wavepacket is $f_v(t) = \tau_v(t / L)$. The physical Volterra convolution is:
$$K_v^{\mathrm{phys}}(x) = 2 \int_0^x f_v(t) f_v(x - t) \, dt = L K_v(x / L).$$
Evaluating at $x = L(1 - \omega_q) = L - \log p^k$:
$$K_v\left(1 - \frac{\log p^k}{L}\right) = \frac{1}{L} K_v^{\mathrm{phys}}(L - \log p^k).$$

Because $f_v$ is even and $L$-periodic, $f_v(L - t) = f_v(t)$. Substituting $s = L - t$ transforms the reflected convolution into the **spatial autocorrelation function**:
$$K_v^{\mathrm{phys}}(L - y) = 2 \int_0^{L - y} f_v(t) f_v(L - y - t) \, dt = 2 \int_y^L f_v(s) f_v(s - y) \, ds \equiv 2 (f_v \ast f_v^\dagger)(y).$$

Therefore, the contribution of each prime power $q = p^k \le c$ to the finite Galerkin quadratic form is:
$$\boxed{Q_{q, N}[v] = -\frac{\log p}{p^{k/2}} \cdot \frac{2}{L} \int_{\log p^k}^L f_v(t) f_v(t - \log p^k) \, dt.}$$

**Epistemic Distinction (Arithmetic Aliasing vs Galerkin Restriction):**
- **Zero Arithmetic Aliasing (Proved):** For every $v \in \mathcal{H}_{c, N}$, the arithmetic shift $\log p^k$ is represented **pointwise and exactly** as the argument of the wavepacket's spatial autocorrelation. The discrete Fourier lattice $\Lambda_{\mathrm{Fourier}} = \{2\pi m / L\}$ does not smear or alias the prime positions $\log p^k$.
- **Finite-Dimensional Restriction (Active):** This exactness does *not* imply that $\mathcal{H}_{c, N}$ spans the full Weil test class. The wavepacket $f_v$ is restricted to band $N$, and the spatial interval is restricted to length $L = \log c$. This finite-dimensional restriction is an approximation of the test space, not an aliasing of the arithmetic points.

---

## 2. The Finite-$N$ Archimedean Term: Mellin Multipliers and Hyperbolic Kernels

### 2.1 The Continuous Archimedean Quadratic Form
The Archimedean piece of the CvS matrix is governed by the smooth multiplier:
$$h_+(r) = \operatorname{Re}\psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi.$$
In the Galerkin basis with Archimedean truncation $T > 0$, the matrix $Q_{\mathrm{arch}, T}$ is the divided-difference matrix of:
$$\psi_{\mathbb{R}, T}(x) = \frac{1}{2\pi^2} \int_{-T}^T h_+(r) \mathcal{S}(r, x, L) \, dr, \qquad \mathcal{S}(r, x, L) \equiv \int_0^L \sin\left(2\pi x \left(1 - \frac{y}{L}\right)\right) \cos(r y) \, dy.$$

Applying Theorem 1.1 to the $y$-integral inside $\psi_{\mathbb{R}, T}$:
$$\langle v, Q_{\mathrm{arch}, T} v \rangle = \frac{1}{\pi} \int_0^T h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr,$$
where the Fourier-side power spectral density is:
$$K_{\mathrm{Fourier}}(v, r, L) = \int_0^L K_v\left(1 - \frac{y}{L}\right) \cos(r y) \, dy = \frac{1}{L} \left| \int_0^L f_v(t) e^{-i r t} \, dt \right|^2 \ge 0.$$

By definition of the Guinand–Weil test function $g_v(r) = \int_{-\Delta}^\Delta \widehat{g}_v(\xi) \cos(2\pi r \xi) d\xi$ with $\widehat{g}_v(\xi) = \pi K_v(1 - |\xi|/\Delta)$, changing variables $\xi = y / (2\pi)$ shows:
$$g_v(r) \equiv K_{\mathrm{Fourier}}(v, r, L).$$
Consequently, the frequency-space identity is exact:
$$\boxed{\lim_{T \to \infty} \langle v, Q_{\mathrm{arch}, T} v \rangle = \frac{1}{\pi}\int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr = \frac{1}{2\pi}\int_{-\infty}^\infty h_+(r) g_v(r) \, dr \equiv W_\infty[g_v].}$$

### 2.2 The Archimedean Normalization Lemma (Exact Coordinate-Space Derivation)

We now prove that the continuous Archimedean quadratic form $\mathcal{Q}_{\mathrm{arch}}(v)$ reconstructs André Weil's regularized Archimedean distribution in coordinate space **all the way through, with exact local constants and counter-terms**.

**Theorem 2.1 (Archimedean Normalization Lemma):**  
*Let $F(y) \equiv K_v(1 - y/L) \mathbf{1}_{[0, L]}(y)$ denote the spatial autocorrelation of the wavepacket $f_v(t) \mathbf{1}_{[0, L]}(t)$, with $F(0) = K_v(1) = 2 \|v\|_2^2$.*  
*Then the cutoff-free Archimedean quadratic form satisfies the exact algebraic identity:*
$$\boxed{\langle v, Q_{\mathrm{arch}, \infty} v \rangle \equiv \frac{1}{2\pi} \int_{-\infty}^\infty h_+(r) g_v(r) \, dr \equiv \frac{1}{2} h_+(0) F(0) + \int_0^\infty \frac{F(0) - F(y)}{2\sinh y} \, dy,}$$
*where:*
$$h_+(0) = \psi(1/4) - \log \pi = -\left(\log \pi + \gamma + \frac{\pi}{2} + 3\log 2\right) \approx -5.37218488.$$
*Decomposed onto the finite interval $[0, L]$, this evaluates identically as:*
$$\boxed{\langle v, Q_{\mathrm{arch}, \infty} v \rangle = h_+(0) \|v\|_2^2 + \int_0^L \frac{K_v(1) - K_v(1 - y/L)}{2\sinh y} \, dy + K_v(1) \operatorname{artanh}(e^{-L}),}$$
*where the local numerator $K_v(1) - K_v(1 - y/L) = \frac{y}{L} K_v'(1) + \mathcal{O}(y^2)$ vanishes at $y = 0$, rendering the integral unconditionally regular and finite.*

*Proof:*  
We prove this in four self-contained steps.

**Step 1: Multiplier Integral Identity.**  
From the Weierstrass partial fraction expansion:
$$h_+(r) = C_{\mathrm{arch}} + \sum_{n=0}^\infty \left[ \frac{1}{n+1} - \frac{2q_n}{q_n^2 + r^2} \right], \qquad C_{\mathrm{arch}} = -\gamma - \log \pi, \quad q_n = 2n + \frac{1}{2}.$$
At $r = 0$:
$$h_+(0) = C_{\mathrm{arch}} + \sum_{n=0}^\infty \left[ \frac{1}{n+1} - \frac{2}{q_n} \right] = -\gamma - \log \pi + \sum_{n=0}^\infty \left( \frac{1}{n+1} - \frac{1}{n + 1/4} \right) = \psi(1/4) - \log \pi.$$
Subtracting $h_+(0)$ from $h_+(r)$, the constant $C_{\mathrm{arch}}$ and the terms $\frac{1}{n+1}$ cancel identically:
$$h_+(r) - h_+(0) = \sum_{n=0}^\infty \left[ \frac{2}{q_n} - \frac{2q_n}{q_n^2 + r^2} \right] = \sum_{n=0}^\infty \frac{2r^2}{q_n(q_n^2 + r^2)}.$$
Using the elementary Laplace transform $\int_0^\infty (1 - \cos(ry)) e^{-q_n y} dy = \frac{1}{q_n} - \frac{q_n}{q_n^2 + r^2} = \frac{1}{2}[\frac{2}{q_n} - \frac{2q_n}{q_n^2 + r^2}]$:
$$h_+(r) - h_+(0) = 2 \sum_{n=0}^\infty \int_0^\infty (1 - \cos(ry)) e^{-(2n + 1/2)y} \, dy.$$
Summing the geometric series $\sum_{n=0}^\infty e^{-(2n + 1/2)y} = \frac{e^{-y/2}}{1 - e^{-2y}} = \frac{1}{2\sinh y}$ under the integral:
$$h_+(r) - h_+(0) = 2 \int_0^\infty \frac{1 - \cos(ry)}{2\sinh y} \, dy = \int_0^\infty \frac{1 - \cos(ry)}{\sinh y} \, dy.$$
This establishes the exact multiplier representation:
$$\boxed{h_+(r) \equiv h_+(0) + \int_0^\infty \frac{1 - \cos(ry)}{\sinh y} \, dy.}$$

**Step 2: Transfer to Coordinate Space.**  
Multiplying by $g_v(r)$ and integrating $\frac{1}{2\pi} \int_{-\infty}^\infty dr$:
$$\frac{1}{2\pi}\int_{-\infty}^\infty h_+(r) g_v(r) \, dr = h_+(0) \left( \frac{1}{2\pi}\int_{-\infty}^\infty g_v(r) \, dr \right) + \int_0^\infty \frac{1}{\sinh y} \left[ \frac{1}{2\pi} \int_{-\infty}^\infty g_v(r) (1 - \cos(ry)) \, dr \right] dy.$$
Because $g_v(r) = \int_0^L F(y) \cos(ry) \, dy$ is the one-sided cosine transform of $F(y) = K_v(1 - y/L) \mathbf{1}_{[0, L]}(y)$, the inverse cosine transform yields:
$$\frac{1}{2\pi} \int_{-\infty}^\infty g_v(r) \, dr = \frac{1}{2} F(0), \qquad \frac{1}{2\pi} \int_{-\infty}^\infty g_v(r) \cos(ry) \, dr = \frac{1}{2} F(y).$$
Since $F(0) = K_v(1) = 2\int_0^1 \tau_v(s)^2 \, ds = 2 \|v\|_2^2$, the first term evaluates to:
$$h_+(0) \left( \frac{1}{2} F(0) \right) = h_+(0) \|v\|_2^2.$$
For the second term, subtracting the cosine inversion gives:
$$\frac{1}{2\pi} \int_{-\infty}^\infty g_v(r) (1 - \cos(ry)) \, dr = \frac{1}{2} [F(0) - F(y)].$$
Interchanging the $r$- and $y$-integrals via Fubini's theorem (as $g_v \in \mathcal{S}(\mathbb{R})$ has rapid decay) yields:
$$\int_0^\infty \frac{1}{\sinh y} \left[ \frac{1}{2\pi}\int_{-\infty}^\infty g_v(r) (1 - \cos(ry)) \, dr \right] dy = \int_0^\infty \frac{F(0) - F(y)}{2\sinh y} \, dy.$$
Adding the two terms yields:
$$\frac{1}{2\pi}\int_{-\infty}^\infty h_+(r) g_v(r) \, dr = \frac{1}{2} h_+(0) F(0) + \int_0^\infty \frac{F(0) - F(y)}{2\sinh y} \, dy = h_+(0) \|v\|_2^2 + \int_0^\infty \frac{F(0) - F(y)}{2\sinh y} \, dy.$$

**Step 3: Interval Splitting and Tail Evaluation.**  
Since $f_v(t)$ is supported on $[0, L]$, its autocorrelation $F(y) = K_v(1 - y/L)$ is identically zero for $y \ge L$. Splitting the integral into $[0, L]$ and $[L, \infty)$:
$$\int_0^\infty \frac{F(0) - F(y)}{2\sinh y} \, dy = \int_0^L \frac{K_v(1) - K_v(1 - y/L)}{2\sinh y} \, dy + \int_L^\infty \frac{K_v(1)}{2\sinh y} \, dy.$$
The tail integral evaluates in closed form:
$$\int_L^\infty \frac{1}{2\sinh y} \, dy = \int_L^\infty \frac{e^{-y}}{1 - e^{-2y}} \, dy = \frac{1}{2}\left[ \log\left( \frac{1 + e^{-y}}{1 - e^{-y}} \right) \right]_L^\infty = \frac{1}{2} \log\left( \frac{1 + e^{-L}}{1 - e^{-L}} \right) = \operatorname{artanh}(e^{-L}).$$

**Step 4: Equivalence with Discrete Cauchy Sum.**  
From Paper 4 (Theorem 4.3), the series expansion is $\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty [\frac{\|v\|_2^2}{n+1} - J(q_n)]$.  
Splitting the summand:
$$\frac{\|v\|_2^2}{n+1} - J(q_n) = \left( \frac{\|v\|_2^2}{n+1} - \frac{2\|v\|_2^2}{q_n} \right) + \left( \frac{2\|v\|_2^2}{q_n} - J(q_n) \right).$$
Summing the first bracket: $\sum_{n=0}^\infty (\frac{\|v\|_2^2}{n+1} - \frac{2\|v\|_2^2}{q_n}) = \|v\|_2^2 (\psi(1/4) + \gamma)$.  
Adding $C_{\mathrm{arch}} \|v\|_2^2 = (-\gamma - \log \pi) \|v\|_2^2$ reproduces $h_+(0) \|v\|_2^2$ identically.  
For the second bracket: substituting $J(q_n) = \int_0^L K_v(1 - y/L) e^{-q_n y} dy$ and $\frac{2\|v\|_2^2}{q_n} = K_v(1) \int_0^\infty e^{-q_n y} dy$ yields the regularized coordinate integral identically. $\blacksquare$

---

## 3. The Finite-$N$ Pole Term: Normalization Audit

### 3.1 Trace of Normalization Factors
We verify the pole term normalization step by step from first principles:

1. **The CvS Pole Source (Connes & van Suijlekom Prop 4.1):**
   $$\psi_{\mathrm{pole}}(x) = \frac{1}{\pi} \int_0^L 2 \cosh(y/2) \sin\left(2\pi x \left(1 - \frac{y}{L}\right)\right) \, dy.$$
2. **Loewner Contraction via Theorem 1.1:**
   Here, the measure on $[0, 1]$ is $d\mu(\omega) = 2 \cosh(L(1-\omega)/2) L d\omega = 2 \cosh(y/2) dy$.
   By Theorem 1.1, the contraction is:
   $$\langle v, Q_{\mathrm{pole}} v \rangle = \int_0^1 K_v(\omega) d\mu(\omega) = \int_0^L 2 \cosh(y/2) K_v\left(1 - \frac{y}{L}\right) \, dy = 2 \int_0^L \cosh(y/2) K_v\left(1 - \frac{y}{L}\right) \, dy.$$
3. **Evaluation of the Guinand–Weil Test Function at $z = i/2$:**
   The test function is defined on $\mathbb{C}$ by:
   $$g_v(z) = \int_{-\Delta}^\Delta \widehat{g}_v(\xi) e^{2\pi i z \xi} \, d\xi, \qquad \Delta = \frac{L}{2\pi}.$$
   Since $\widehat{g}_v(\xi) = \pi K_v(1 - |\xi|/\Delta)$ is even in $\xi$:
   $$g_v(i/2) = \int_{-\Delta}^\Delta \widehat{g}_v(\xi) e^{-\pi \xi} \, d\xi = \int_{-\Delta}^\Delta \widehat{g}_v(\xi) \cosh(\pi \xi) \, d\xi = 2 \int_0^\Delta \widehat{g}_v(\xi) \cosh(\pi \xi) \, d\xi.$$
   Substituting $\xi = \frac{y}{2\pi}$ ($d\xi = \frac{dy}{2\pi}$, $\pi \xi = \frac{y}{2}$, and $\widehat{g}_v(\xi) = \pi K_v(1 - y/L)$):
   $$g_v(i/2) = 2 \int_0^L \pi K_v\left(1 - \frac{y}{L}\right) \cosh(y/2) \frac{dy}{2\pi} = \int_0^L K_v\left(1 - \frac{y}{L}\right) \cosh(y/2) \, dy.$$

### 3.2 The Exact Identity
Comparing Step 2 and Step 3:
$$\langle v, Q_{\mathrm{pole}} v \rangle = 2 \left[ \int_0^L K_v\left(1 - \frac{y}{L}\right) \cosh(y/2) \, dy \right] = 2 g_v(i/2).$$

In the explicit formula of Guinand and Weil (Weil 1952, Bombieri 2000 Section 2), the residue of $\zeta(s)$ at $s = 0, 1$ produces the pole term:
$$\mathcal{W}_{\mathrm{pole}}[g] = g(i/2) + g(-i/2).$$
For any even test function, $g(-i/2) = g(i/2)$, so:
$$\mathcal{W}_{\mathrm{pole}}[g_v] = 2 g_v(i/2).$$

Therefore:
$$\boxed{\langle v, Q_{\mathrm{pole}}^{(N)} v \rangle \equiv 2 g_v(i/2) \equiv g_v(i/2) + g_v(-i/2) \qquad \forall N \ge 1, \; \forall c > 1.}$$
The normalization is exact; the factor of 4 in early working notes was an algebraic typo arising from double-counting the even-reflection factor of 2.

---

## 4. Finite-$c$ Identity, Prime Tail, and the Density Question

### 4.1 The Subspace Identification Theorem
Let $\mathcal{H}_{c, N} \subset L^2([0, L])$ denote the $(N+1)$-dimensional Galerkin subspace of even trigonometric polynomials of band $N$, and let $g_v$ denote the induced test function for $v \in \mathcal{H}_{c, N}$.

Define the **prime-truncated Weil functional** $\mathcal{W}_c$ at cutoff $c$:
$$\mathcal{W}_c[g] \equiv g(i/2) + g(-i/2) - \frac{1}{\pi} \sum_{q = p^k \le c} \frac{\Lambda(q)}{\sqrt{q}} \widehat{g}\left(\frac{\log q}{2\pi}\right) + \frac{1}{2\pi} \int_{-\infty}^\infty h_+(r) g(r) \, dr.$$

**Theorem 4.1 (Subspace Weil Identity):**  
*On the finite Galerkin subspace $\mathcal{H}_{c, N}$, the cutoff-free matrix $Q_{c, N}^{(\infty)} = Q_{\mathrm{prime}}^{(c)} + Q_{\mathrm{pole}} + Q_{\mathrm{arch}, \infty}$ satisfies:*
$$\langle v, Q_{c, N}^{(\infty)} v \rangle \equiv \mathcal{W}_c[g_v] \qquad \forall v \in \mathcal{H}_{c, N}.$$
*Furthermore, because $\operatorname{supp}\widehat{g}_v \subseteq [-\frac{\log c}{2\pi}, \frac{\log c}{2\pi}]$, every prime power $q = p^k > c$ satisfies $\frac{\log q}{2\pi} > \Delta$, so $\widehat{g}_v(\frac{\log q}{2\pi}) = 0$. Consequently, the omitted prime sum vanishes identically on the subspace:*
$$\boxed{\langle v, Q_{c, N}^{(\infty)} v \rangle \equiv \mathcal{W}[g_v] \equiv \sum_{\rho \in Z_\zeta} g_v\left(\frac{\rho - 1/2}{i}\right) \qquad \forall v \in \mathcal{H}_{c, N}.}$$

### 4.2 Discrepancy Hierarchy for Arbitrary Test Functions
For an arbitrary test function $g \in \mathcal{S}_{\mathrm{Weil}}$ in Weil's admissible class (not restricted to $\mathcal{H}_{c, N}$), the error between the finite matrix evaluation and $\mathcal{W}[g]$ decomposes into three distinct components:

$$\mathcal{W}[g] - \langle v, Q_{c, N}^{(T)} v \rangle = \mathcal{R}_{\mathrm{cutoff}}(T; v) + \mathcal{R}_{\mathrm{prime}}(c; g) + \mathcal{R}_{\mathrm{Galerkin}}(N, c; g),$$

where:
1. **The Archimedean Cutoff Tail $\mathcal{R}_{\mathrm{cutoff}}(T; v)$:**
   $$\mathcal{R}_{\mathrm{cutoff}}(T; v) = \frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr.$$
   Vanishes as $T \to \infty$ with geometric rate $(a_N/T)^2$ for $T > a_N = 2\pi N / L$ (Paper 4 Theorem 5.1).
2. **The Prime Truncation Tail $\mathcal{R}_{\mathrm{prime}}(c; g)$:**
   $$\mathcal{R}_{\mathrm{prime}}(c; g) = -\frac{1}{\pi} \sum_{p^k > c} \frac{\log p}{p^{k/2}} \widehat{g}\left(\frac{\log p^k}{2\pi}\right).$$
   For an arbitrary test function $g$, $\mathcal{R}_{\mathrm{prime}}(c; g) \ne 0$ at finite $c$. Recovery of the full prime series strictly requires the prime cutoff limit $c \to \infty$.
3. **The Galerkin Approximation Defect $\mathcal{R}_{\mathrm{Galerkin}}(N, c; g)$:**
   $$\mathcal{R}_{\mathrm{Galerkin}}(N, c; g) = \mathcal{W}_c[g - g_v].$$
   The distance between $g$ and the finite-rank subspace in the Weil functional topology.

### 4.3 The Gate 2/3 Boundary: Admissible Topology for Density
A recurring risk in the literature is asserting that "trigonometric polynomials are dense in $L^2([0, L])$, therefore Galerkin test functions are dense in Weil's class."
This statement conflates distinct topological vector spaces:
- The space of Guinand–Weil test functions $\mathcal{S}_{\mathrm{Weil}}$ consists of entire functions of exponential type with uniform decay on horizontal strips:
  $$g(x + i y) = \mathcal{O}\big((1 + |x|)^{-(1 + \delta)}\big), \qquad \delta > 0.$$
- The explicit formula $\mathcal{W}[g] = \sum_\rho g(\gamma_\rho)$ is continuous only with respect to a topology ensuring the absolute convergence of the zero sum (e.g. inductive limit topology on Fréchet spaces of Paley–Wiener type).
- **Open Gate 2/3 Question:** Establishing that $\bigcup_{c > 1} \bigcup_{N \ge 1} \{g_v : v \in \mathcal{H}_{c, N}\}$ is dense in $\mathcal{S}_{\mathrm{Weil}}$ requires an explicit approximation theorem in the horizontal-strip topology as $(N, c) \to (\infty, \infty)$.
This is an active mathematical obligation at the Gate 2 / Gate 3 boundary and is deliberately kept open.

---

## 5. Summary Audit Matrix

| Audit Item | Theoretical Question | Status | Finding / Resolution |
| :--- | :--- | :---: | :--- |
| **1. Prime Sampling** | Does $Q_{\mathrm{prime}}^{(c)}$ alias prime locations $\log p^k$? | **RESOLVED (EXACT)** | Exact pointwise autocorrelation evaluation; zero aliasing. |
| **2. Pole Normalization** | Factor of 2 vs 4 in $\langle v, Q_{\mathrm{pole}} v \rangle$? | **RESOLVED (EXACT)** | Exactly $2 g_v(i/2) = g_v(i/2) + g_v(-i/2)$, matching Guinand–Weil residue. |
| **3. Archimedean Kernel** | Exact local counter-term and constant in $Q_{\mathrm{arch}}$? | **RESOLVED (EXACT)** | Proved: $h_+(0)\|v\|^2 + \int_0^L \frac{K_v(1) - K_v(1-y/L)}{2\sinh y} dy + K_v(1)\operatorname{artanh}(e^{-L})$. |
| **4. Subspace Identity** | Does $Q_{c, N}^{(\infty)} = \mathcal{W}_c$ on $\mathcal{H}_{c, N}$? | **RESOLVED (EXACT)** | Exact identity $\langle v, Q_{c, N}^{(\infty)} v \rangle \equiv \mathcal{W}_c[g_v] = \sum_\rho g_v(\gamma_\rho)$. |
| **5. Finite-$c$ Prime Tail** | Does $Q_{c, N}^{(\infty)}$ equal full $\mathcal{W}$ for all $g$? | **ISOLATED** | No; full $\mathcal{W}$ requires $c \to \infty$ to include primes $p^k > c$. |
| **6. Galerkin Density** | Does $L^2$ density imply Weil density? | **OPEN (GATE 2/3)** | Requires horizontal-strip topology and joint $(N, c) \to \infty$ analysis. |

$$\boxed{\textbf{GATE 3 FINITE-SUBSPACE IDENTITY: PROVED EXACTLY ON } \mathcal{H}_{c, N}}$$
$$\boxed{\textbf{GATE 2/3 DENSITY AND CONTINUUM EXTENSION: OPEN}}$$

**Strategic Takeaway:**  
No structural obstruction has been found; every term of $Q_{c, N}^{(\infty)}$ on the Galerkin subspace $\mathcal{H}_{c, N}$ matches André Weil's explicit quadratic functional identically, with exact normalizations and regularizations.
We can now pivot to **Gate 1 Work Item B (Analytic Reduction of the Joint-Limit Tail Extinction Target)** with complete confidence in the arithmetic substrate.
