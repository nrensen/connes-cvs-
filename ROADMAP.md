# Project Roadmap: Finite-Rank Connes–van Suijlekom Galerkin Positivity and the Weil Explicit Formula

**Document Version:** 1.0  
**Date:** September 2026  
**Status:** Canonical Project Strategy & Master Architecture  
**Associated Manuscripts:** 
- Paper 4: *An Exact Resolvent and Commutator Toolkit* (Locked Toolkit Baseline)
- Paper 4B: *The Dirichlet Continuum Limit, Barrier Mechanics, and Asymptotic Weil Positivity* (Asymptotic Programme)

---

## 1. Executive Summary & The Core Paradigm Shift

The Connes–van Suijlekom (2025) and Connes–Consani–Moscovici (2026) framework projects André Weil's explicit quadratic functional of prime number theory onto finite-rank Galerkin subspaces of frequency band $N$ on a logarithmic scaling interval $[0, L] = [0, \log c]$:

$$\mathcal{Q}_{\mathrm{Weil}}(v) = \mathcal{Q}_{\mathrm{pole}}(v) + \mathcal{Q}_{\mathrm{prime}}(v) + \mathcal{Q}_{\mathrm{arch}}(v).$$

In Paper 4, we proved the unconditional global non-negativity of the Fourier-side kernel:
$$K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0 \qquad \forall r \in \mathbb{R}, \; \forall v \in \mathbb{R}^{N+1}.$$

However, **kernel non-negativity alone is not sufficient to prove positivity of the Archimedean form or the Weil form**, because the Archimedean multiplier:
$$h_+(r) = \operatorname{Re}\psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi$$
is **strictly negative** on a bounded low-frequency interval:
$$h_+(0) = -\gamma - 2\log 2 - \log \pi \approx -5.37218 < 0,$$
possessing a unique positive root at:
$$r_* \approx 6.28984.$$

For $r > r_*$, Stirling's formula guarantees that $h_+(r) = \log(r/(2\pi)) + \mathcal{O}(r^{-2}) > 0$ is strictly positive and monotonically increasing. Consequently, the Archimedean quadratic form:
$$\mathcal{Q}_{\mathrm{arch}}(v) = \frac{1}{\pi} \int_0^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$$
is intrinsically a **signed weighted spectral energy**, not an independently positive quadratic form.

### The Fundamental Paradigm Shift
We must not attempt to "defeat $h_+$" by demanding more positivity from the Fourier kernel $K_{\mathrm{Fourier}}$. We have already extracted the maximum algebraic positivity that the Fourier kernel contains. Instead, the central question of the research programme is:

$$\boxed{\textbf{Why does the complete Weil functional remain positive despite the negative spectral weight of } h_+?}$$

The entire mathematical challenge reduces to determining whether the negative spectral weight on $[0, r_*]$ is controlled by the rest of the Weil functional.

---

## 2. Problem Formulation: The Operator Dominance Principle

Decompose the Archimedean multiplier into positive and negative parts:
$$h_+(r) = h_+^+(r) - h_+^-(r),$$
where:
$$h_+^-(r) \equiv |h_+(r)| \cdot \mathbf{1}_{[0, r_*]}(r) \ge 0, \qquad h_+^+(r) \equiv h_+(r) \cdot \mathbf{1}_{[r_*, \infty)}(r) \ge 0.$$

Crucially, **$h_+^-(r)$ is supported strictly on the compact, bounded low-frequency interval $[0, r_*]$**. 

This induces a corresponding operator splitting of the Archimedean matrix:
$$\mathcal{Q}_{\mathrm{arch}} = \mathcal{Q}_{\mathrm{arch}}^{(+)} - \mathcal{Q}_{\mathrm{arch}}^{(-)},$$
where both $\mathcal{Q}_{\mathrm{arch}}^{(+)}$ and $\mathcal{Q}_{\mathrm{arch}}^{(-)}$ are **manifestly positive semi-definite quadratic forms**:
$$\mathcal{Q}_{\mathrm{arch}}^{(-)}(v) = \frac{1}{\pi} \int_0^{r_*} |h_+(r)| K_{\mathrm{Fourier}}(v, r, L) \, dr \succeq 0.$$

The total Weil quadratic form then reads:
$$\mathcal{Q}_{\mathrm{Weil}} = \underbrace{\mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}}^{(+)}}_{\mathcal{Q}_{\mathrm{positive\ side}}} - \mathcal{Q}_{\mathrm{arch}}^{(-)}.$$

Weil positivity on the Galerkin subspace is therefore mathematically equivalent to the **operator dominance inequality**:

$$\boxed{\mathcal{Q}_{\mathrm{positive\ side}} \succeq \mathcal{Q}_{\mathrm{arch}}^{(-)} \qquad \Longleftrightarrow \qquad \lambda_{\min}\big(\mathcal{Q}_{\mathrm{positive\ side}}, \mathcal{Q}_{\mathrm{arch}}^{(-)}\big) \ge 1.}$$

### Decoupling from the Ground State
This formulation is entirely independent of the ground-state eigenvector $v_0$. It poses Weil positivity not as an asymptotic barrier-tunneling phenomenon, but as a global operator inequality between the positive terms of the explicit formula and a finite-band, compact negative perturbation.

---

## 3. Resolvent Structure of the Negative Operator

Because the negative interval $[0, r_*]$ is compact and bounded, we have:
$$\mathcal{Q}_{\mathrm{arch}}^{(-)}(v) = \frac{1}{\pi} \int_0^{r_*} |h_+(r)| \Phi_v(r)^2 \, dr.$$
From Paper 4 (Theorem 2.1 & Proposition 2.2), the Fourier amplitude $\Phi_v(r)$ admits the exact boundary resolvent representation:
$$\Phi_v(r) = \frac{2}{\sqrt{L}} \frac{\sin(rL/2)}{r} D\left(-\frac{1}{r^2}\right),$$
where $D(z) = v_0 + \sqrt{2}\sum_{m=1}^N \frac{v_m}{1 + z a_m^2}$ is the canonical scalar boundary resolvent.
Substituting this into $\mathcal{Q}_{\mathrm{arch}}^{(-)}$ yields the exact formula:

$$\boxed{\mathcal{Q}_{\mathrm{arch}}^{(-)}(v) = \frac{4}{\pi L} \int_0^{r_*} |h_+(r)| \frac{\sin^2(rL/2)}{r^2} D\left(-\frac{1}{r^2}\right)^2 \, dr.}$$

### Algebraic Simplification & Resolvent Pole-Zero Geometry
This identity provides a structured algebraic framework:
1. **Domain of the Inverted Spectral Variable:** Under the change of variables $z = -1/r^2$, the compact integration domain $r \in [0, r_*]$ is mapped onto the semi-infinite negative ray:
   $$z \in (-\infty, z_*], \qquad z_* = -\frac{1}{r_*^2} \approx -\frac{1}{(6.28984)^2} \approx -0.025277.$$
2. **Poles of the Boundary Resolvent:** The rational function $D(z)$ has isolated simple poles located on the negative real axis at the inverse squared Fourier frequencies:
   $$z_m = -\frac{1}{a_m^2} = -\frac{L^2}{4\pi^2 m^2}, \qquad m = 1, 2, \dots, N.$$
3. **Intersection with the Negative Domain:** For our primary benchmark $c = 13$ ($L = \log 13 \approx 2.56495$):
   $$z_1 = -\frac{1}{a_1^2} \approx -0.16664 \in (-\infty, z_*],$$
   $$z_2 = -\frac{1}{a_2^2} \approx -0.04166 \in (-\infty, z_*],$$
   $$z_3 = -\frac{1}{a_3^2} \approx -0.01852 > z_* \quad (\text{outside the negative zone!}).$$
   Thus, the transformed negative-frequency interval intersects only the first two resolvent pole locations $(z_1, z_2)$. All higher poles $z_m$ ($m \ge 3$) lie outside this interval.
4. **Removable Pole-Zero Geometry vs. Singularities:** At the Fourier lattice frequencies $r = a_m$, $D(-1/r^2)$ has a pole, but simultaneously the envelope factor $\sin(rL/2) = \sin(\pi m) = 0$ vanishes. Because the full amplitude $\Phi_v(r)$ is an entire function, these apparent singularities are **completely removable**. The negative Archimedean integral is therefore **not** an integral governed by genuine singularities. Rather, for $c = 13$, the problem is a finite-dimensional **resolvent interpolation problem around two removable pole-zero pairs**. Understanding the algebraic structure of these removable pairs provides a concrete pathway toward bounding the negative Archimedean operator.

---

## 4. Analysis of Methodological Routes

We identify three primary routes of increasing mathematical depth to establish operator dominance, supplemented by a real-space transfer method:

```
[Route A: Archimedean Positivity via Uncertainty] (Deprecated / Unlikely)
  Attempt to show high-frequency mass from Paley-Wiener forces Q_arch >= 0 independently.
        |  (Rejected: Explicit formula treats arch/prime/pole as coupled distribution)
        v
[Route B: Exact Weierstrass Decomposition & Resolvent Pairing] (Primary Analytical Route)
  Decompose h+(r) into Weierstrass poles -> Q_arch = C||v||^2 + sum [||v||^2/(n+1) - J(q_n)].
  Pair negative J(q_n) resolvents against explicit pole/prime boundary terms.
        |
        +---> [Advanced Extension: DLMF Completely Monotone Physical-Space Transfer]
        |       Decompose h+(r) via digamma integral into completely monotone density;
        |       pull negative well into physical coordinate space against prime Dirac deltas.
        v
[Route C: Integrated Bivariate Kernel & Dual Lattice Sampling] (Grand Unified Architecture)
  Represent Q_Weil as a single spatial integral against K_Weil(x, y).
  Analyze the interaction between the Fourier lattice {a_m} and Arithmetic lattice {log p^k}.
  Formulate as de Branges / Beurling space positivity problem.
```

### Route A: Archimedean Positivity via Uncertainty (Deprecated)
*Concept:* Ask whether the finite bandwidth of the trigonometric space forces enough high-frequency mass $\int_{r_*}^\infty h_+ K_v dr$ to unconditionally overwhelm $\int_0^{r_*} |h_+| K_v dr$.  
*Assessment:* **Unlikely to succeed.** In Weil's explicit formula, the gamma-factor, prime-power, and pole terms are complementary distributions of an integrated adelic trace. Demanding that $\mathcal{Q}_{\mathrm{arch}}$ be positive in isolation ignores the physical role of the pole and prime forms in stabilizing the functional.

### Route B: Exact Weierstrass Decomposition and Pole Pairing (Primary Analytical Route)
*Concept:* From Paper 4 (Theorem 5.1 & Corollary 5.3), the Weierstrass partial fraction expansion of the digamma function gives:
$$h_+(r) = -\gamma - \log \pi + \sum_{n=0}^\infty \left[ \frac{1}{n+1} - \frac{2q_n}{q_n^2 + r^2} \right], \qquad q_n = 2n + \frac{1}{2}.$$
Integrating against $K_{\mathrm{Fourier}}$ yields the exact series:
$$\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty \frac{\|v\|_2^2}{n+1} - \sum_{n=0}^\infty J(q_n),$$
where every $J(q_n) = \frac{1}{\pi}\int_0^\infty \frac{2q_n}{q_n^2 + r^2} K_{\mathrm{Fourier}}(v, r, L) dr > 0$ is a **positive algebraic resolvent form** with exact closed-form evaluation:
$$J(q) = \frac{2 v_0^2}{q} + \sum_{m=1}^N \frac{2 q v_m^2}{q^2 + a_m^2} + \text{boundary leakage term } B_q(v).$$
The negative terms in the Archimedean functional are precisely the discrete sequence $-J(q_n)$.

#### The Remarkable Cancellation Hypothesis (Pairing $-J(q_n)$ with $\mathcal{Q}_{\mathrm{pole}}$)
The Archimedean quadratic form has the schematic structure:
$$\mathcal{Q}_{\mathrm{arch}}(v) = C_{\mathrm{arch}} \|v\|_2^2 + \sum_{n=0}^\infty \frac{\|v\|_2^2}{n+1} - \sum_{n=0}^\infty J(q_n).$$
The negative part is solely generated by the sequence $-\sum_{n=0}^\infty J(q_n)$.  
Now, recall that the pole contribution $\mathcal{Q}_{\mathrm{pole}}$ in Weil's explicit formula evaluates the test function at the poles of the Riemann zeta function ($s = 1$ and $s = 0$, or $s = \pm 1/2$ in the critical strip centering). In the boundary representation of Paper 4, the pole term introduces explicit boundary and resolvent evaluations at the endpoints $t = 0$ and $t = L$.  
The central hypothesis of Route B is that **the negative resolvent sequence $-\sum_{n=0}^\infty J(q_n)$ contains precisely the algebraic partners needed to pair with and absorb $\mathcal{Q}_{\mathrm{pole}}$**, reorganizing the total functional into:

$$\boxed{\mathcal{Q}_{\mathrm{total}} = \text{manifestly positive quadratic terms} + \text{prime-power correction terms}.}$$

Once this pairing is achieved, the pole terms no longer act as an independent hazard, and the prime-power terms $\mathcal{Q}_{\mathrm{prime}}$ remain as the sole structured correction to be bounded.

### Advanced Analytical Route: Completely Monotone Representations & Physical-Space Transfer (DLMF Formulation)
Beyond the discrete Weierstrass series, the digamma function $\psi(z)$ admits the classical Binet / Gauss integral representation (DLMF §5.9.12):
$$\psi(z) = \log z + \int_0^\infty \left( \frac{1}{t} - \frac{1}{1 - e^{-t}} \right) e^{-t z} \, dt \qquad (\operatorname{Re} z > 0).$$
Setting $z = \frac{1}{4} + \frac{i r}{2}$ and taking the real part:
$$h_+(r) = \operatorname{Re}\psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi = \log\left|\frac{1}{4} + \frac{ir}{2}\right| - \log \pi + \int_0^\infty \left( \frac{1}{t} - \frac{1}{1 - e^{-t}} \right) e^{-t/4} \cos\left(\frac{r t}{2}\right) \, dt.$$
The integrand factor $\left(\frac{1}{1 - e^{-t}} - \frac{1}{t}\right)$ is strictly positive and completely monotone on $(0, \infty)$.  
This representation decomposes $h_+(r)$ into an asymptotic positive logarithm plus an oscillatory cosine transform of a completely monotone density.

If this representation is split as:
$$h_+(r) = h_{\mathrm{positive}}(r) - h_{\mathrm{elementary}}(r),$$
where $h_{\mathrm{elementary}}(r)$ captures the compact negative well on $[0, r_*]$ and possesses an explicit, tractable inverse Fourier transform, then the integral:
$$\frac{1}{\pi} \int_0^\infty h_{\mathrm{elementary}}(r) K_{\mathrm{Fourier}}(v, r, L) \, dr$$
can be pulled back directly into physical coordinate space $[0, L]$ via the reflected autocorrelation identity (Paper 4, Theorem 2.3):
$$K_{\mathrm{Fourier}}(v, r, L) = \frac{1}{L} \int_0^L K_v^{\mathrm{phys}}(L - y) \cos(r y) \, dy.$$
In physical coordinate space, the prime-power contributions:
$$\mathcal{Q}_{\mathrm{prime}}(v) = -\frac{1}{L} \sum_{p^k \le c} \frac{\log p}{p^{k/2}} K_v^{\mathrm{phys}}(L - \log p^k)$$
are naturally localized Dirac evaluations at points $y = \log p^k$.  
Mapping the negative Archimedean piece into physical space would allow a **pointwise, spatial comparison** between the negative Archimedean density and the discrete prime delta masses, completely bypassing the oscillatory Gibbs phenomenon of the Fourier transform.

### Route C: Integrated Bivariate Kernel & Dual Lattice Sampling (Grand Unified Theory)
*Concept:* Formulate the complete Weil quadratic form as a single bivariate Fredholm integral on $[0, L] \times [0, L]$:
$$\mathcal{Q}_{\mathrm{Weil}}(v) = \iint_{[0, L]^2} T_v(x) \, \mathcal{K}_{\mathrm{Weil}}(x, y) \, T_v(y) \, dx \, dy.$$
Then Weil positivity becomes the positive semi-definiteness of the Fredholm kernel $\mathcal{K}_{\mathrm{Weil}} \succeq 0$ on the appropriate function space.

#### The Dual Lattice Sampling Mechanism: Fourier Spectral Lattice vs Arithmetic Prime Lattice
A crucial structural discovery in Paper 4 is the exact lattice sampling property of the Fourier kernel:
$$K_{\mathrm{Fourier}}(v, a_m, L) = \frac{L}{2} v_m^2 \qquad \text{for } a_m = \frac{2\pi m}{L}.$$
Thus, the Fourier-side kernel acts as an orthogonal sampling functional on the **Fourier spectral lattice**:
$$\Lambda_{\mathrm{Fourier}} = \left\{ a_m = \frac{2\pi m}{L} : m \in \mathbb{Z} \right\}.$$
Meanwhile, the prime contribution in Weil's explicit formula:
$$\mathcal{Q}_{\mathrm{prime}}(v) = -\sum_{p, k} \frac{\log p}{p^{k/2}} \big[ F_v(\log p^k) + F_v(-\log p^k) \big]$$
is built entirely from discrete evaluations at the **arithmetic prime lattice**:
$$\Lambda_{\mathrm{arith}} = \left\{ \log(p^k) : p \text{ prime}, k \ge 1 \right\}.$$
Therefore, the entire Weil positivity problem is fundamentally an inquiry into the **interaction between these two incommensurate discrete sampling structures**:
- The periodic Fourier lattice $\Lambda_{\mathrm{Fourier}}$ (generated by the box size $L = \log c$).
- The rigid arithmetic lattice $\Lambda_{\mathrm{arith}}$ (generated by the primes $\log p$).

Because the spatial wave $T_v(x)$ is an entire function of exponential type $L/2$, both sampling functionals act on the same Paley–Wiener / de Branges space. Weil positivity is thereby transformed from an intractable general quadratic-form inequality into a statement about a positive discrete spectral measure minus an arithmetic discrete measure acting on entire functions of exponential type. This connects the Connes–CvS Galerkin truncation directly to the axiomatic positivity theory of Louis de Branges and Arne Beurling.

---

## 5. The Two Distinct Research Programmes

To maintain absolute mathematical discipline, the project strictly bifurcates the overall objective into two independent programmes:

```
========================================================================================
PROGRAMME 1: FINITE-RANK GALERKIN POSITIVITY
Goal: Prove that for every finite prime cutoff c > 1 and every band N >= 1,
      the Galerkin matrix is strictly positive definite:
      Q_{c, N} >= 0.
Method: Finite-dimensional operator algebra, resolvent pairings, generalized eigenvalue bounds.
========================================================================================
                                     |
                                     | [Requires independent density theorem]
                                     v
========================================================================================
PROGRAMME 2: CONTINUUM DENSITY AND THE WEIL CRITERION
Goal: Prove that the union of finite Galerkin test spaces is dense in Weil's
      admissible test class on the idele class group:
      union_{c, N} H_{c, N} dense in W_admissible  ==>  Weil Positivity  ==>  RH.
Method: Semiclassical continuum limits, Sobolev regularizations, adelic harmonic analysis.
========================================================================================
```

Paper 4 and Paper 4B provide the exact finite-rank toolkit and the analytical programme for Programme 1. Finite-rank positivity $\mathcal{Q}_{c, N} \succeq 0$ remains an open target of the present investigation and does not require prior resolution of Programme 2.

---

## 6. Low-Mode Arithmetic Localization for $c = 13$

For the primary benchmark cutoff $c = 13$, an arithmetic coincidence severely restricts the dimension of the dangerous subspace:
$$L = \log 13 \approx 2.56495, \qquad a_m = \frac{2\pi m}{L} \approx 2.4496 \cdot m.$$

Comparing the discrete Fourier frequencies $a_m$ against the zero of $h_+(r)$ ($r_* \approx 6.28984$):
- Mode $m = 0$: $r = 0 < r_*$ ($h_+(0) \approx -5.372$, inside negative region)
- Mode $m = 1$: $a_1 \approx 2.4496 < r_*$ ($h_+(a_1) \approx -2.571$, inside negative region)
- Mode $m = 2$: $a_2 \approx 4.8992 < r_*$ ($h_+(a_2) \approx -0.738$, inside negative region)
- Mode $m = 3$: $a_3 \approx 7.3488 > r_*$ ($h_+(a_3) \approx +0.421$, **outside negative region**)
- Modes $m \ge 4$: $a_m \gg r_*$ (strictly positive, $h_+(a_m) \sim \log m$)

### The Three-Mode Sector Hypothesis
For $c = 13$, the discrete Fourier lattice frequencies satisfy $a_0 = 0 < r_*$, $a_1 < r_*$, and $a_2 < r_*$, whereas $a_3 \approx 7.3488 > r_*$. Thus, the **Fourier lattice samples** $a_0, a_1, a_2$ lie inside the negative region of $h_+(r)$, while $a_m$ for $m \ge 3$ do not.

However, the negative Archimedean form:
$$\mathcal{Q}_{\mathrm{arch}}^{(-)}(v) = \frac{1}{\pi} \int_0^{r_*} |h_+(r)| \Phi_v(r)^2 \, dr$$
integrates the squared amplitude continuously across the entire interval $r \in [0, r_*]$, rather than sampling strictly at discrete lattice points. Every basis amplitude $\Phi_m(r)$ has non-trivial support throughout $0 < r < r_*$. Consequently, modes $m \ge 3$ can and generally do participate in the negative operator through **off-lattice cross terms**.

The proposed block decomposition:
$$\mathbb{R}^{N+1} = \mathcal{V}_{\mathrm{low}} \oplus \mathcal{V}_{\mathrm{high}}, \qquad \mathcal{V}_{\mathrm{low}} = \operatorname{span}\{e_0, e_1, e_2\}, \quad \mathcal{V}_{\mathrm{high}} = \operatorname{span}\{e_3, \dots, e_N\},$$
is therefore a **natural structural and numerical ansatz to investigate**, rather than an established mathematical reduction to a $3 \times 3$ problem:
$$\mathcal{Q}_{\mathrm{Weil}} = \begin{pmatrix} \mathcal{Q}_{\mathrm{low}} & \mathcal{Q}_{\mathrm{cross}} \\ \mathcal{Q}_{\mathrm{cross}}^T & \mathcal{Q}_{\mathrm{high}} \end{pmatrix}.$$

The lattice sampling suggests a natural three-mode low-frequency sector. The extent to which the negative Archimedean operator is effectively confined to, or dominated by, this sector is an open question to be tested directly by Cell 63. The generalized eigenvectors of $(\mathcal{Q}_{\mathrm{positive}}, \mathcal{Q}_{\mathrm{arch}}^{(-)})$ will reveal whether the dangerous subspace actually collapses onto $\operatorname{span}\{v_0, v_1, v_2\}$ or requires broader modal participation.

---

## 7. Concrete Computational Suites: Reconnaissance (Cell 63) and High-Precision Schur Decoupling (Cell 64)

Before attempting formal proofs, we execute targeted computational suites to inspect the negative operator directly and verify finite-rank positivity at high precision:

### Objectives & Mathematical Protocol of `cell63.py` (Reconnaissance)
1. **Explicit Matrix Construction of $\mathcal{Q}_{\mathrm{arch}}^{(-)}$:**
   Compute the $(N+1) \times (N+1)$ positive semi-definite Gram matrix:
   $$\big[\mathcal{Q}_{\mathrm{arch}}^{(-)}\big]_{mn} = \frac{1}{\pi} \int_0^{r_*} |h_+(r)| \Phi_m(r) \Phi_n(r) \, dr,$$
   in the canonical basis $v$, where $\Phi_m(r)$ are the canonical Fourier basis amplitudes.
2. **Explicit Construction of $\mathcal{Q}_{\mathrm{positive\ side}}$:**
   Compute $\mathcal{Q}_{\mathrm{positive}} = \mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}}^{(+)}$, where:
   $$\mathcal{Q}_{\mathrm{arch}}^{(+)} = \mathcal{Q}_{\mathrm{arch}} + \mathcal{Q}_{\mathrm{arch}}^{(-)}.$$
3. **Gram Matrix Definiteness & Nullspace Audit:**
   Establish whether $\mathcal{Q}_{\mathrm{arch}}^{(-)}$ is strictly positive definite on the finite Galerkin basis or possesses an effective numerical nullspace. If $\mathcal{Q}_{\mathrm{arch}}^{(-)} \succ 0$, the generalized eigenproblem is standard. If $\mathcal{Q}_{\mathrm{arch}}^{(-)}$ has near-zero singular values, restrict the generalized problem to $(\ker \mathcal{Q}_{\mathrm{arch}}^{(-)})^\perp$ while verifying independent positivity $\mathcal{Q}_{\mathrm{positive}} \succ 0$ on $\ker \mathcal{Q}_{\mathrm{arch}}^{(-)}$.
4. **Generalized Spectral Decomposition:**
   Solve the generalized eigenvalue problem:
   $$\mathcal{Q}_{\mathrm{positive}} x = \lambda \mathcal{Q}_{\mathrm{arch}}^{(-)} x.$$
   - **Numerical Dominance Test:** If the computed generalized spectrum satisfies $\lambda_{\min} > 1$ with a margin substantially larger than numerical error, this provides strong finite-$N$ numerical evidence for $\mathcal{Q}_{\mathrm{Weil}} \succeq 0$. A rigorous certificate would require certified eigenvalue / interval error bounds (planned as a subsequent Cell 64).
   - **Dangerous State Identification:** The generalized eigenvector $x_{\min}$ corresponding to $\lambda_{\min}$ defines the exact profile of the "most dangerous test vector" challenging Weil positivity.
5. **Three Core Diagnostic Suites:**
   - **Diagnostic 1 (Spectrum & Effective Rank of $\mathcal{Q}_-$):** Tabulate the full spectrum of $\mathcal{Q}_{\mathrm{arch}}^{(-)}$ to determine its conditioning and effective dimensional rank across $N \in \{4, 8, 12, 16, 20, 24\}$.
   - **Diagnostic 2 (Modal Projection of Dangerous State):** Compute the fractional energy projection of $x_{\min}$ onto the 3-mode sector $\operatorname{span}\{e_0, e_1, e_2\}$ vs. its orthogonal complement $\operatorname{span}\{e_3, \dots, e_N\}$ to test the Three-Mode Sector Hypothesis directly.
   - **Diagnostic 3 (Schur-Complement Test for High Modes):** Partition the full Weil matrix as $\mathcal{Q}_{\mathrm{Weil}} = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$, where $A$ is $3 \times 3$ (low modes) and $C$ is $(N-2) \times (N-2)$ (high modes). Verify whether $C \succ 0$ and evaluate the spectrum of the Schur complement:
     $$S_{\mathrm{low}} = A - B C^{-1} B^T$$
     to determine whether high modes can be eliminated while rigorously preserving positivity.

---

### 7.1 Key Findings of Cell 63 Reconnaissance (`cell63.out`)

The execution of `cell63.py` across $N \in \{4, 8, 12, 16, 20, 24\}$ at 50-digit precision has produced a decisive structural diagnostic of the finite-dimensional operator dominance problem:

1. **Severe Ill-Conditioning of the Gram Matrix $\mathcal{Q}_-^{(N)}$:**
   - As a Gram matrix of Fourier basis functions integrated over the fixed compact interval $[0, r_*]$, $\mathcal{Q}_-^{(N)}$ develops rapidly collapsing singular values:
     $$\sigma_{\min}(\mathcal{Q}_-) = 4.03 \times 10^{-5} \,(N=4), \quad 1.16 \times 10^{-18} \,(N=8), \quad 7.63 \times 10^{-36} \,(N=12),$$
     falling below the 50-digit precision floor ($\approx -3.5 \times 10^{-51}$) for $N \ge 16$.
   - At 50-digit precision, the effective numerical rank saturates at **13** for $N \ge 12$, with the remaining directions falling below the numerical resolution threshold (the apparent nullspace is consistent with finite-precision resolution rather than an exact mathematical kernel).
2. **Breakdown of Generalized Eigenvalue Whitening at $N \ge 12$:**
   - The standard whitening transformation $\mathcal{Q}_-^{-1/2} \mathcal{Q}_{\mathrm{pos}} \mathcal{Q}_-^{-1/2}$ violently magnifies floating-point noise along the vanishing singular directions of $\mathcal{Q}_-$.
   - The eigenpair residual $\|\mathcal{Q}_{\mathrm{pos}} x_{\min} - \lambda_{\min} \mathcal{Q}_- x_{\min}\|_2$ degrades by 40 orders of magnitude:
     $$3.8 \times 10^{-51} \,(N=4), \quad 9.8 \times 10^{-51} \,(N=8) \quad \longrightarrow \quad 1.1 \times 10^{-21} \,(N=12) \quad \longrightarrow \quad 1.9 \times 10^{-11} \,(N=24).$$
   - Consequently, the apparent dominance margins $\mu_{\min} = \lambda_{\min} - 1 \sim 10^{-22}$ at $N \ge 12$ lie within the numerical error cloud of the whitening breakdown and **must not be interpreted as physical proof of operator dominance**.
3. **Striking Verification of the Three-Mode Sector Hypothesis:**
   - Despite the ill-conditioning of the generalized spectrum, the dangerous state $x_{\min}$ exhibits an extraordinary and stable modal energy concentration in the 3-mode sector $\mathcal{V}_{\mathrm{low}} = \operatorname{span}\{e_0, e_1, e_2\}$:
     $$E_{\mathrm{low}} / E_{\mathrm{total}} = 99.59\% \,(N=4), \quad 97.85\% \,(N=8), \quad 97.50\% \,(N=12), \quad 97.82\% \,(N=16), \quad 98.30\% \,(N=20), \quad 98.41\% \,(N=24).$$
   - The coordinates of $x_{\min}$ converge stably to:
     $$x_{\min} \approx (0.59, -0.69, 0.38, -0.13, 0.02, -0.002, 0.0, \dots),$$
     showing a strongly localized and apparently rapidly decaying high-mode tail. The physical challenge to Weil positivity is empirically localized in the 3-mode sector across tested dimensions.
4. **Schur Complement High-Mode Decoupling:**
   - In the block partition $\mathcal{Q}_{\mathrm{Weil}} = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$, the high-mode block is strictly positive definite ($C \succ 0$) across all dimensions ($\lambda_{\min}(C) = 1.15 \times 10^{-3}$ at $N=4$ down to $2.01 \times 10^{-27}$ at $N=24$).
   - The $3 \times 3$ Schur complement $S_{\mathrm{low}} = A - B C^{-1} B^T$ remains strictly positive definite across all dimensions ($\lambda_{\min}(S_{\mathrm{low}}) = 8.86 \times 10^{-15}$ at $N=4$ down to $5.84 \times 10^{-36}$ at $N=24$).
   - The numerical results indicate that Schur-complement elimination is a viable and substantially better-conditioned candidate route than generalized-eigenvalue whitening, avoiding the singular inversion of $\mathcal{Q}_-$.

---

### Objectives & Mathematical Protocol of `cell64.py` (High-Precision Schur Decoupling)

Motivated by the ill-conditioning of the Gram matrix $\mathcal{Q}_-^{(N)}$ diagnosed in Cell 63, `cell64.py` formulates an alternative, numerically well-conditioned pathway to verify finite-rank Weil positivity without inverting $\mathcal{Q}_-$:
1. **$LDL^T$ Diagonal Pivot Audit at 80 dps:** Rather than performing generalized eigenvalue whitening, compute the symmetric $LDL^T$ factorization of the high-mode block $C$ and the low-mode Schur complement $S_{\mathrm{low}} = A - B C^{-1} B^T$ at 80 decimal digits of precision (`dps = 80`). Strict positivity of all diagonal pivots $D_{ii} > 0$ verifies positive definiteness of the computed matrix with backward reconstruction residual $\|M - L D L^T\|_\infty / \|M\|_\infty \le 10^{-81}$.
2. **Three-Mode Schur Decoupling Test:** Partition $\mathcal{Q}_{\mathrm{Weil}}$ into the low-frequency sector $\mathcal{V}_{\mathrm{low}} = \operatorname{span}\{e_0, e_1, e_2\}$ and high-frequency sector $\mathcal{V}_{\mathrm{high}} = \operatorname{span}\{e_3, \dots, e_N\}$. Verify $C \succ 0$ and $S_{\mathrm{low}} \succ 0$ across $N \in \{4, 8, 12, 16, 20, 24\}$, testing $\mathcal{Q}_{\mathrm{Weil}} \succ 0$ via the classical symmetric Schur complement criterion.
3. **Cutoff Sensitivity Sweep at $N=24$:** Sweep $m_{\mathrm{cut}} \in \{1, 2, 3, 4, 6, 8, 12 = N/2\}$ to test the stability of Schur elimination across the semiclassical barrier.

---

### 7.2 Key Findings of Cell 64 High-Precision Positivity Verification (`cell64.out`)

The execution of `cell64.py` across $N \in \{4, 8, 12, 16, 20, 24\}$ at 80-digit precision (`dps = 80`) has achieved high-precision computational verification of finite-rank Weil positivity via symmetric Schur complement block decoupling, completely bypassing the ill-conditioned Gram inversion:

1. **High-Precision Factorization and Backward Stability:**
   - Evaluated using an exact self-contained $LDL^T$ symmetric factorization algorithm at 80 decimal digits of precision.
   - The relative backward errors $\|M - L D L^T\|_{\infty} / \|M\|_{\infty}$ for both the high-mode block $C$ and the $3 \times 3$ Schur complement $S_{\mathrm{low}}$ remain stably bounded between $10^{-81}$ and $10^{-82}$ across all tested dimensions ($N \in \{4, 8, 12, 16, 20, 24\}$), demonstrating backward stability of the factorization over 35 decimal orders below the physical ground state $\lambda_0 \sim 10^{-43}$ at $N=24$. (This accurately measures reconstruction error of the computed matrix; rigorous enclosure of the continuous mathematical operator is reserved for subsequent interval methods).

2. **Strict Pivot Positivity Across All Dimensions:**
   - **High-Mode Block Definiteness:** All diagonal pivots $D_{ii}(C) > 0$ are strictly positive across all dimensions $N \in \{4, \dots, 24\}$, confirming $C \succ 0$. The minimum pivot $\min(D_C) \approx 1.46 \times 10^{-6}$ stabilizes asymptotically for $N \ge 12$.
   - **Low-Mode Schur Complement Definiteness:** All diagonal pivots $D_{ii}(S_{\mathrm{low}}) > 0$ are strictly positive across all dimensions $N \in \{4, \dots, 24\}$, confirming $S_{\mathrm{low}} \succ 0$.
   - **Finite-Rank Positivity Verified:** By the classical symmetric Schur complement criterion ($\mathcal{Q}_{\mathrm{Weil}} \succ 0 \iff C \succ 0 \text{ and } S_{\mathrm{low}} \succ 0$), 80-digit numerical verification of positivity for the tested finite matrices $N \in \{4, 8, 12, 16, 20, 24\}$ at $c = 13$ is achieved without inverting $\mathcal{Q}_-$ (rigorous operator enclosure deferred to interval methods).

3. **Ground-State Scale Preservation and Spectral Separation in the 3-Mode Schur Complement:**
   - The smallest eigenvalue $\lambda_{\min}(S_{\mathrm{low}})$ tracks the full matrix ground state $\lambda_0(\mathcal{Q}_{\mathrm{Weil}})$ within $0.4\%$–$5.3\%$ across 30 decimal orders of magnitude:
     $$\begin{array}{c|c|c|c|c}
     N & \lambda_0(\mathcal{Q}_{\mathrm{Weil}}) & \lambda_{\min}(S_{\mathrm{low}}) & \text{Relative Shift} & \lambda_{\min}(C) \\
     \hline
     4  & 8.8274 \times 10^{-15} & 8.8634 \times 10^{-15} & +0.41\% & 1.151 \times 10^{-3} \\
     8  & 6.7109 \times 10^{-23} & 6.8583 \times 10^{-23} & +2.20\% & 2.012 \times 10^{-10} \\
     12 & 1.7825 \times 10^{-29} & 1.8457 \times 10^{-29} & +3.55\% & 6.229 \times 10^{-16} \\
     16 & 7.1184 \times 10^{-35} & 7.4309 \times 10^{-35} & +4.39\% & 1.621 \times 10^{-20} \\
     20 & 1.3232 \times 10^{-39} & 1.3886 \times 10^{-39} & +4.94\% & 1.409 \times 10^{-24} \\
     24 & 2.5335 \times 10^{-43} & 2.6669 \times 10^{-43} & +5.26\% & 2.0126 \times 10^{-27}
     \end{array}$$
   - **Structural Significance:** After exactly integrating out all high modes $m \ge 3$, the resulting 3D effective operator $S_{\mathrm{low}} = A - B C^{-1} B^T$ captures essentially the entire exponentially small ground-state tunneling scale, while the high sector remains strongly positive relative to it ($\lambda_{\min}(C) \gg \lambda_0(\mathcal{Q}_{\mathrm{Weil}})$, being $\sim 16$ orders of magnitude larger at $N=24$). The physical near-zero direction lives squarely inside the low-mode effective operator.
   - The condition number of $S_{\mathrm{low}}$ is bounded by $\kappa(S_{\mathrm{low}}) \approx 4.22 \times 10^{13}$ at $N=24$, completely resolving the whitening conditioning collapse of Cell 63 ($\kappa(\mathcal{Q}_-) > 10^{50}$).

4. **Cutoff Hierarchy Sweep at $N = 24$:**
   - Sweeping the partition threshold $m_{\mathrm{cut}} \in \{1, 2, 3, 4, 6, 8, 12 = N/2\}$ demonstrates that all pivots remain strictly positive ($D_{ii} > 0$) for every cutoff.
   - Furthermore, once $m_{\mathrm{cut}} \ge 6$, the Schur complement eigenvalue $\min\operatorname{eig}(S)$ matches $\lambda_0(\mathcal{Q}_{\mathrm{Weil}}) = 2.53348484008 \times 10^{-43}$ to 12 significant figures, establishing a clean hierarchy:
     $$m \ge 12 \quad \longrightarrow \quad m \ge 8 \quad \longrightarrow \quad m \ge 6 \quad \text{(increasingly negligible back-reaction)}$$
     rather than the decoupling being an artifact of $m_{\mathrm{cut}} = 3$.

---

## 8. The 7-Stage Strategic Project Roadmap

```
========================================================================================
STAGE I: GLOBAL FOURIER KERNEL POSITIVITY
Status: COMPLETED (Paper 4, Theorem 4.1)
Result: K_Fourier(v, r, L) = |Phi_v(r)|^2 >= 0 proven algebraically on R.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE II: EXACT FINITE-N RESOLVENT & COMMUTATOR TOOLKIT
Status: COMPLETED (Paper 4, Theorems 3.1, 5.1, 6.1, 7.3)
Result: Closed-form Cauchy transform J(q), Weierstrass pole series, rank-2k commutators,
        and exact K^2-resolvent formula for D_1/D_0.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE III: OPERATOR DOMINANCE RECONNAISSANCE & SCHUR DECOUPLING
Status: COMPLETED (Cell 63 & Cell 64)
Results:
  1. Q_- Gram matrix condition collapse diagnosed; generalized eigenvalue whitening
     identified as ill-posed for N >= 12 (Cell 63).
  2. Three-Mode Sector Hypothesis confirmed empirically (98% modal energy in {e0, e1, e2}) (Cell 63).
  3. 80-digit numerical verification of positivity for tested finite matrices via symmetric LDL^T Schur complement
     decoupling (all pivots D_{ii}(C) > 0 and D_{ii}(S_low) > 0 verified with backward error <= 2.6e-81 across all N in {4, 8, 12, 16, 20, 24}) (Cell 64).
  4. Ground-state scale capture established: \lambda_min(S_low) tracks \lambda_0(Q_Weil) within 5%
     across 30 orders of magnitude, while \lambda_min(C) >> \lambda_0(Q_Weil) (16 orders higher at N=24).
Next Target: Stage IV (Analytical Pairing of Weierstrass Resolvents) & Stage V (Three-Mode Schur Reduction)
========================================================================================
                                     |
                                     v
========================================================================================
STAGE IV: ANALYTICAL PAIRING OF WEIERSTRASS RESOLVENTS
Status: PLANNED
Tasks:
  1. Pair the negative terms -J(q_n) against Q_pole and Q_prime.
  2. Isolate exact algebraic cancellations between discrete lattice poles and prime sums.
  3. Formulate manifestly positive quadratic sub-blocks.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE V: THREE-MODE SECTOR REDUCTION AND EFFECTIVE HAMILTONIAN ANALYSIS
Status: COMPLETED COMPUTATIONAL INVESTIGATION (Cell 65: Loewner Monotonicity & Resolvent Decoupling)
Tasks:
  1. Formulate the finite-dimensional effective Hamiltonian S_low(N) = A_N - \Sigma_low(N)
     with self-energy \Sigma_low(N) = B_N C_N^{-1} B_N^T.
     - Advance: Computed bare A and self-energy \Sigma_low(N) at 80 dps across N in {4, 8, 12, 16, 20, 24} (Cell 65).
     - Advance: Verified strict Loewner Monotonicity \Delta \Sigma(N) \succ 0 across every step at 80 dps (\min eig(\Delta \Sigma) > 0), demonstrating that \Sigma_low(N) increases monotonically in the Loewner ordering across tested dimensions and is bounded above by A (Cell 65).
     - Advance: Geometric increment collapse: ||\Sigma(N) - \Sigma(N-4)||_inf plunges from 4.0e-5 down to 1.7e-26, providing strong numerical evidence for operator monotone convergence \Sigma_low(N) \uparrow \Sigma_\infty \preceq A across tested dimensions and formulating the analytical convergence target (Cell 65).
     - Advance: High-mode continuum decoupling: modes m >= 12 above the barrier top carry only 0.014% of the self-energy norm, decaying by ~ 10^{-4} per 4 modes (falling to 4.77e-14 at m=24) (Cell 65).
  2. Prove block definiteness reduction: Q_N > 0 <=> C_k(N) > 0 and S_k(N) > 0 for fixed small k.
  3. Establish quantitative analytical bounds on the Schur self-energy ||B_N C_N^{-1} B_N^T||
     and investigate what asymptotic structure the 3x3 Schur complement inherits from the
     high-mode resolvent C_N^{-1}, connecting to the D_0, D_1 / first-jet machinery.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE VI: CONTINUUM ASYMPTOTICS & SUBSPACES DENSITY
Status: COMPLETED ANALYTICAL ADVANCE (Paper 4B Propositions 8.7–8.17) & SUBSTANTIVE HYPOTHESIS AUDIT (Cell 66 & Cell 67)
Tasks:
  1. Prove polynomial bound |D_1/D_0| <= C N^p via sector-decomposed resolvents.
     - Advance: Exact rank-two commutator identity [K, Q] = \psi d^T - d \psi^T proved (Proposition 8.8).
     - Advance: D_0-free bound on odd overlaps a_j and reduction of M_1 established (Proposition 8.8).
     - Advance: Large-N asymptotic linear coercivity W[F_N] = 2\pi(1 + \log(\log c / 2)) N + O(\log N) proved via Archimedean boundary layer and explicit prime quadratic form algebra, yielding M_d^{(1)} >= c_0 N with c_0 = \pi(1 + \log(\log c / 2)) > 0 (Proposition 8.9).
     - Advance (Part I Unconditional): Exact commutator doublet cancellation a_0 = -((\mu_0 - \lambda)/D_0) b_00, Parseval excited odd enclosure \sum_{j>=1} a_j^2 <= 2 N \Lambda_*^2, and Cauchy–Schwarz moment coercivity \langle d, R_even d \rangle >= (4/C_R) N established unconditionally (Proposition 8.10).
     - Advance (Part II Conditional): Polynomial decoupling |D_1/D_0| <= C N^p (p = \max(\eta + \gamma_e, \eta + \gamma_{\mathrm{filt}}, \max(2, q)) < \infty, with \eta = \max(2, \gamma), non-sharp) established under WKB tunneling-flux (H1), even-sector transmission growth (H2), odd-sector resolvent bound (H2_odd), and excited second-moment condition (D_0^2 M_{2,exc} <= C_2 N^q) hypotheses, with M_1 = O(N^{\max(2, \gamma)}) and D_0^2 M_2 = O(N^{\max(2, q)}) (Proposition 8.10).
     - Advance (Substantive Mechanism of H1): Semiclassical flux-matching relation R_tun(N) = (\mu_0 - \lambda) / D_0^2 \in [2.4, 6.0] confirmed across 20 decimal orders of magnitude (N in [8, 24]), providing strong numerical evidence for the semiclassical flux-matching invariant underlying H1 (Landau-Lifshitz QM §50).
     - Advance (Substantive Mechanism of H2 & H2_odd): Audited the bound-state ladder beneath the barrier top. The bare polynomial-gap formulation is not supported by the data: the relevant bare gaps collapse exponentially (E_1 - \lambda \sim 10^{-37}, \mu_2 - E_1 \sim 10^{-34} at N=24). However, mode-by-mode transmission cancellation keeps the weighted resolvent and filtering ratios O(1) individually (d_k^2 / (\mu_{k+1}-E_k) <= 6.85, d_k^2 / (E_k - \lambda) \approx 580, a_k^2 / (\mu_k - \lambda) \approx 4 - 6), distinguishing the analytic growth hypotheses (H2, H2_odd) from the underlying semiclassical mechanism conjecture.
     - Advance (Substantive Mechanism of H3): Provided strong numerical evidence for exponential boundary extinction D_0^2 \le C_0 e^{-\sigma N}, with effective decay rate \sigma(N) = -(1/N)\log(D_0^2) converging to \sigma_WKB = (\pi/2)\log c \approx 4.029 within 5.0% at N=24.
     - Advance (Conditional Second-Moment Reduction via Relative Gap Control): Proved that assuming H2_odd alongside the relative tunneling-gap condition H2_gap: D_0^2 / (\mu_1 - \lambda) <= C_gap N^{q_0}, the second-moment condition D_0^2 M_{2,exc} <= C_2 N^q holds rigorously with q = \gamma + \max(0, q_0) (Proposition 8.11). Established the exact convex combination \rho_2^{exc} \equiv D_0^2 M_{2,exc} / M_1^{exc} = \sum_{j>=1} w_j (D_0^2 / (\mu_j - \lambda)) <= R_{gap}^{max}. Formulated candidate semiclassical barrier thinning with its five open analytical bridges, identifying the diagnostic for the proposed asymptotic mechanism \Delta \sigma_N^{gap} = -(1/N)\log R_{gap}^{max} >= 0.610 across N in {8, 12, 16, 20, 24} (Cell 67).
     - Advance (Exact Coordinate Parseval Identity and Dipole Relative Gap Reduction): Proved the exact finite-$N$ Parseval coordinate identity $D_0^2 M_2 = \|Kc\|^2 = \sum_{m=1}^N m^2 v_m^2$ unconditionally (with $\|Kc\|^2 \le N^2$), sharpening to $D_0^2 M_2 \le B_\infty(c) < \infty$ conditional on uniform solitary decay $|v_{N, m}| \le C_{\mathrm{sol}} q^m$. Proved the exact transition dipole representation $b_{0j}^2 = \frac{D_0^2 a_j^2}{(\mu_j - \lambda)^2}$ and excited decomposition $D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2 = \|Kc\|^2 \sin^2 \theta(Kc, u_0)$. Established the structural reduction $R_{\mathrm{gap}} \le \|P_{\perp u_0} Kc\|^2 / [a_1^2 / (\mu_1 - \lambda)]$, providing a structural explanation for the observed suppression via the $99.99987\%$ directional alignment observed at $N = 24$ between the coordinate-derivative wavepacket $Kc$ and the lowest odd mode $u_0$. Elevated Route B (weighted relative-gap theorem $\rho_2^{\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2 / M_1^{\mathrm{exc}}$) as the primary analytical target, with Route A conditional on a uniform first-mode transmission lower bound (Proposition 8.12).
     - Advance (Exact Even-Resolvent Operator Representation, Ground-Doublet Pole Cancellation, Unified Mode Sum, and Route B Operator Bound): Proved exact finite-dimensional second-commutator identity $(Q_{\mathrm{even}} - \lambda I) K^2 c = - D_0 (K \boldsymbol\psi + M_1 d)$ on $\mathbb{R}^{2N+1}$, established exact ground-doublet pole cancellation $K v_{\mathrm{exc}} = \|v_{\mathrm{exc}}\|^2 c + D_0 \Xi$ via the First Resolvent Identity, and proved the unified mode-by-mode representation $\Xi = \sum_{k \ge 1} d_k [\sum_{j \ge 1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}] u_k$ with exact $j=0$ mode cancellation. Proved the exact unconditional Route B operator inequality $\rho_2^{\mathrm{exc}} \le \frac{D_0^2 \|\Xi\|^2}{M_1^{\mathrm{exc}}(1 - D_0^2 M_{2,\mathrm{exc}})}$ and dual variational Rayleigh quotient representation $\rho_2^{\mathrm{exc}} = \frac{D_0^2}{\mathcal{R}_{Q_{\mathrm{odd}} - \lambda}(v_{\mathrm{exc}})}$ (Proposition 8.13). This eliminates the independent excited-odd-gap hypothesis $\mathrm{H2}_{\mathrm{gap}}$ from the algebraic reduction and transfers the remaining control problem entirely to even-sector resolvent estimates.
     - Advance (Stieltjes Difference Quotient Representation and Spectral Enclosures): Proved exact Stieltjes difference quotient identity $A_k = \frac{F_{\mathrm{exc}}(E_k) - F_{\mathrm{exc}}(\lambda)}{E_k - \lambda}$ for each mode coefficient of $\Xi = \sum_{k \ge 1} d_k A_k u_k$. Proved smooth non-singular derivative enclosure for mode $k=1$ ($M_{2,\mathrm{exc}} \le A_1 \le \frac{M_1^{\mathrm{exc}}}{\mu_1 - E_1}$ on the pole-free interval $[\lambda, E_1]$), established high-mode transmission enclosure $A_k^{(> k)} \le \frac{M_1^{\mathrm{exc}}}{\mu_{k+1} - E_k}$, and identified the resonant parity-doublet transmission factor $\frac{d_k a_k}{\mu_k - E_k}$ as the remaining modewise quantity requiring analytical control in a mode-by-mode approach (Proposition 8.14).
     - Advance (Global Even-Resolvent Operator Representation, Parseval Identity, and Resonant Absorption): Proved exact global operator identity $\Xi = g(Q_{\mathrm{even}}) d_{\mathrm{exc}} = \int \frac{d\nu_{\mathrm{exc}}(\mu)}{\mu - \lambda} (\mu I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}} = M_1^{\mathrm{exc}} \sum_{j \ge 1} w_j (\mu_j I - Q_{\mathrm{even}})^{-1} d_{\mathrm{exc}}$. Established exact finite-$N$ Parseval even-resolvent identity $a_j^2 \sum_{k=1}^{N-1} \frac{d_k^2}{(\mu_j - E_k)^2} = \|K u_j\|^2 - b_{0j}^2 \le N^2$ for every excited odd mode without eigenvalue denominators (with $a_j \ne 0$ proved by strict spectral interlacing). Proved second-generation coordinate-derivative wavepacket geometry $\Xi = \frac{1}{D_0} P_{\perp c} K v_{\mathrm{exc}} \implies \|\Xi\|^2 = \frac{1}{D_0^2} \|P_{\perp c} K v_{\mathrm{exc}}\|^2$. Established global Stieltjes divided-difference kernel representation $\|\Xi\|^2 = \iint [\frac{G_{\mathrm{even}}(\mu') - G_{\mathrm{even}}(\mu)}{\mu - \mu'}] \frac{d\nu_{\mathrm{exc}}(\mu)}{\mu - \lambda} \frac{d\nu_{\mathrm{exc}}(\mu')}{\mu' - \lambda}$ and closed-form upper bound $\|\Xi\|^2 \le M_1^{\mathrm{exc}} \sum_{j \ge 1} \frac{\|K u_j\|^2 - b_{0j}^2}{\mu_j - \lambda}$, automatically absorbing the resonant contribution and reducing the second-moment bound directly to $\rho_2^{\mathrm{exc}} \le \frac{N^3}{1 - D_0^2 M_{2,\mathrm{exc}}} R_{\mathrm{gap}}^{\max}(N)$ without modewise $\mathcal{O}(1)$ assumptions (Proposition 8.15).
     - Advance (Weighted Coordinate-Energy Measure, Exact Trace Enclosure, and High-Energy Decoupling): Proved exact odd coordinate trace identity $\sum_{j=0}^{N-1} \|K u_j\|^2 = \operatorname{Tr}_{H_{\mathrm{odd}}}(K^2) = \frac{N(N+1)(2N+1)}{6} = \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6}$, with excited bound $\sum_{j=1}^{N-1} \|K u_j\|^2 < \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6} = \mathcal{O}(N^3)$. Defined the coordinate-energy spectral measure $d\eta_{\mathrm{coord}}(\mu) = \sum_{j \ge 1} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2 \delta_{\mu_j}$ and proved that the exact weighted coordinate-energy sum evaluates to the Stieltjes transform $\mathcal{S}_{\mathrm{coord}}(N) = D_0^2 \int \frac{d\eta_{\mathrm{coord}}(\mu)}{\mu - \lambda}$. Proved that the high-energy (above-barrier) sector $\mu_j \ge V_* > \sup_N \lambda_N$ decouples exponentially fast as $\mathcal{O}(e^{-\sigma N} N^3) \to 0$ under Hypothesis H3 and uniform barrier separation $c_{\mathrm{high}} > 0$ (unconditional with respect to relative-gap hypotheses), isolating the active analytical target strictly to the finite bound-state ladder (Proposition 8.16).
     - Advance (Exact Operator Trace Identity, Bound-State Spectral Decomposition, and Semiclassical Localization Bridges B1 & B2): Proved exact operator trace identity on $H_{\mathrm{odd}}$: $\mathcal{S}_{\mathrm{coord}}(N) = D_0^2 \operatorname{Tr}_{H_{\mathrm{odd}}}[(Q_{\mathrm{odd}} - \lambda I)^{-1} P_{\perp u_0} K P_{\perp u_0^{\mathrm{even}}} K] \le D_0^2 \sum_{m=1}^N m^2 [R_{\mathrm{odd},\perp}(\lambda)]_{mm}$, and established exact orthogonal spectral decomposition $\mathcal{S}_{\mathrm{coord}} = \mathcal{S}_{\mathrm{bound}} + \mathcal{S}_{\mathrm{high}}$ with exponential high-energy decoupling $\mathcal{S}_{\mathrm{high}} = \mathcal{O}(e^{-\sigma N} N^3) \to 0$ and finite-$N$ bound-state inequality $\mathcal{S}_{\mathrm{bound}} \le R_{\mathrm{gap}}^{\max}(N) \sum_{j=1}^{N_{\mathrm{bound}}} \|K u_j\|^2$ unconditionally (Part I). Formulated the conditional continuum reduction under explicit Bridges B1 (uniform bound-state cardinality $\bar{N}_{\mathrm{bound}} < \infty$) and B2 (uniform discrete coordinate kinetic energy $\sup_N \|K u_j^{(N)}\|^2 \le \mathcal{C}_j < \infty$), proving $\mathcal{S}_{\mathrm{bound}} \le \mathcal{C}_{\mathrm{bound}} R_{\mathrm{gap}}^{\max}(N)$ and establishing that $R_{\mathrm{gap}}^{\max}(N) \to 0$ is sufficient for continuum decoupling $\rho_2^{\mathrm{exc}} \to 0$ (with $R_{\mathrm{gap}}^{\max} = \mathcal{O}(1)$ yielding uniform boundedness), eliminating the polynomial $N^3$ rate requirement (Proposition 8.17).
     - Advance (Global Resolvent Commutator Identity, Weighted Coordinate-Trace Duality, and Barrier-Thinning Quenching): Proved exact resolvent commutator identity $[K, R(z)] = - R(z) [K, Q] R(z) = - (R(z)\boldsymbol\psi)(d^T R(z)) + (R(z)d)(\boldsymbol\psi^T R(z))$ and evaluated the nested commutator $[K, [K, R(z)]]$ to establish the exact parity trace difference $\operatorname{Tr}_{H_{\mathrm{odd}}}[K^2 R_{\mathrm{odd}}(z)] - \operatorname{Tr}_{H_{\mathrm{even}}}[K^2 R_{\mathrm{even}}(z)] = \langle R_{\mathrm{odd}}(z)\boldsymbol\psi, K R_{\mathrm{even}}(z) d \rangle$ unconditionally (Proposition 8.18). Proved the exact equivalence between the Fourier trace $\sum_{m=1}^N m^2 [R_{\mathrm{odd},\perp}]_{mm}$ and the coordinate-energy sum $\sum_{j=1}^{N-1} \frac{\|Ku_j\|^2}{\mu_j - \lambda} = \mathcal{T}_N$. Refuted the bare polynomial trace hypothesis ($\mathcal{T}_N \le C N^p$) by showing that the first excited odd mode tunneling gap $\mu_1 - \lambda \sim e^{-\sigma_1 N}$ forces the bare trace to diverge exponentially: $\mathcal{T}_N \ge \mathcal{C}_1 e^{+\sigma_1 N} \to \infty$. Proved that conditional on Bridges B1 and B2 and the barrier-thinning hypothesis $\sigma_0 > \sigma_1$ ($\Delta\sigma \equiv \sigma_0 - \sigma_1 > 0$), the physical scaled trace is exponentially quenched: $D_0^2 \mathcal{T}_N \le \mathcal{C}_{\mathrm{bound}} R_{\mathrm{gap}}^{\max}(N) + \mathcal{O}(e^{-2\sigma_0 N} N^3) \sim \mathcal{O}(e^{-(\sigma_0 + \Delta\sigma)N}) \to 0$, establishing that the relative tunneling gap $R_{\mathrm{gap}}^{\max} = \frac{D_0^2}{\mu_1 - \lambda}$ is the true invariant       - Advance (Exact Dipole Factorization of the First Relative Tunneling Gap and Even Resolvent Cancellation): Proved the exact finite-$N$ algebraic factorization $R_{\mathrm{gap}}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda} = \frac{b_{01}^2}{\mathcal{R}_1}$ unconditionally (Proposition 8.19), with numerical values $\mathcal{R}_1 \in [3.91, 8.00]$ consistent with an $\mathcal{O}(1)$ limit across tested dimensions. Established coordinate wavepacket upper enclosure $b_{01}^2 \le \|P_{\perp u_0} Kc\|^2 = \|v_{\mathrm{exc}}\|^2 \implies R_{\mathrm{gap}}(N) \le \frac{\|v_{\mathrm{exc}}\|^2}{\mathcal{R}_1}$, and proved the exact even-resolvent cancellation identity $b_{01}^2 = \|Ku_1\|^2 - a_1^2 \sum_{k=1}^N \frac{d_k^2}{(\mu_1 - E_k)^2}$, showing that the asymptotic equivalence $R_{\mathrm{gap}}(N) \to 0 \iff b_{01}^2 \to 0$ is conditional on $\mathcal{R}_1$ remaining bounded away from zero.
      - Advance (Exact Even-Resolvent Representation of $K u_1$, Stieltjes Derivative Architecture, and Algebraic Resolution of the Dipole Cancellation): Proved the exact finite-$N$ resolvent operator representation $K u_1 = a_1 (Q_{\mathrm{even}} - \mu_1 I)^{-1} d$ on $H_{\mathrm{even}}$ unconditionally, exploiting $\mu_1 \notin \sigma(Q_{\mathrm{even}})$ from strict parity interlacing (Proposition 8.20). Proved the exact Stieltjes derivative identity $\|K u_1\|^2 = a_1^2 G_d'(\mu_1) = a_1^2 \sum_{k=0}^N \frac{d_k^2}{(\mu_1 - E_k)^2}$. Proved that the observed cancellation $\|K u_1\|^2 - \mathcal{E}_{\mathrm{even}, 1} = b_{01}^2$ is identically the exact algebraic projection of the even resolvent vector isolating the $k=0$ ground-state pole $E_0 = \lambda$, representing an exact Pythagorean decomposition on $H_{\mathrm{even}}$: $\|(Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2 = \|P_c (Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2 + \|P_{\perp c} (Q_{\mathrm{even}} - \mu_1 I)^{-1} d\|^2$. Established that $b_{01} = - \langle P_{\perp u_0} K c, u_1 \rangle$ is suppressed because $Kc$ is aligned with the ground state $u_0$ to $> 99.99987\%$, so that the excited wavepacket residual $v_{\mathrm{exc}} = P_{\perp u_0} K c$ has norm $\|v_{\mathrm{exc}}\|^2 \sim 10^{-6}$, carrying $99.9995\%$ of its weight in the first excited dipole mode $u_1$.
    2. Control the bound-state contribution to the weighted coordinate-energy measure $\mathcal{S}_{\mathrm{bound}}(N) = \sum_{j=1}^{N_{\mathrm{bound}}} \frac{D_0^2}{\mu_j - \lambda} \|P_{\perp u_0^{\mathrm{even}}} K u_j\|^2$; establish Bridges B1 (uniform cardinality) and B2 (uniform coordinate kinetic energy) from the discrete Galerkin matrices (Milestone M12).
    3. Establish exponential boundary-defect decoupling D(N) -> 0.
       - Advance: D(N) <= C' e^{-\sigma N} (1 + N^p)^2 = O(e^{-\sigma N} N^{2p}) ---> 0 conditionally established under Hypotheses H1–H3 (with H2_odd) for T >= 1, demonstrating exponential boundary decoupling (Proposition 8.10).
    4. Formulate formal WKB double-well potential and prove density of union_{c, N} H_{c, N}.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE VII: INVOCATION OF THE WEIL CRITERION
Status: LONG-TERM OBJECTIVE
Tasks:
  1. Combine finite-rank positivity with continuum density.
  2. Deduce unconditional Weil positivity W(g) >= 0 on the idele class group.
  3. Conclude the Riemann Hypothesis.
========================================================================================
```

---

## 9. Operational Milestones

| Milestone | Action Item | Target Artifact / Script | Status / Deliverable |
| :---: | :--- | :--- | :--- |
| **M1** | Implement the finite-band negative operator and generalized eigenvalue suite | `cell63.py` | **COMPLETED** (`cell63.out`: Gram conditioning & modal localization) |
| **M2** | Audit the generalized spectrum conditioning across $N \in \{4, 8, 12, 16, 20, 24\}$ | Analytical Review | **COMPLETED** (Diagnosed whitening breakdown at $N \ge 12$; residual loss) |
| **M3** | Analyze coordinates and modal energy of the dangerous vector $x_{\min}$ | Diagnostic Report | **COMPLETED** (Confirmed $98\%$ energy in $\{e_0, e_1, e_2\}$ across all $N$) |
| **M4** | Formulate high-precision positivity suite via Schur complement block decoupling | `cell64.py` | **COMPLETED** (`cell64.out`: High-precision numerical verification via $LDL^T$ pivots $D_{ii}(C) > 0, D_{ii}(S_{\mathrm{low}}) > 0$ at 80 dps across all $N$; backward error $\le 2.6 \times 10^{-81}$; ground-state scale capture within 5%) |
| **M5** | Algebraically pair $J(q_n)$ with the pole and prime representations | Paper 4B Section Update | Exact positive block formulation (Stage IV) |
| **M6** | Lower-bound the Dirichlet kernel Weil functional $\mathcal{W}[F_N] = \langle d, Q d \rangle$ | Paper 4B Proposition 8.9 | **COMPLETED** (Closed-form boundary layer $\mathcal{W}_{\mathrm{arch}} = \mathcal{C}_{\mathbb{R}} N$, $\mathcal{W}_{\mathrm{prime}} = \mathcal{O}(\log N)$, $\mathcal{W}_{\mathrm{pole}} = \mathcal{O}(1)$, proving linear coercivity $\mathcal{M}_d^{(1)} \ge c_0 N$ with $c_0 = \pi(1 + \log\frac{\log c}{2}) > 0$) |
| **M7** | Bound the overlap-weighted odd resolvent moment $M_1 = \sum_j \frac{a_j^2}{\mu_j - \lambda}$ and prove polynomial decoupling $|D_1/D_0| \le C N^p$ | Paper 4B Proposition 8.10 | **CONDITIONAL PROPOSITION ESTABLISHED ALGEBRAICALLY; ASYMPTOTIC HYPOTHESES REMAIN UNPROVED** (Part I unconditional finite-$N$ doublet cancellation $\frac{a_0^2}{\mu_0 - \lambda} = \frac{\mu_0 - \lambda}{D_0^2} b_{00}^2$, Parseval excited bound, and Cauchy–Schwarz coercivity $\langle d, R_{\mathrm{even}} d \rangle \ge \frac{4}{\mathcal{C}_{\mathbb{R}}} N$; Part II polynomial control $|D_1/D_0| \le C N^p$ with $p = \max(\eta + \gamma_e, \eta + \gamma_{\mathrm{filt}}, \max(2, q))$ under Hypotheses H1, H2, H2$_{\mathrm{odd}}$, and technical second-moment condition $D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$ with $M_1 = \mathcal{O}(N^{\max(2, \gamma)})$ and $D_0^2 M_2 = \mathcal{O}(N^{\max(2, q)})$; exponential boundary decoupling $\mathcal{D}(N) = \mathcal{O}(e^{-\sigma N} N^{2p}) \to 0$ adding Hypothesis H3) |
| **M8** | Audit the substantive mechanisms of Hypotheses H1, H2, H2$_{\mathrm{odd}}$, and H3 (flux matching & transmission cancellation) | `cell66.py` | **COMPLETED** (`cell66.out`: verified H1 ratio $\mathcal{R}_{\mathrm{tun}} \in [2.41, 5.92]$, demonstrated exponential collapse of bare gaps refuted naive H2, confirmed mode-by-mode transmission cancellation $\frac{d_k^2}{\mu_{k+1}-E_k} \le \mathcal{O}(1)$, $\frac{a_k^2}{\mu_k-\lambda} \le \mathcal{O}(1)$, $\Sigma_{\mathrm{filt}}(24) = 94.66$, and $\sigma(N) \to \sigma_{\mathrm{WKB}} \approx 4.029$ within 5.0%) |
| **M9** | Audit relative tunneling-scale ordering, barrier thinning, and excited second resolvent moment $D_0^2 M_{2,\mathrm{exc}}$ | Paper 4B Proposition 8.11 & `cell67.py` | **COMPLETED** (`cell67.py`: audited $\mathrm{H2}_{\mathrm{gap}}$ via $R_{\mathrm{gap}}^{\max} = \frac{D_0^2}{\mu_1 - \lambda} \le 4.68 \times 10^{-6}$ decaying to $4.39 \times 10^{-7}$ at $N=24$; confirmed gap decay action $\Delta \sigma_N^{\mathrm{gap}} \ge 0.610$; proved convex combination $\rho_2^{\mathrm{exc}} = \sum_j w_j r_j \le R_{\mathrm{gap}}^{\max}$; formulated conditional reduction $\mathrm{H2}_{\mathrm{odd}} + \mathrm{H2}_{\mathrm{gap}} \implies D_0^2 M_{2,\mathrm{exc}} \le C_2 N^q$ with candidate barrier thinning mechanism) |
| **M10** | Formulate exact coordinate Parseval identity $D_0^2 M_2 = \|Kc\|^2$, excited dipole decomposition $D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2$, and reduce relative tunneling gap suppression to coordinate projection residual | Paper 4B Proposition 8.12 | **COMPLETED UNCONDITIONALLY FOR IDENTITIES; CONDITIONAL ON SOLITARY BOUND FOR $\mathcal{O}(1)$** (Exact finite-$N$ identities $D_0^2 M_2 = \|Kc\|^2$ and $D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2$ proved unconditionally; $\|Kc\|^2 \le B_\infty < \infty$ proved conditional on uniform solitary decay; relative gap reduced to $R_{\mathrm{gap}} \le \|P_{\perp u_0} Kc\|^2 / [a_1^2/(\mu_1-\lambda)]$, structurally connecting suppression to observed $99.99987\%$ alignment at $N=24$; Route B elevated as primary analytical target) |
| **M11** | Exact even-resolvent representation of odd coordinate residual, ground-doublet pole cancellation, and Route B operator bound | Paper 4B Proposition 8.13 | **COMPLETED — EXACT FINITE-$N$ ALGEBRAIC REDUCTION** (Exact finite-dimensional second-commutator identity $(Q_{\mathrm{even}} - \lambda I) K^2 c = - D_0 (K \boldsymbol\psi + M_1 d)$, ground-doublet pole cancellation via First Resolvent Identity, and unified mode-by-mode representation $\Xi = \sum_{k \ge 1} d_k [\sum_{j \ge 1} \frac{a_j^2}{(\mu_j - \lambda)(\mu_j - E_k)}] u_k$ established; exact master operator inequality $\rho_2^{\mathrm{exc}} \le \frac{D_0^2 \|\Xi\|^2}{M_1^{\mathrm{exc}}(1 - D_0^2 M_{2,\mathrm{exc}})}$ and dual variational Rayleigh quotient representation proved; eliminates reliance on excited odd tunneling gaps by transferring the remaining bound to even-sector resolvent estimates) |
| **M12** | Control the bound-state contribution to the weighted coordinate-energy measure $\mathcal{S}_{\mathrm{coord}}$ and establish whether its suppression follows from coordinate-energy localization rather than the worst-case relative gap | Paper 4B Propositions 8.14, 8.15, 8.16 & 8.17 | **ADVANCED — BOUND-STATE ENCLOSURE & EXACT OPERATOR TRACE ESTABLISHED; BRIDGES B1 & B2 IDENTIFIED** (Proved exact Stieltjes difference quotient identity $A_k = \frac{F_{\mathrm{exc}}(E_k) - F_{\mathrm{exc}}(\lambda)}{E_k - \lambda}$, exact Parseval even-resolvent identity $a_j^2 \sum_{k=1}^{N-1} \frac{d_k^2}{(\mu_j - E_k)^2} = \|K u_j\|^2 - b_{0j}^2 \le N^2$, and global closed-form quadratic enclosure $\|\Xi\|^2 \le M_1^{\mathrm{exc}} \sum_{j=1}^{N-1} \frac{\|K u_j\|^2 - b_{0j}^2}{\mu_j - \lambda}$ automatically absorbing the resonant contribution. Proved exact odd coordinate trace identity $\sum_{j=0}^{N-1} \|K u_j\|^2 = \frac{N(N+1)(2N+1)}{6} = \frac{N^3}{3} + \frac{N^2}{2} + \frac{N}{6}$, and established exponential decoupling of the high-energy (above-barrier) sector $\mathcal{S}_{\mathrm{high}} = \mathcal{O}(e^{-\sigma N} N^3) \to 0$ under H3 and uniform separation $c_{\mathrm{high}} > 0$. Proved exact operator trace identity $\mathcal{S}_{\mathrm{coord}} = D_0^2 \operatorname{Tr}[(Q_{\mathrm{odd}}-\lambda)^{-1} P_{\perp u_0} K P_{\perp u_0^{\mathrm{even}}} K] \le D_0^2 \sum m^2 [R_{\mathrm{odd},\perp}]_{mm}$, and exact bound-state inequality $\mathcal{S}_{\mathrm{bound}} \le R_{\mathrm{gap}}^{\max} \sum_{j=1}^{N_{\mathrm{bound}}} \|Ku_j\|^2$ unconditionally. Identified the remaining analytical agenda as proving from the discrete Galerkin operators: (1) Bridge B1: uniform bound-state cardinality $\bar{N}_{\mathrm{bound}} < \infty$, and (2) Bridge B2: uniform discrete coordinate kinetic energy $\sup_N \|Ku_j^{(N)}\|^2 < \infty$. Under B1 and B2, proved $\mathcal{S}_{\mathrm{bound}} \le \mathcal{C}_{\mathrm{bound}} R_{\mathrm{gap}}^{\max}(N)$, establishing that $R_{\mathrm{gap}}^{\max}(N) \to 0$ is sufficient for continuum decoupling $\rho_2^{\mathrm{exc}} \to 0$, while $R_{\mathrm{gap}}^{\max} = \mathcal{O}(1)$ yields uniform boundedness.) |
| **M13** | Audit global weighted resolvent trace $\mathcal{T}_N$, coordinate kinetic energies, and exponential quenching $D_0^2 \mathcal{T}_N$ | Paper 4B Proposition 8.18 & `cell68.py` | **COMPLETED AS DIAGNOSTIC — POLYNOMIAL TRACE REFUTED** (`cell68.py`: verified exact equality $\mathcal{T}_N^{\mathrm{trace}} = \mathcal{T}_N^{\mathrm{spec}}$ to $10^{-45}$; proved bare trace $\mathcal{T}_N$ is exponentially divergent as $e^{+\sigma_1 N}$, mathematically refuting the bare polynomial trace hypothesis; verified discrete coordinate kinetic energies $\|Ku_j\|^2 = \mathcal{O}(1)$ for bound states and total trace identity $\sum \|Ku_j\|^2 = \frac{N(N+1)(2N+1)}{6}$; proved exact parity trace difference $\operatorname{Tr}_{H_{\mathrm{odd}}}(K^2 R_{\mathrm{odd}}) - \operatorname{Tr}_{H_{\mathrm{even}}}(K^2 R_{\mathrm{even}}) = \langle R_{\mathrm{odd}}\boldsymbol\psi, K R_{\mathrm{even}} d \rangle$; confirmed that physical quenching $D_0^2 \mathcal{T}_N \to 0$ is governed by the relative tunneling gap $R_{\mathrm{gap}}^{\max}$ conditional on B1+B2 and $\sigma_0 > \sigma_1$) |
| **M14** | Exact reduction and asymptotic control of the first relative tunneling gap $R_{\mathrm{gap}} = \frac{D_0^2}{\mu_1 - \lambda}$ | Paper 4B Proposition 8.19 & `cell69.py` | **COMPLETED** (`cell69.py`: verified exact factorization $R_{\mathrm{gap}} = b_{01}^2 / \mathcal{R}_1$ and excited dipole share $> 99.99\% \to 99.9995\%$ across $N \in \{8, \dots, 24\}$; verified even-resolvent cancellation identity $b_{01}^2 = \|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1}$ to $10^{-18}$ at 50 dps; verified transmission ratio $\mathcal{R}_1 \in [3.91, 8.00]$ is $\mathcal{O}(1)$) |
| **M15** | Analytical resolution of transition dipole cancellation and even-resolvent Stieltjes derivative architecture | Paper 4B Proposition 8.20 | **COMPLETED** (Proved $Ku_1 = a_1 (Q_{\mathrm{even}} - \mu_1 I)^{-1} d$, $\|Ku_1\|^2 = a_1^2 G_d'(\mu_1)$, and established that $\|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1} = b_{01}^2$ is the exact Pythagorean isolation of the $k=0$ ground-state pole; proved $b_{01} = -\langle P_{\perp u_0} Kc, u_1 \rangle$ connects dipole suppression to the $99.99987\%$ ground-state alignment of $Kc$) |

---

*Document approved for implementation as the canonical guiding architecture for the Connes–CvS research programme.*

