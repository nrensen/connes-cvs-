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

The lattice sampling suggests a natural three-mode low-frequency sector. The extent to which the negative Archimedean operator is effectively confined to, or dominated by, this sector is an open question to be tested directly by Cell 63a. The generalized eigenvectors of $(\mathcal{Q}_{\mathrm{positive}}, \mathcal{Q}_{\mathrm{arch}}^{(-)})$ will reveal whether the dangerous subspace actually collapses onto $\operatorname{span}\{v_0, v_1, v_2\}$ or requires broader modal participation.

---

## 7. Concrete Experimental Suite: Cell 63a Specification

Before attempting formal proofs, we execute a targeted computational reconnaissance to construct and inspect the negative operator directly:

### Objectives & Mathematical Protocol of `cell63a.py`
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
   - **Numerical Dominance Test:** If the computed generalized spectrum satisfies $\lambda_{\min} > 1$ with a margin substantially larger than numerical error, this provides strong finite-$N$ numerical evidence for $\mathcal{Q}_{\mathrm{Weil}} \succeq 0$. A rigorous certificate would require certified eigenvalue / interval error bounds (planned as a subsequent Cell 63b).
   - **Dangerous State Identification:** The generalized eigenvector $x_{\min}$ corresponding to $\lambda_{\min}$ defines the exact profile of the "most dangerous test vector" challenging Weil positivity.
5. **Three Core Diagnostic Suites:**
   - **Diagnostic 1 (Spectrum & Effective Rank of $\mathcal{Q}_-$):** Tabulate the full spectrum of $\mathcal{Q}_{\mathrm{arch}}^{(-)}$ to determine its conditioning and effective dimensional rank across $N \in \{4, 8, 12, 16, 20, 24\}$.
   - **Diagnostic 2 (Modal Projection of Dangerous State):** Compute the fractional energy projection of $x_{\min}$ onto the 3-mode sector $\operatorname{span}\{e_0, e_1, e_2\}$ vs. its orthogonal complement $\operatorname{span}\{e_3, \dots, e_N\}$ to test the Three-Mode Sector Hypothesis directly.
   - **Diagnostic 3 (Schur-Complement Test for High Modes):** Partition the full Weil matrix as $\mathcal{Q}_{\mathrm{Weil}} = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$, where $A$ is $3 \times 3$ (low modes) and $C$ is $(N-2) \times (N-2)$ (high modes). Verify whether $C \succ 0$ and evaluate the spectrum of the Schur complement:
     $$S_{\mathrm{low}} = A - B C^{-1} B^T$$
     to determine whether high modes can be eliminated while rigorously preserving positivity.


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
STAGE III: LOW-FREQUENCY OPERATOR DECOMPOSITION & CERTIFICATION
Status: ACTIVE TARGET (Cell 63a)
Tasks:
  1. Construct explicit matrix for Q_arch^{(-)} on [0, r_*].
  2. Compute generalized spectrum of (Q_positive, Q_arch^{(-)}).
  3. Verify lambda_min >= 1 across dimensions N = 4 ... 24.
  4. Characterize the extremal "most dangerous" test vector x_min.
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
STAGE V: THREE-MODE SECTOR TESTING AND OPERATOR COUPLING
Status: PLANNED
Tasks:
  1. Inspect the principal components of the extremal eigenvector x_min from Cell 63a.
  2. Test the dominance of the 3-mode sector V_low = span{v_0, v_1, v_2} vs off-lattice coupling.
  3. Bound the high-frequency tail via the discrete lattice inequality h_+(a_m) > 0.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE VI: CONTINUUM ASYMPTOTICS & SUBSPACES DENSITY
Status: PLANNED (Paper 4B)
Tasks:
  1. Prove polynomial bound |D_1/D_0| <= C N^p via sector-decomposed resolvents.
  2. Establish exponential boundary-defect decoupling D(N) -> 0.
  3. Formulate formal WKB double-well potential and prove density of union_{c, N} H_{c, N}.
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

## 9. Immediate Operational Milestones

| Milestone | Action Item | Target Artifact / Script | Deliverable |
| :---: | :--- | :--- | :--- |
| **M1** | Implement the finite-band negative operator and generalized eigenvalue suite | `cell63a.py` | `cell63a.out` (Certificate $\lambda_{\min}$) |
| **M2** | Audit the generalized spectrum across $N \in \{8, 12, 16, 20, 24\}$ | Analytical Review | Determination of the dominance margin |
| **M3** | Analyze the coordinates and spatial wave profile of the dangerous vector $x_{\min}$ | Diagnostic Report | Identification of the physical obstruction |
| **M4** | Algebraically pair $J(q_n)$ with the pole and prime representations | Paper 4B Section Update | Exact positive block formulation |

---

*Document approved for implementation as the canonical guiding architecture for the Connes–CvS research programme.*

