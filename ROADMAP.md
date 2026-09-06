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

## 3. Analysis of Methodological Routes

We identify three distinct routes of increasing mathematical depth to establish operator dominance:

```
[Route A: Archimedean Positivity via Uncertainty] (Deprecated / Unlikely)
  Attempt to show high-frequency mass from Paley-Wiener forces Q_arch >= 0 independently.
        |  (Rejected: Explicit formula treats arch/prime/pole as coupled distribution)
        v
[Route B: Exact Weierstrass Decomposition & Resolvent Pairing] (Primary Analytical Route)
  Decompose h+(r) into Weierstrass poles -> Q_arch = C||v||^2 + sum [||v||^2/(n+1) - J(q_n)].
  Pair negative J(q_n) resolvents against explicit pole/prime boundary terms.
        |
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
*Objective:* Pair the negative resolvent sequence $-\sum_{n=0}^\infty J(q_n)$ directly with the pole form $\mathcal{Q}_{\mathrm{pole}}$ and prime form $\mathcal{Q}_{\mathrm{prime}}$ to produce manifestly positive blocks.

### Route C: Integrated Bivariate Kernel & Dual Lattice Sampling (Grand Unified Theory)
*Concept:* Formulate the complete Weil quadratic form as a single bivariate Fredholm integral on $[0, L] \times [0, L]$:
$$\mathcal{Q}_{\mathrm{Weil}}(v) = \iint_{[0, L]^2} T_v(x) \, \mathcal{K}_{\mathrm{Weil}}(x, y) \, T_v(y) \, dx \, dy.$$
The operator structure is governed by the collision of **two distinct discrete lattices**:
1. **The Fourier Spectral Lattice:** $a_m = \frac{2\pi m}{L}$, where $K_{\mathrm{Fourier}}(a_m) = \frac{L}{2} v_m^2$ samples the wave orthogonally.
2. **The Arithmetic Prime Lattice:** $\log(p^k)$, where prime-power Dirac masses act as discrete shifts on the spatial wave.

In this framework, Weil positivity transforms into a discrete measure perturbation problem on an entire function of exponential type, establishing a rigorous connection to the spectral theory of de Branges and Beurling spaces.

---

## 4. The Two Distinct Research Programmes

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

Paper 4 and Paper 4B are dedicated to Programme 1 and the bridge to Programme 2, respectively. Finite-rank positivity $\mathcal{Q}_{c, N} \succeq 0$ is a self-contained, mathematically rigorous theorem that does not require prior resolution of Programme 2.

---

## 5. Low-Mode Arithmetic Localization for $c = 13$

For the primary benchmark cutoff $c = 13$, an arithmetic coincidence severely restricts the dimension of the dangerous subspace:
$$L = \log 13 \approx 2.56495, \qquad a_m = \frac{2\pi m}{L} \approx 2.4496 \cdot m.$$

Comparing the discrete Fourier frequencies $a_m$ against the zero of $h_+(r)$ ($r_* \approx 6.28984$):
- Mode $m = 0$: $r = 0 < r_*$ ($h_+(0) \approx -5.372$, inside negative region)
- Mode $m = 1$: $a_1 \approx 2.4496 < r_*$ ($h_+(a_1) \approx -2.571$, inside negative region)
- Mode $m = 2$: $a_2 \approx 4.8992 < r_*$ ($h_+(a_2) \approx -0.738$, inside negative region)
- Mode $m = 3$: $a_3 \approx 7.3488 > r_*$ ($h_+(a_3) \approx +0.421$, **outside negative region**)
- Modes $m \ge 4$: $a_m \gg r_*$ (strictly positive, $h_+(a_m) \sim \log m$)

### The 3-Mode Dangerous Subspace
Only the **first three modes** $(v_0, v_1, v_2)$ interact directly with the negative spectral support of $h_+(r)$! For all modes $m \ge 3$, the lattice sampling energy $v_m^2 h_+(a_m) > 0$ is strictly positive.

This induces a block decomposition of the Galerkin space:
$$\mathbb{R}^{N+1} = \mathcal{V}_{\mathrm{low}} \oplus \mathcal{V}_{\mathrm{high}}, \qquad \mathcal{V}_{\mathrm{low}} = \operatorname{span}\{e_0, e_1, e_2\}, \quad \mathcal{V}_{\mathrm{high}} = \operatorname{span}\{e_3, \dots, e_N\},$$
decomposing the Weil matrix into:
$$\mathcal{Q}_{\mathrm{Weil}} = \begin{pmatrix} \mathcal{Q}_{\mathrm{low}} & \mathcal{Q}_{\mathrm{cross}} \\ \mathcal{Q}_{\mathrm{cross}}^T & \mathcal{Q}_{\mathrm{high}} \end{pmatrix}.$$
Because $\mathcal{Q}_{\mathrm{high}}$ is dominated by positive lattice weights $h_+(a_m) > 0$, Weil positivity reduces to controlling a **$3 \times 3$ low-frequency block** and its cross-coupling via Schur complements!

---

## 6. Concrete Experimental Suite: Cell 63a Specification

Before attempting formal proofs, we execute a targeted computational experiment to construct and inspect the negative operator directly:

### Objectives of `cell63a.py`
1. **Explicit Matrix Construction of $\mathcal{Q}_{\mathrm{arch}}^{(-)}$:**
   Compute the $(N+1) \times (N+1)$ positive semi-definite matrix:
   $$\big[\mathcal{Q}_{\mathrm{arch}}^{(-)}\big]_{mn} = \frac{1}{\pi} \int_0^{r_*} |h_+(r)| \Phi_m(r) \Phi_n(r) \, dr,$$
   in the canonical basis $v$, where $\Phi_m(r)$ are the canonical Fourier basis amplitudes.
2. **Explicit Construction of $\mathcal{Q}_{\mathrm{positive\ side}}$:**
   Compute $\mathcal{Q}_{\mathrm{positive}} = \mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}}^{(+)}$, where:
   $$\mathcal{Q}_{\mathrm{arch}}^{(+)} = \mathcal{Q}_{\mathrm{arch}} + \mathcal{Q}_{\mathrm{arch}}^{(-)}.$$
3. **Generalized Spectral Decomposition:**
   Solve the generalized eigenvalue problem:
   $$\mathcal{Q}_{\mathrm{positive}} x = \lambda \mathcal{Q}_{\mathrm{arch}}^{(-)} x.$$
   - **Positivity Certificate:** If $\lambda_{\min} \ge 1$, then $\mathcal{Q}_{\mathrm{Weil}} \succeq 0$ is rigorously certified on $\mathbb{R}^{N+1}$.
   - **Dangerous State Identification:** The generalized eigenvector $x_{\min}$ corresponding to $\lambda_{\min}$ defines the exact profile of the "most dangerous test vector" challenging Weil positivity.
4. **Dimension Sweep ($N \in \{4, 8, 12, 16, 20, 24\}$):**
   Track the trajectory of $\lambda_{\min}(N)$ to establish whether the dominance margin $\lambda_{\min} - 1$ is bounded away from zero or approaches 1.

---

## 7. The 7-Stage Strategic Project Roadmap

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
STAGE V: 3-MODE REDUCTION AND SCHUR COMPLEMENT BOUNDS
Status: PLANNED
Tasks:
  1. Project onto the 3-mode dangerous subspace V_low = span{v_0, v_1, v_2}.
  2. Prove positivity of the 3x3 low-frequency Schur complement.
  3. Bound the high-frequency tail via the discrete lattice inequality h_+(a_m) > 0.
========================================================================================
                                     |
                                     v
========================================================================================
STAGE VI: SEMICLASSICAL CONTINUUM LIMIT & SUBSPACES DENSITY
Status: PLANNED (Paper 4B)
Tasks:
  1. Prove polynomial bound |D_1/D_0| <= C N^p via sector-decomposed resolvents.
  2. Establish exponential boundary-defect decoupling D(N) -> 0.
  3. Prove strong resolvent convergence and density of union_{c, N} H_{c, N}.
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

## 8. Immediate Operational Milestones

| Milestone | Action Item | Target Artifact / Script | Deliverable |
| :---: | :--- | :--- | :--- |
| **M1** | Implement the finite-band negative operator and generalized eigenvalue suite | `cell63a.py` | `cell63a.out` (Certificate $\lambda_{\min}$) |
| **M2** | Audit the generalized spectrum across $N \in \{8, 12, 16, 20, 24\}$ | Analytical Review | Determination of the dominance margin |
| **M3** | Analyze the coordinates and spatial wave profile of the dangerous vector $x_{\min}$ | Diagnostic Report | Identification of the physical obstruction |
| **M4** | Algebraically pair $J(q_n)$ with the pole and prime representations | Paper 4B Section Update | Exact positive block formulation |

---

*Document approved for implementation as the canonical guiding architecture for the Connes–CvS research programme.*
