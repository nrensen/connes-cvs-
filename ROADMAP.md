# Project Roadmap: Finite-Rank Connes–van Suijlekom Galerkin Positivity and the Weil Explicit Formula

**Document Version:** 2.0 (Refocused Strategic Architecture)  
**Date:** September 2026  
**Status:** Canonical Project Strategy & Master Mathematical Roadmap  
**Repository Architecture:** [AGENTS.md](file:///c:/data/github/connes-cvs-/AGENTS.md) (Operating Standards) | [cell_history_map.md](file:///c:/data/github/connes-cvs-/cell_history_map.md) (Comprehensive Historical Research Notebook)  
**Associated Manuscripts:**
- Paper 4: *An Exact Resolvent and Commutator Toolkit* (Locked Toolkit Baseline)
- Paper 4B: *The Dirichlet Continuum Limit, Barrier Mechanics, and Asymptotic Weil Positivity* (Asymptotic Programme)
- Paper 5: *Continuum Operator Limits and the Weil Criterion* (In Preparation)

---

## 1. Executive Summary & Current Position

The Connes–van Suijlekom (2025) and Connes–Consani–Moscovici (2026) framework projects André Weil's explicit quadratic functional of prime number theory onto finite-rank Galerkin subspaces of frequency band $N$ on a logarithmic scaling interval $[0, L] = [0, \log c]$:

$$\mathcal{Q}_{\mathrm{Weil}}(v) = \mathcal{Q}_{\mathrm{pole}}(v) + \mathcal{Q}_{\mathrm{prime}}(v) + \mathcal{Q}_{\mathrm{arch}}(v).$$

In Paper 4, we proved the unconditional global non-negativity of the Fourier-side kernel:
$$K_{\mathrm{Fourier}}(v, r, L) = \Phi_v(r)^2 \ge 0 \qquad \forall r \in \mathbb{R}, \; \forall v \in \mathbb{R}^{N+1}.$$

However, **kernel non-negativity alone does not prove positivity of the Archimedean form or the Weil form**, because the Archimedean multiplier:
$$h_+(r) = \operatorname{Re}\psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi$$
is **strictly negative** on a bounded low-frequency interval $[0, r_*]$ with $r_* \approx 6.28984$ ($h_+(0) \approx -5.37218$). The Archimedean form $\mathcal{Q}_{\mathrm{arch}}$ is therefore intrinsically a signed quadratic form.

> ### Current Position Statement
> The near-term programme is defined by a single mathematical target: proving the sufficient joint-limit estimate for asymptotic tail extinction,
> $$\lim_{L \to \infty} \limsup_{N \to \infty} \Delta_j(N) R_{\mathrm{spec}}(N, L) = 0 \quad \Longrightarrow \quad \lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1.$$
> Gate 1 does not require proving a particular spectral growth law; individual growth models ($p < 2$, WKB exponential barrier penetration, polynomial bounds, coordinate transmission cancellation) are strictly **candidate analytical routes toward this core proposition**, not standalone objectives.
> In parallel, a controlled **early Gate 3 structural reconnaissance** is initiated to audit the algebraic fidelity of $Q_{c, N} \stackrel{?}{\longleftrightarrow} \mathcal{W}$ and the reconciliation between the Fourier lattice $\Lambda_{\mathrm{Fourier}} = \{2\pi m / L\}$ and prime arithmetic lattice $\Lambda_{\mathrm{arith}} = \{\log p^k\}$.

---

## 2. The Strategic Structure: Near-Term Dual Tracks & The Unbroken Chain

To maintain absolute mathematical integrity and prevent confusing computational progress with analytical proof, the near-term programme executes two tightly coordinated parallel tracks leading into the continuum chain of obligations:

```
+----------------------------------------------------------------------------------------------------+
| NEAR-TERM DUAL TRACKS (PARALLEL EXECUTION)                                                         |
|                                                                                                    |
| Track A (Gate 1 Core): Asymptotic Tail Enclosure                                                  |
|            lim   limsup  Delta_j(N) R_spec(N, L) = 0   ===>   Pi_{j, tail}(N, L) ---> 1            |
|           L->oo   N->oo                                                                            |
|                                                                                                    |
| Track B (Early Gate 3): Arithmetic Structural Reconnaissance                                       |
|            Q_{c, N}  <--?-->  W[g]    |    Lambda_Fourier {2pi m/L}  <--?-->  Lambda_arith {log p^k}|
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  | (Passage to Continuum Operator & Full Gate 3)
                                                  v
+----------------------------------------------------------------------------------------------------+
| LONG-TERM MATHEMATICAL HORIZON: CONTINUUM WEIL POSITIVITY & RH                                     |
|                                                                                                    |
|    Q_{c, N}  ======>  Q_{c, oo}  ======>  W[g] >= 0 on Idele Class Group  ======>  RH (Weil 1952)  |
+----------------------------------------------------------------------------------------------------+
```

### The Unbroken Chain of Obligations
Every analytical and computational effort in this repository must advance a specific link in the four-stage chain of obligations:

$$\boxed{\textbf{Finite Galerkin Matrix } \mathcal{Q}_{c, N} \;\xrightarrow{\textbf{Gates 1 \& 2}}\; \textbf{Continuum Operator } \mathcal{Q}_\infty \;\xrightarrow{\textbf{Gate 3}}\; \textbf{Weil Functional } \mathcal{W} \;\xrightarrow{\textbf{Gates 4 \& 5}}\; \textbf{Riemann Hypothesis}}$$

1. **Finite Galerkin Matrix ($\mathcal{Q}_{c, N} \succ 0$):** High-precision verification via $LDL^T$ Schur complement decoupling (verified at 80 dps across tested dimensions in `cell64.py`) and exact finite-$N$ identities.
2. **Continuum Operator ($\mathcal{Q}_\infty$):** Rigorous definition of the limiting operator/form $Q_\infty = \lim_{N \to \infty} Q_{c, N}$, its self-adjoint domain, boundary condition specification, and preservation of low-energy bound states.
3. **The Weil Bridge ($\mathcal{Q}_\infty \stackrel{?}{=} \mathcal{W}$):** Proving that the continuum quadratic form coincides identically with André Weil's explicit quadratic functional on the idele class group. **This is a hard arithmetic falsification gate.** Early reconnaissance ($\text{Gate 1 work} \parallel \text{early Gate 3 structural audit}$) verifies algebraic viability before full continuum construction.
4. **The Weil Criterion ($\mathcal{W}[g] \ge 0 \implies \mathrm{RH}$):** Proving density of the Galerkin test class in Weil's admissible class and deducing the Riemann Hypothesis via Weil (1952).

---

## 3. Mathematical Foundations & Established Architecture

The repository's permanent mathematical foundation is established across Paper 4 and Paper 4B:

### 3.1 Operator Dominance & Schur Decoupling
Splitting the Archimedean multiplier into positive and negative parts $h_+(r) = h_+^+(r) - h_+^-(r)$ where $h_+^-$ is supported strictly on $[0, r_*]$ induces the operator splitting:
$$\mathcal{Q}_{\mathrm{Weil}} = \underbrace{\mathcal{Q}_{\mathrm{pole}} + \mathcal{Q}_{\mathrm{prime}} + \mathcal{Q}_{\mathrm{arch}}^{(+)}}_{\mathcal{Q}_{\mathrm{positive}}} - \mathcal{Q}_{\mathrm{arch}}^{(-)}.$$
Positivity is mathematically equivalent to the operator dominance inequality:
$$\mathcal{Q}_{\mathrm{positive}} \succeq \mathcal{Q}_{\mathrm{arch}}^{(-)} \qquad \Longleftrightarrow \qquad \lambda_{\min}(\mathcal{Q}_{\mathrm{positive}}, \mathcal{Q}_{\mathrm{arch}}^{(-)}) \ge 1.$$
Because generalized eigenvalue whitening suffers from Gram matrix ill-conditioning ($\sigma_{\min}(\mathcal{Q}_-) \to 0$), positivity is certified via symmetric Schur complement block decoupling:
$$\mathcal{Q}_{\mathrm{Weil}} = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}, \qquad C \succ 0, \quad S_{\mathrm{low}} = A - B C^{-1} B^T \succ 0.$$
For $c = 13$, Fourier lattice frequencies $a_m = 2\pi m / L$ lie inside $[0, r_*]$ only for $m \in \{0, 1, 2\}$. The $3 \times 3$ effective operator $S_{\mathrm{low}}$ captures the full ground-state scale ($\lambda_0 \sim 10^{-43}$ at $N=24$) within $5\%$, while the high-mode block $C$ is decoupled and strongly positive ($\lambda_{\min}(C) \sim 10^{-27} \gg \lambda_0$).

### 3.2 Coordinate Wavepacket Geometry & Dipole Cancellation
By the exact finite-$N$ commutator identity $[K, Q] = \boldsymbol\psi d^T - d \boldsymbol\psi^T$, the coordinate wavepacket $Kc$ satisfies:
$$D_0^2 M_2 = \|Kc\|^2 = \sum_{m=1}^N m^2 v_m^2, \qquad D_0^2 M_{2,\mathrm{exc}} = \|P_{\perp u_0} Kc\|^2.$$
At $N=24$, $Kc$ is aligned with the ground state $u_0$ to $> 99.99987\%$. The relative tunneling gap factors as:
$$R_{\mathrm{gap}}(N) \equiv \frac{D_0^2}{\mu_1 - \lambda} = \frac{b_{01}^2}{\mathcal{R}_1}, \qquad \mathcal{R}_1 \in [3.91, 8.00] = \Theta(1),$$
where $b_{01}^2 = \|Ku_1\|^2 - \mathcal{E}_{\mathrm{even}, 1}$ is the exact Pythagorean projection isolating the ground-state pole $E_0 = \lambda$.

### 3.3 Regularized Stieltjes Function & Two-Pole Clustering
Defining the regularized Stieltjes derivative $H(\mu) \equiv (\mu - \lambda)^2 G_d'(\mu) = D_0^2 + \sum_{k=1}^N d_k^2 (\frac{\mu - \lambda}{E_k - \mu})^2$, spectral tail filtering is governed by:
$$\mathcal{Q}_j = \frac{H(\mu_1)}{H(\mu_j)}, \qquad T_j = \mathcal{K}_j \mathcal{Q}_j.$$
By strict parity interlacing $E_j < \mu_j < E_{j+1}$, the two adjacent poles $\{E_j, E_{j+1}\}$ carry $> 99.997\%$ of the total spectral mass of $H(\mu_j)$. The pole asymmetry identity:
$$\frac{H_{j+1}(\mu_j)}{H_j(\mu_j)} = \alpha_j \left(\frac{L_j}{R_j}\right)^2, \qquad \alpha_j \equiv \frac{d_{j+1}^2}{d_j^2}, \quad L_j = \mu_j - E_j, \quad R_j = E_{j+1} - \mu_j,$$
reveals that boundary amplification $\alpha_2 \approx 9.21 \times 10^4$ perfectly balances geometric gap asymmetry $(R_2/L_2)^2 \approx 6.03 \times 10^4$, yielding the $\mathcal{O}(1)$ pole balance $H_3/H_2 \approx 1.528 = \Theta(1)$.

### 3.4 Stieltjes Product Representation & Universal Interlacing Tail Bound
The consecutive boundary weight ratio factors as $\alpha_j = \zeta_j \Pi_j$, where $\zeta_j = \frac{E_{j+1} - z_j^*}{z_j^* - E_j}$ and $\Pi_j = \prod_{\ell \ne j} \omega_{j, \ell}$. The pairwise factor satisfies the exact identity:
$$\omega_{j, \ell} - 1 = \frac{\Delta_j \delta_\ell}{|E_j - z_\ell^*| |E_{j+1} - E_\ell|}, \qquad \delta_\ell = z_\ell^* - E_\ell.$$
From Stieltjes interlacing $E_\ell < z_\ell^* < E_{\ell+1}$, the displacement is unconditionally bounded by $\delta_\ell < \Delta_\ell \equiv E_{\ell+1} - E_\ell$, yielding the **Universal Interlacing Tail Bound (Lemma 8.27)**:
$$0 < \omega_{j, \ell} - 1 < \eta_{\mathrm{inter}}(j, \ell) \equiv \frac{\Delta_j \Delta_\ell}{(E_\ell - E_j)(E_\ell - E_{j+1})}, \qquad \forall \ell \ge j+2.$$
This completely eliminates boundary weights $d_\ell^2$, ladder ratios $\alpha_\ell$, and sign ratios $\varepsilon_\ell$ from the asymptotic tail.

### 3.5 Operator-Norm Enclosure & The Dual-Regime Hierarchy
Summing the interlacing bound for $\ell > L$ yields the telescoping operator-norm bound:
$$\sum_{\ell > L} (\omega_{j, \ell} - 1) < \mathcal{S}_{\mathrm{inter}}(N; L) \le \mathcal{E}_j^{\mathrm{exact}}(N; L) \le \mathcal{E}_j^{\mathrm{op}}(N; L) \equiv \Delta_j \frac{\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{(E_{L+1} - E_j)(E_{L+1} - E_{j+1})}.$$
Auditing in `cell96.out` and `cell97.out` established the exact hierarchy:
$$\boxed{S_{\mathrm{dev}} < S_{\mathrm{inter}} \le \mathcal{E}_2^{\mathrm{exact}} < \mathcal{E}_2^{\mathrm{op}} \quad \text{without exception across all } N, L.}$$
The slack factor $\frac{E_N}{E_N - E_{L+1}}$ is $\le 1.0003$ for $L \in \{4, 8\}$ and $1.52$ for $L=11, N=24$. Tail extinction operates across two distinct mathematical regimes:
- **Regime A (Fixed Archimedean Cutoff $T$):** By Loewner divided-difference smoothness, $\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T) < \infty$ uniformly in $N$. The spectrum does not escape to infinity. Tail extinction follows unconditionally once the continuum base stabilizes: $\inf_N (E_{L+1}^{(N)} - E_{j+1}^{(N)}) > 0$.
- **Regime B (Resolution-Preserving Scaling $T(N) > \alpha_N = \frac{2\pi N}{L}$):** Even if the upper spectrum $E_N$ grows, the bound-state tunneling splitting $\Delta_j(N) \to 0$ collapses exponentially ($\Delta_2 \approx 1.37 \times 10^{-26}$ at $N=24$), driving $\mathcal{E}_2^{\mathrm{op}} \to 3.06 \times 10^{-26}$. Tail extinction is governed by the product $\Delta_j(N) R_{\mathrm{spec}}(N, L) \to 0$.

---

## 4. The Five Conceptual Mathematical Gates

Rather than tracking sequential computational scripts, the research programme is organized around five strictly sequenced **Mathematical Gates**. Progress is measured exclusively by resolving the analytical obstructions defining each gate.

```
+----------------------------------------------------------------------------------------------------+
| GATE 1: FINITE-N SPECTRAL MECHANISM & ASYMPTOTIC TAIL EXTINCTION (ACTIVE GATE)                     |
| Obstruction: Rigorous bound on tunneling splitting Delta_j(N) and relative growth R_spec(N, L).   |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
| GATE 2: CONTINUUM LIMITING OPERATOR & QUADRATIC FORM                                               |
| Obstruction: Self-adjoint domain, Dirichlet boundary layer, resolvent and form convergence.        |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
| GATE 3: ARITHMETIC RECONNECTION & THE WEIL BRIDGE (HARD FALSIFICATION GATE)                        |
| Obstruction: Proving Q_oo[f] = W[f] on the idele class group via dual lattice sampling.          |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
| GATE 4: POSITIVITY EQUIVALENCE                                                                     |
| Obstruction: Equivalence between Q_oo >= 0 and Weil positivity W[g] >= 0.                          |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
| GATE 5: INVOCATION OF THE WEIL CRITERION & THE RIEMANN HYPOTHESIS                                  |
| Obstruction: Density of Galerkin subspaces in Weil's test class and deduction of RH.               |
+----------------------------------------------------------------------------------------------------+
```

---

### Gate 1: Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction (ACTIVE GATE)

- **Central Mathematical Proposition (Gate 1 Target):**
  Prove that the remote Stieltjes product converges to unity in the joint limit:
  $$\boxed{\lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{E}_j^{\mathrm{exact}}(N, L) = 0 \quad \Longrightarrow \quad \lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1.}$$
  Via the operator-norm enclosure (Paper 4B Lemma 8.28), a sufficient condition for this proposition is:
  $$\boxed{\lim_{L \to \infty} \limsup_{N \to \infty} \Delta_j(N) R_{\mathrm{spec}}(N, L) = 0.}$$
  *Gate 1 does not require proving a particular spectral growth law or enforcing a specific exponent. All asymptotic models ($p < 2$, WKB tunneling, polynomial gap bounds, coordinate wavepacket transmission cancellation) are candidate routes toward this proposition, not standalone objectives.*

- **Central Obstruction:** Transitioning from discrete finite-$N$ numerical certification (where $\mathcal{E}_2^{\mathrm{op}} \le 3.06 \times 10^{-26}$ at $N=24$) to an unconditional analytical proof that $\Delta_j(N) R_{\mathrm{spec}}(N, L) \to 0$.

- **Competing Analytical Routes & Component Estimates:**
  1. **Route 1A (Exponential Tunneling Splitting with Soft Growth — Dominant Route):**
     - *Mechanics:* Establish an analytical bound on bound-state tunneling splitting $\Delta_j(N) \le C_j e^{-\sigma_j N}$ with $\sigma_j > 0$ derived from the discrete Galerkin potential barrier, coupled with any soft/polynomial bound on the upper spectrum $R_{\mathrm{spec}}(N, L) = o(e^{-\sigma_j N})^{-1} = o(e^{\sigma_j N})$.
     - *Status:* Empirically confirmed ($\Delta_2 \sim 10^{-26}$ at $N=24$), analytical WKB barrier modeling active in Paper 4B Section 8.27.
  2. **Route 1B (Fixed Cutoff $T$ Uniform Loewner Boundedness):**
     - *Mechanics:* Under a fixed Archimedean cutoff $T$, Loewner divided-difference smoothness ensures $\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T) < \infty$ uniformly in $N$. The spectrum cannot escape to infinity ($R_{\mathrm{spec}}$ bounded). Tail extinction follows unconditionally once the continuum Ritz base stabilizes: $\inf_N (E_{L+1}^{(N)} - E_{j+1}^{(N)}) \ge g_* > 0$.
     - *Status:* High-$T$ recovery sweeps (`cell95.py`) ongoing to certify $g_{11} \approx 0.418$.
  3. **Route 1C (Coordinate Wavepacket Transmission Cancellation):**
     - *Mechanics:* Bounding the excitation residual $\|P_{\perp u_0} Kc\|^2$ and demonstrating that mode-by-mode destructive phase interference quenches oscillatory transmission across the barrier top.
     - *Status:* Exact Pythagorean decomposition established in Paper 4; reduces the problem to dipole excitation residual and mode transmission factors.
  4. **Route 1D (Direct Operator-Theoretic / Resolvent Bypass):**
     - *Mechanics:* Bypassing mode-by-mode product estimates entirely via direct trace-class resolvent convergence or relative compactness of the perturbation $(Q_{\mathrm{even}}^{(N)} - \mu)^{-1} - (Q_{\mathrm{even}}^{(\infty)} - \mu)^{-1}$.
     - *Status:* Candidate analytical formulation under investigation for Paper 5.

---

### Gate 2: Continuum Limiting Operator & Quadratic Form

- **Mathematical Target:** Construct the limiting operator $Q_\infty = \lim_{N \to \infty} Q_{c, N}$ as a self-adjoint, semibounded operator on an appropriate Hilbert space $\mathcal{H}$, and establish strong resolvent and form convergence.
- **Central Obstruction:** The finite Galerkin matrices $Q_{c, N}$ are defined via discrete trigonometric projections with Dirichlet boundary conditions. In the continuum limit, the boundary layer at $t = 0, L$ must be reconciled with the domain of the singular integral operator.
- **Concrete Sub-Objectives:**
  1. **Objective 2.1 (Domain & Boundary Conditions):** Define the dense domain $\mathcal{D}(Q_\infty) \subset L^2([0, L])$ and identify the boundary conditions satisfied by continuum wavefunctions (Dirichlet vs Robin vs free endpoint).
  2. **Objective 2.2 (Loewner Monotonicity & Friedrichs Extension):** Exploit the established finite-$N$ Loewner monotonicity ($\Delta \Sigma(N) \succ 0$) to construct $Q_\infty$ via monotone quadratic form convergence or Friedrichs extension.
  3. **Objective 2.3 (Bound-State Preservation):** Prove that the discrete bound-state ladder beneath the barrier top ($E_0, \dots, E_{11}$) converges to genuine discrete eigenvalues of $Q_\infty$ with finite cardinality $\bar{N}_{\mathrm{bound}} < \infty$ (Bridges B1 and B2).

---

### Gate 3: Arithmetic Reconnection & The Weil Bridge (HARD FALSIFICATION GATE)

- **Mathematical Target:** Prove that the continuum quadratic form evaluates identically to André Weil's explicit quadratic functional:
  $$\langle f, Q_\infty f \rangle \stackrel{?}{=} \mathcal{W}[f] \qquad \forall f \in \mathcal{S}_{\mathrm{Weil}}.$$
- **Central Obstruction:** Bridging the two incommensurate sampling structures: the periodic Fourier lattice $\Lambda_{\mathrm{Fourier}} = \{2\pi m / L\}$ and the rigid arithmetic prime lattice $\Lambda_{\mathrm{arith}} = \{\log p^k\}$.

- **Early Gate 3 Structural Reconnaissance (Parallel Track):**
  > [!NOTE]
  > **Reconnaissance Principle ($\text{Gate 1 work} \parallel \text{early Gate 3 structural audit}$):**
  > Rather than deferring the arithmetic identification until the full continuum operator $Q_\infty$ is constructed, a small, tightly controlled structural investigation is conducted in parallel with Gate 1.
  > 
  > **The Central Arithmetic Questions:**
  > 1. *Exact Algebraic Decomposition:* What exact structural feature of the finite Galerkin matrix $Q_{c, N}$ is supposed to become the prime lattice contribution $\sum_{p^k \le c} \frac{\log p}{p^{k/2}} g(\log p^k)$, and what becomes the gamma/pole contribution?
  > 2. *Lattice Incommensurability:* The discrete Fourier lattice $\Lambda_{\mathrm{Fourier}} = \{2\pi m / L\}$ and arithmetic prime lattice $\Lambda_{\mathrm{arith}} = \{\log p^k\}$ represent fundamentally different geometries. Does the Galerkin trigonometric projection preserve arithmetic prime-power localization, or does it introduce an irreversible continuum aliasing error?
  >
  > **Early Falsification Criterion:** If the discrete finite-$N$ decomposition of $Q_{c, N}$ cannot be mapped algebraically to the Weil prime-plus-archimedean structure with vanishing residue as $N \to \infty$, the Galerkin truncation represents an isolated toy model disconnected from prime number theory.

- **Hard Falsification Protocol:**
  > [!CAUTION]
  > **Arithmetic Falsification Protocol:** If the continuum quadratic form $Q_\infty$ fails to match André Weil's explicit functional $\mathcal{W}$ on the idele class group, the finite-rank Galerkin construction represents an isolated toy model and cannot prove the Riemann Hypothesis.
  > If falsification occurs at Gate 3, the research programme must **immediately halt**, record the exact algebraic discrepancy, and pivot to reformulating the Galerkin discretization. Extending numerical sweeps past a failed arithmetic gate is strictly prohibited.

- **Concrete Long-Term Sub-Objectives:**
  1. **Objective 3.1 (DLMF Physical-Space Transfer):** Evaluate the Binet/Gauss integral representation of $h_+(r)$ pulled back into coordinate space against the prime Dirac comb $-\frac{1}{L}\sum_{p^k \le c} \frac{\log p}{p^{k/2}} \delta(y - \log p^k)$.
  2. **Objective 3.2 (Pole Absorption Proof):** Rigorously prove that the discrete Weierstrass resolvent sequence $-\sum_{n=0}^\infty J(q_n)$ pairs with and absorbs the zeta pole functional $\mathcal{Q}_{\mathrm{pole}}$ identically in the continuum limit.

---

### Gate 4: Positivity Equivalence

- **Mathematical Target:** Establish that operator non-negativity of the continuum operator is mathematically equivalent to the non-negativity of Weil's functional on test functions:
  $$Q_\infty \succeq 0 \qquad \Longleftrightarrow \qquad \mathcal{W}[g] \ge 0 \quad \forall g \in C_c^\infty(C_{\mathbb{Q}}).$$
- **Central Obstruction:** Resolving the Archimedean sign problem in the continuum: demonstrating that the negative well on $[0, r_*]$ is absorbed by operator dominance $\mathcal{Q}_{\mathrm{positive}} \succeq \mathcal{Q}_{\mathrm{arch}}^{(-)}$ without relying on ground-state asymptotics.
- **Concrete Sub-Objectives:**
  1. **Objective 4.1 (Infinite-Dimensional Schur Decoupling):** Extend the finite-$N$ Schur complement positivity $S_{\mathrm{low}} \succ 0$ and $C \succ 0$ to a rigorous closed-subspace factorization on $\mathcal{H}_{\mathrm{low}} \oplus \mathcal{H}_{\mathrm{high}}$.
  2. **Objective 4.2 (Spectral Measure Positivity):** Prove that the positive distribution $\mathcal{W}[g]$ inherits positive definiteness from $Q_\infty \succeq 0$ across the full test class.

---

### Gate 5: Invocation of the Weil Criterion & The Riemann Hypothesis

- **Mathematical Target:** Conclude the Riemann Hypothesis via André Weil's 1952 criterion:
  $$\mathcal{W}[g \ast g^*] \ge 0 \quad \forall g \in \mathcal{S}(C_{\mathbb{Q}}) \qquad \Longrightarrow \qquad \text{All non-trivial zeros of } \zeta(s) \text{ have } \operatorname{Re}(s) = \frac{1}{2}.$$
- **Central Obstruction:** Density of the Galerkin truncation spaces $\bigcup_{c > 1} \bigcup_{N \ge 1} \mathcal{H}_{c, N}$ in Weil's admissible class on the idele class group $C_{\mathbb{Q}} = \mathbb{A}_{\mathbb{Q}}^\times / \mathbb{Q}^\times$.
- **Concrete Sub-Objectives:**
  1. **Objective 5.1 (Two-Parameter Density Theorem):** Prove that the double family of finite-rank test functions is dense in Weil's test space in the inductive limit topology as $c \to \infty$ and $N \to \infty$.
  2. **Objective 5.2 (Unconditional Adelic Positivity):** Deduce $\mathcal{W}[g] \ge 0$ unconditionally for all admissible $g$.
  3. **Objective 5.3 (Final RH Conclusion):** Apply Weil's theorem to establish the Riemann Hypothesis.

---

## 5. Operational Principles for Computational Cells & Strategic Governance

To maintain repository standards and prevent `ROADMAP.md` from degenerating into a cell log:

1. **Subordination of Computation to Mathematical Gates:**
   Computational cells are subordinate diagnostic instruments designed to test, verify, or falsify specific mathematical conjectures within an active Gate. **A computational cell is never a roadmap milestone in itself.**
2. **Pre-Flight Formulation Requirement:**
   Before authoring any new computational script (`cell*.py`), the author must document:
   - **Target Gate:** Which of Gates 1–5 does this script serve?
   - **Target Proposition / Hypothesis:** What exact mathematical statement is being tested?
   - **Verification / Falsification Criterion:** What precise quantitative outcome confirms or refutes the conjecture?
   - **Bridge to Continuum:** How does the computed quantity advance the chain toward the Weil form?
3. **Repository Separation of Responsibilities:**
   - [cell_history_map.md](file:///c:/data/github/connes-cvs-/cell_history_map.md): Records every cell chronologically, detailing hypotheses, methods, raw outputs, established results, and refuted conjectures.
   - [ROADMAP.md](file:///c:/data/github/connes-cvs-/ROADMAP.md): Defines active mathematical gates, analytical obstructions, and the strategic chain of obligations.
   - [paper4.md](file:///c:/data/github/connes-cvs-/paper4_exact_resolvent_and_dirichlet_limit.md) & [paper4b.md](file:///c:/data/github/connes-cvs-/paper4b_dirichlet_continuum_limit_and_barrier_mechanics.md): Capture vetted, permanent mathematics with strict epistemic labeling.
4. **Roadmap Update Granularity & Strategic Lifecycle:**
   `ROADMAP.md` operates at the level of gates and propositions, not computational steps or routes. `ROADMAP.md` should **rarely change after a cell succeeds**. It changes when our **belief about the research strategy** changes:
   - *No Roadmap Update:* "Cell 98 confirms exponential tunneling to $N=32$" is a historical computational datum logged in `cell_history_map.md`.
   - *Justified Roadmap Update:* "An analytical proof of the tunneling bound is completed, resolving the high-spectrum growth obstruction and closing Gate 1" represents a change in strategic state that justifies a roadmap update.
   - Routes can change rapidly as analytical ideas evolve; gates and target propositions change only when our fundamental mathematical understanding changes.

---

## 6. Active Strategic Milestones (Gate-Oriented)

The following operational milestones define the active analytical and computational pipeline:

| Milestone | Target Gate | Mathematical Objective | Target Artifact / Script | Status |
| :---: | :---: | :--- | :--- | :---: |
| **M-G1.0** | Gate 1 | **Central Target Proposition 1.0:** Prove joint-limit tail extinction $\lim_{L \to \infty} \limsup_{N \to \infty} \Delta_j(N) R_{\mathrm{spec}}(N, L) = 0 \implies \Pi_{j, \mathrm{tail}} \to 1$ | Paper 4B Proposition 8.32 | **ACTIVE** |
| **M-G1.1** | Gate 1 | Route 1A Analysis: Analytical upper bound on tunneling splitting $\Delta_j(N) \le C_j e^{-\sigma_j N}$ from Galerkin barrier potential | Paper 4B Section 8.27 | **IN PROGRESS** |
| **M-G1.2** | Gate 1 | Route 1A Analysis: Minimal high-spectrum growth bound $R_{\mathrm{spec}}(N, L) = o(e^{\sigma_j N})$ under resolution scaling $T(N) > \alpha_N$ | Paper 4B Section 8.27 | **PLANNED** |
| **M-G1.3** | Gate 1 | Route 1B Analysis: High-$T$ Recovery Sweep across $T \in [200, 1600]$ at $N \in \{160, 176, 192\}$ certifying boundary gap $g_{11} \approx 0.418$ | `cell95.py` / `cell95.out` | **IN PROGRESS** |
| **M-G3.0** | Gate 3 | **Early Gate 3 Structural Reconnaissance:** Algebraic audit of finite-$N$ Galerkin matrix $Q_{c, N}$ vs prime lattice comb $\Lambda_{\mathrm{arith}}$ and pole absorption | Research Note / Paper 5 Section 1 | **ACTIVE (PARALLEL)** |
| **M-G2.1** | Gate 2 | Formulation of the continuum limiting quadratic form $Q_\infty$ and Friedrichs domain $\mathcal{D}(Q_\infty)$ | Paper 5 Section 2 | **PLANNED** |
| **M-G2.2** | Gate 2 | Proof of strong resolvent and form convergence $Q_{c, N} \to Q_\infty$ via Loewner monotonicity $\Delta \Sigma(N) \succ 0$ | Paper 5 Section 3 | **PLANNED** |
| **M-G2.3** | Gate 2 | Rigorous proof of Bridges B1 ($\bar{N}_{\mathrm{bound}} < \infty$) and B2 ($\sup_N \|Ku_j^{(N)}\|^2 < \infty$) for discrete Galerkin operators | Paper 4B Proposition 8.33 | **PLANNED** |
| **M-G3.1** | Gate 3 | Analytical evaluation of $Q_\infty$ against prime Dirac delta comb and proof of Weil functional identity $Q_\infty[f] = \mathcal{W}[f]$ | Paper 5 Section 4 | **FALSIFICATION GATE** |
| **M-G3.2** | Gate 3 | DLMF Binet integral coordinate-space transfer and dual lattice sampling ($\Lambda_{\mathrm{Fourier}}$ vs $\Lambda_{\mathrm{arith}}$) | Paper 5 Section 5 | **PLANNED** |
| **M-G4.1** | Gate 4 | Proof of continuum operator dominance $\mathcal{Q}_{\mathrm{positive}} \succeq \mathcal{Q}_{\mathrm{arch}}^{(-)}$ on $\mathcal{H}$ | Paper 5 Section 6 | **PLANNED** |
| **M-G5.1** | Gate 5 | Double density theorem ($\bigcup_{c, N} \mathcal{H}_{c, N}$ dense in $\mathcal{W}_{\mathrm{admissible}}$) and deduction of RH via Weil (1952) | Paper 5 Section 7 | **LONG-TERM TARGET** |

---

*Document approved as the canonical guiding architecture for the Connes–CvS research programme. All future computational cells and manuscript revisions must trace directly to an active Gate defined herein.*
