# Gate 1 Analytical Reduction: The Tail Extinction Theorem

**Status:** Working analytical note (Gate 1 Work Item B)
**Target:** [ROADMAP.md Gate 1 — M-G1.B](file:///c:/data/github/connes-cvs-/ROADMAP.md)
**Dependencies:** Paper NR2 Propositions 8.29, 8.30, 8.31 (all established)
**Date:** 2026-09-10

---

## Executive Summary

This note formalizes the analytical reduction of the Gate 1 tail extinction proposition:
$$\mathbf{H}_{\mathrm{tail}}(j): \quad \lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{E}_j^{\mathrm{exact}}(N, L) = 0 \quad \Longrightarrow \quad \lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1.$$

| Analytical Route | Required Hypotheses | Status of Each Hypothesis | Power |
| :--- | :--- | :--- | :--- |
| **Route 1B** (Fixed-$T$) | (I) $\mathrm{H}_{\mathrm{gap}}(J)$: uniform spectral gap; (II) $\mathrm{H}_{\mathrm{collapse}}(j)$: bound-state Ritz collapse | Both are active numerical conjectures with strong empirical support | Fixed-$T$ regime only |
| **Route 1A** (Tunneling) | Exponential splitting $\Delta_j \le C e^{-\sigma N}$ and $R_{\mathrm{spec}} = o(e^{\sigma N})$ | Empirical WKB only; no rigorous proof | Handles $T(N) \to \infty$ |
| **Route 1C** (Wavepacket) | Dipole cancellation bound on $\|P_{\perp u_0} Kc\|^2$ | Pythagorean identity proved (Paper NR1) | Complementary mechanism |
| **Route 1D** (Resolvent) | Direct trace-class resolvent convergence | Formulation stage | Bypasses mode-by-mode |

**Main Result:** Route 1B yields a *rigorous conditional theorem* (Proposition 8.37 in Paper NR2) reducing the infinite-dimensional tail extinction problem to two explicit, testable hypotheses about the finite-dimensional Galerkin spectrum.

---

## 1. The Tail Extinction Problem

### 1.1 Precise Statement

Let $Q_{\mathrm{even}}^{(N)}$ denote the even-parity Galerkin matrix of rank $N+1$ with ordered eigenvalues $0 \le E_0^{(N)} < E_1^{(N)} < \cdots < E_N^{(N)}$ and spectral gaps $\Delta_k^{(N)} \equiv E_{k+1}^{(N)} - E_k^{(N)}$.

Let $j \ge 0$ be a fixed focus mode (the bound-state index whose stability under the continuum limit is being analysed). The consecutive boundary weight ratio factors as:
$$\alpha_j = \zeta_j \cdot \Pi_{j, \mathrm{core}}(L) \cdot \Pi_{j, \mathrm{tail}}(L),$$
where the remote tail product is:
$$\Pi_{j, \mathrm{tail}}(L) \equiv \prod_{\ell = L+1}^{N-1} \omega_{j, \ell}, \qquad \omega_{j, \ell} \equiv \frac{(E_{j+1} - z_\ell^*)(z_{\ell+1}^* - E_j)}{(z_\ell^* - E_j)(E_{j+1} - z_{\ell+1}^*)},$$
with $z_\ell^*$ the Stieltjes interlacing zeros satisfying $E_\ell < z_\ell^* < E_{\ell+1}$.

**The Gate 1 Target Proposition ($\mathbf{H}_{\mathrm{tail}}(j)$):**
$$\boxed{\lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{E}_j^{\mathrm{exact}}(N, L) = 0 \quad \Longrightarrow \quad \lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1.}$$

The implication follows from the unconditional enclosure $1 < \Pi_{j, \mathrm{tail}} \le \exp(\mathcal{E}_j^{\mathrm{exact}})$ (Proposition 8.31).

### 1.2 The Existing Enclosure Chain

The exact telescoping exponent (Proposition 8.31, Paper NR2) provides the rigorous finite-$N$ bounding hierarchy:
$$\mathcal{S}_{\mathrm{dev}} < \mathcal{S}_{\mathrm{inter}} \le \mathcal{E}_j^{\mathrm{exact}}(N, L) \le \mathcal{E}_j^{\mathrm{op}}(N, L),$$
where:
$$\mathcal{E}_j^{\mathrm{exact}}(N, L) = \frac{\Delta_j^{(N)} \big(E_N^{(N)} - E_{L+1}^{(N)}\big)}{(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)})}, \qquad \mathcal{E}_j^{\mathrm{op}}(N, L) = \frac{\Delta_j^{(N)} \cdot \|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)})}.$$

The factored form is:
$$\mathcal{E}_j^{\mathrm{op}}(N, L) = \Delta_j^{(N)} \cdot R_{\mathrm{spec}}(N, L), \qquad R_{\mathrm{spec}}(N, L) \equiv \frac{\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)})}.$$

### 1.3 The Three Regimes

The product $\Delta_j R_{\mathrm{spec}} \to 0$ can be achieved via three distinct mechanisms:
- **Regime A:** $R_{\mathrm{spec}} \to 0$ independently (spectral growth slower than denominator).
- **Regime B (empirically dominant):** $\Delta_j \to 0$ while $R_{\mathrm{spec}} = \Theta(1)$.
- **Regime C:** Both factors contribute; $\Delta_j$ decays exponentially, $R_{\mathrm{spec}}$ grows subexponentially.

At $N = 24$, $T = 400$, $L = 11$: $\Delta_2 \approx 1.37 \times 10^{-26}$, $R_{\mathrm{spec}} \approx 2.23$ — firmly in Regime B.

---

## 2. Route 1B — Fixed-$T$ Analytical Reduction (Strongest Available)

### 2.1 The Analytical Ingredients

Route 1B assembles three established results from Paper NR2 into a conditional tail-extinction theorem:

**Ingredient 1 — Uniform Operator-Norm Bound (Proposition 8.29, Rigorous):**
At fixed cutoff parameters $c > 1$ and $T \ge 1$, the even Galerkin operator norm is uniformly bounded:
$$\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T) < \infty \qquad (\forall N \ge 1).$$
This confines all eigenvalues to the compact interval $[0, M(c, T)]$.

**Ingredient 2 — Monotone Ritz Convergence (Proposition 8.30, Rigorous):**
By the Courant–Fischer–Weyl min-max principle on nested Galerkin subspaces $V_N \subset V_{N+1}$:
$$E_k^{(N+1)} \le E_k^{(N)} \qquad (\forall N \ge 1, \; \forall k),$$
and the Ritz limits exist unconditionally:
$$E_k^{(\infty)} \equiv \lim_{N \to \infty} E_k^{(N)} = \inf_{N \ge k} E_k^{(N)} \ge 0 \qquad (\forall k \ge 0).$$

**Ingredient 3 — Telescoping Enclosure (Proposition 8.31, Rigorous):**
The remote tail is bounded by the exact telescoping exponent:
$$1 < \Pi_{j, \mathrm{tail}}(N; L) \le \exp\big(\mathcal{E}_j^{\mathrm{exact}}(N, L)\big) \le \exp\big(\mathcal{E}_j^{\mathrm{op}}(N, L)\big).$$

### 2.2 The Two Hypotheses

The Route 1B reduction requires two hypotheses about the Galerkin spectrum:

**Hypothesis $\mathrm{H}_{\mathrm{gap}}(J)$ (Uniform Boundary-Gap Separation):**
There exist an integer $J \ge j + 1$, a constant $\eta_J > 0$, and a finite threshold $N_0 \ge 1$ such that:
$$E_{J+1}^{(N)} - E_J^{(N)} \ge \eta_J \qquad (\forall N \ge N_0).$$

> **Empirical Status:** At $J = 11$, the spectral gap $g_{11}(N) = E_{12}^{(N)} - E_{11}^{(N)}$ satisfies $g_{11} \in [0.54, 0.57]$ across $N \in \{28, \ldots, 44\}$ (`cell90`–`cell91`). The gap shows remarkable numerical invariance while individual eigenvalues $E_{11}$ and $E_{12}$ continue drifting downward. Active high-$N$ exploration is ongoing in `cell95`.

**Hypothesis $\mathrm{H}_{\mathrm{collapse}}(j)$ (Bound-State Ritz Collapse):**
$$\Delta_j^{(\infty)} \equiv \lim_{N \to \infty} \big[E_{j+1}^{(N)} - E_j^{(N)}\big] = 0.$$
That is, the continuum Ritz limits of modes $j$ and $j+1$ coincide: $E_{j+1}^{(\infty)} = E_j^{(\infty)}$.

> **Empirical Status:** For the deep tunneling modes ($j \le 10$), the gaps collapse exponentially: $\Delta_2(N) \approx 1.37 \times 10^{-26}$ at $N = 24$ (Table 8.25.29 in Paper NR2), decreasing by 6–7 orders of magnitude per $\Delta N = 4$. The 11-mode tunneling cluster converges toward a degenerate ground-state energy $E_k^{(\infty)} \approx 0$ for $k \le 10$.

### 2.3 The Conditional Reduction Theorem

**Proposition 8.37 (Route 1B Conditional Reduction of $\mathbf{H}_{\mathrm{tail}}$; see Paper NR2).**
*Let $c > 1$ and $T \ge 1$ be fixed cutoff parameters. Let $j \ge 0$ be a fixed focus mode. Assume:*
1. *$\mathrm{H}_{\mathrm{gap}}(J)$ holds for some $J \ge j + 1$.*
2. *$\mathrm{H}_{\mathrm{collapse}}(j)$ holds: $\Delta_j^{(\infty)} = 0$.*

*Then the tail extinction condition holds:*
$$\boxed{\lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{E}_j^{\mathrm{op}}(N, L) = 0,}$$
*and consequently:*
$$\lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1.$$

**Proof.** Fix any $L \ge J$. For all $N > L$:

*Step 1 (Numerator bound).* By the operator-norm ceiling (Proposition 8.29):
$$E_N^{(N)} - E_{L+1}^{(N)} \le E_N^{(N)} \le \|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T).$$

*Step 2 (Denominator bound).* Since $j + 1 \le J$ and eigenvalues are ordered:
$$E_{j+1}^{(N)} \le E_J^{(N)}.$$
By $\mathrm{H}_{\mathrm{gap}}(J)$:
$$E_{J+1}^{(N)} - E_J^{(N)} \ge \eta_J \qquad (\forall N \ge N_0),$$
so $E_{J+1}^{(N)} \ge E_J^{(N)} + \eta_J \ge E_{j+1}^{(N)} + \eta_J$.

Since $L \ge J$ and eigenvalues are ordered, $E_{L+1}^{(N)} \ge E_{J+1}^{(N)}$. Therefore:
$$E_{L+1}^{(N)} - E_{j+1}^{(N)} \ge E_{J+1}^{(N)} - E_{j+1}^{(N)} \ge \eta_J.$$
Similarly:
$$E_{L+1}^{(N)} - E_j^{(N)} \ge E_{L+1}^{(N)} - E_{j+1}^{(N)} \ge \eta_J.$$

Combining:
$$(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)}) \ge \eta_J^2.$$

*Step 3 (Envelope bound).* From Steps 1–2 and the telescoping identity:
$$\mathcal{E}_j^{\mathrm{op}}(N, L) = \frac{\Delta_j^{(N)} \cdot \|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)})} \le \frac{M(c, T)}{\eta_J^2} \cdot \Delta_j^{(N)}.$$

*Step 4 (Limit).* By $\mathrm{H}_{\mathrm{collapse}}(j)$:
$$\limsup_{N \to \infty} \mathcal{E}_j^{\mathrm{op}}(N, L) \le \frac{M(c, T)}{\eta_J^2} \cdot \Delta_j^{(\infty)} = 0.$$

Since this bound is uniform in $L \ge J$:
$$\lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{E}_j^{\mathrm{op}}(N, L) = 0.$$

Because $\mathcal{E}_j^{\mathrm{exact}} \le \mathcal{E}_j^{\mathrm{op}}$ (Proposition 8.31), $\mathbf{H}_{\mathrm{tail}}(j)$ follows, and:
$$1 < \Pi_{j, \mathrm{tail}}(N, L) \le \exp\big(\mathcal{E}_j^{\mathrm{op}}(N, L)\big) \longrightarrow 1. \qquad \blacksquare$$

### 2.4 Structural Interpretation

The theorem cleanly separates the tail extinction into two independent spectral conditions targeting different parts of the spectrum:

| Hypothesis | Controls | Spectral Region | Physical Meaning |
| :--- | :--- | :--- | :--- |
| $\mathrm{H}_{\mathrm{gap}}(J)$ | Denominator (spectral separation) | Barrier boundary ($E_J, E_{J+1}$) | Macroscopic gap separates tunneling modes from continuum |
| $\mathrm{H}_{\mathrm{collapse}}(j)$ | Numerator (tunneling splitting) | Deep bound states ($E_j, E_{j+1}$) | Consecutive bound states converge to degenerate Ritz limit |

**Neither hypothesis alone suffices:**
- Without $\mathrm{H}_{\mathrm{gap}}$: the denominator could vanish (as $E_{L+1}^{(N)} \to E_{j+1}^{(N)}$), producing a $0/0$ indeterminate form.
- Without $\mathrm{H}_{\mathrm{collapse}}$: if $\Delta_j^{(\infty)} > 0$, then $\mathcal{E}_j^{\mathrm{op}} \ge \Delta_j^{(\infty)} / \eta_J^2 > 0$ for all $L$, and the tail does not close.

### 2.5 Epistemic Status

> **Rigorous conditional theorem, pending two hypotheses.**
> The logical chain Proposition 8.29 + Proposition 8.30 + Proposition 8.31 $\implies$ (Proposition 8.37 conditional on $\mathrm{H}_{\mathrm{gap}} + \mathrm{H}_{\mathrm{collapse}}$) is a complete proof. No additional analytical tools are imported; the entire argument operates within the established Paper NR2 framework.
>
> Neither $\mathrm{H}_{\mathrm{gap}}(J)$ nor $\mathrm{H}_{\mathrm{collapse}}(j)$ is currently proved. Both have strong numerical support across all tested dimensions.

---

## 3. Route 1A — Conditional Tunneling Reduction

### 3.1 The Exponential Tunneling Hypothesis

**Hypothesis $\mathrm{H}_{\mathrm{tunnel}}(j)$:**
There exist constants $C_j < \infty$ and $\sigma_j > 0$ such that:
$$\Delta_j(N) \le C_j \, e^{-\sigma_j N} \qquad (\forall N \ge 1).$$

> **Empirical Status:** The WKB barrier computation (Paper NR2 Section 4.1) yields the semiclassical scaling:
> $$\mathcal{S}_{\mathrm{WKB}}(N, c) \approx \frac{\pi N}{4} \log c,$$
> predicting $\sigma_j \approx \frac{\pi}{4} \log c$. At $c = 13$: $\sigma_j \approx \frac{\pi}{4} \log 13 \approx 2.014$.
> Numerically, the actual boundary suppression exceeds the WKB prediction by $\sim 5\%$ (`cell47`). However, no rigorous analytical proof of exponential splitting exists for the discrete Galerkin operator.

### 3.2 Conditional Tail Extinction under Exponential Tunneling

**Proposition (Route 1A, Conditional).**
*If $\mathrm{H}_{\mathrm{tunnel}}(j)$ holds with exponent $\sigma_j > 0$, and if:*
$$R_{\mathrm{spec}}(N, L) = o\big(e^{\sigma_j N}\big) \qquad (N \to \infty, \text{ for each fixed } L),$$
*then $\mathbf{H}_{\mathrm{tail}}(j)$ holds.*

**Proof.** Immediate from the factored product:
$$\mathcal{E}_j^{\mathrm{op}}(N, L) = \Delta_j(N) \cdot R_{\mathrm{spec}}(N, L) \le C_j \, e^{-\sigma_j N} \cdot R_{\mathrm{spec}}(N, L) \to 0. \qquad \blacksquare$$

### 3.3 Comparison with Route 1B

Route 1A is strictly more powerful than Route 1B in one important respect: it handles the **resolution-preserving regime** $T(N) > \alpha_N = \frac{2\pi N}{L}$ where the operator norm $\|Q_{\mathrm{even}}^{(N)}\|$ may grow with $N$. As long as the growth is subexponential in $N$, the exponential tunneling damping defeats it.

However, Route 1A requires proving two much harder analytical statements:
1. **The exponential splitting bound $\mathrm{H}_{\mathrm{tunnel}}$** — this requires a rigorous WKB-type argument for a discrete Galerkin matrix, which is a non-trivial problem in semiclassical analysis. The effective potential $V_{\mathrm{conf}}(t) = T''(t)/T(t) + E$ is defined through the ground-state wavefunction itself, creating a circular dependence that must be resolved.
2. **The subexponential growth of $R_{\mathrm{spec}}$** — under fixed $T$, this is automatic ($R_{\mathrm{spec}}$ is bounded); under scaling $T(N)$, it requires understanding how the upper spectrum $E_N^{(N)}$ grows with the Archimedean bandwidth.

### 3.4 The WKB Analytical Gap

Converting the empirical WKB observation into a proof requires:
1. An analytical characterization of the effective Schrödinger potential as a function of the Galerkin matrix entries.
2. A rigorous discrete semiclassical tunneling estimate — not a continuous ODE Sturm–Liouville theorem applied naively to a discrete operator (per AGENTS.md: "Do not import continuous ODE theorems to discrete operators").
3. Uniform estimates across all $N$, not just pointwise verification at finitely many dimensions.

This is an open hard analytical problem (M-G1.1 in [ROADMAP.md](file:///c:/data/github/connes-cvs-/ROADMAP.md)).

---

## 4. Route Comparison Matrix

| | **Route 1B** (Fixed-$T$ Loewner) | **Route 1A** (Exponential Tunneling) | **Route 1C** (Wavepacket Cancellation) | **Route 1D** (Resolvent Bypass) |
| :--- | :--- | :--- | :--- | :--- |
| **Key Hypothesis** | $\mathrm{H}_{\mathrm{gap}}(J) + \mathrm{H}_{\mathrm{collapse}}(j)$ | $\mathrm{H}_{\mathrm{tunnel}}(j)$ + soft growth | Dipole residual bound | Trace-class perturbation |
| **Paper NR2 Tools Used** | Prop 8.29, 8.30, 8.31 | Prop 8.31, Section 4.1 | Commutator identity (Paper NR1) | — |
| **Regime** | Fixed $T$ only | All regimes (incl. $T(N) \to \infty$) | Mode-by-mode | Global |
| **Analytical Difficulty** | Low (conditional theorem is elementary) | High (WKB for discrete operators) | Medium (coordinate bounds) | High (functional analysis) |
| **Proof Status** | ✅ Conditional theorem proved | ❌ Conditional + unproved hypothesis | ⚠️ Pythagorean identity proved, bounds open | ❌ Formulation stage |
| **Closes Gate 1?** | Yes, if hypotheses proved | Yes, if hypotheses proved | Complementary information | Yes, if proved |
| **Connection to Gate 2** | $\mathrm{H}_{\mathrm{gap}}$ related to continuum spectral structure | Independent of Gate 2 | Independent | Subsumes Gate 2 |

---

## 5. Road to Closing the Hypotheses

### 5.1 Closing $\mathrm{H}_{\mathrm{gap}}(J)$

The uniform boundary-gap hypothesis requires proving that the spectral gap $g_J(N) = E_{J+1}^{(N)} - E_J^{(N)}$ remains bounded away from zero for all sufficiently large $N$.

**Possible analytical approaches:**
1. **Gate 2 spectral structure:** If the continuum operator $Q_\infty$ (constructed via Gate 2) has a spectral gap between its $J$-th and $(J+1)$-th eigenvalues, then by standard Rayleigh–Ritz convergence theory, the Galerkin gap $g_J(N) \to g_J(\infty) > 0$. This connects $\mathrm{H}_{\mathrm{gap}}$ to the continuum spectral problem.
2. **Direct monotonicity argument:** If $g_J(N)$ can be shown to be eventually non-decreasing (monotone from below), then $g_J(\infty) = \lim_{N \to \infty} g_J(N) \ge g_J(N_0) > 0$. However, the empirical data shows $g_{11}$ fluctuating within a narrow band ($0.54$–$0.57$), not strictly increasing.
3. **Operator-theoretic bound:** A lower bound on $g_J(N)$ could potentially be derived from the divided-difference structure of the Galerkin matrix via Schur complement or spectral perturbation arguments.

**Active investigation:** `cell95.py` is currently exploring $g_{11}(N)$ at high dimensions $N \in \{128, 160, 176, 192\}$ with calibrated Archimedean cutoffs.

### 5.2 Closing $\mathrm{H}_{\mathrm{collapse}}(j)$

The bound-state Ritz collapse hypothesis requires proving $E_{j+1}^{(\infty)} = E_j^{(\infty)}$ for the deep tunneling modes.

**Possible analytical approaches:**
1. **Continuum spectral degeneracy:** If $Q_\infty$ has a degenerate ground-state energy $\lambda_0 = 0$ with multiplicity $\ge j+2$ (i.e., the entire tunneling cluster collapses to a single point), then $E_k^{(\infty)} = 0$ for $k \le j+1$ and $\Delta_j^{(\infty)} = 0$. This is the Gate 2 approach.
2. **Direct exponential splitting via Route 1A:** If $\mathrm{H}_{\mathrm{tunnel}}(j)$ is proved (exponential decay $\Delta_j \le C e^{-\sigma N}$), then $\Delta_j^{(\infty)} = 0$ follows trivially. Thus Route 1A *implies* $\mathrm{H}_{\mathrm{collapse}}$.
3. **Loewner monotone convergence:** If the Galerkin sequence $Q_{\mathrm{even}}^{(N)}$ satisfies a Loewner monotonicity property (the increments $Q_{\mathrm{even}}^{(N+1)} - Q_{\mathrm{even}}^{(N)}$ are positive semidefinite in an appropriate embedding), then the convergence of eigenvalues can be strengthened beyond Ritz monotonicity.

### 5.3 The Logical Hierarchy of Routes

The relationships between routes clarify the strategic landscape:

$$\mathrm{H}_{\mathrm{tunnel}}(j) \implies \mathrm{H}_{\mathrm{collapse}}(j),$$
since $\Delta_j(N) \le C e^{-\sigma N} \to 0$.

$$\mathrm{H}_{\mathrm{gap}}(J) + \mathrm{H}_{\mathrm{tunnel}}(j) \implies \mathrm{H}_{\mathrm{gap}}(J) + \mathrm{H}_{\mathrm{collapse}}(j) \implies \mathbf{H}_{\mathrm{tail}}(j),$$
so Route 1A (with gap) is strictly stronger than Route 1B's hypotheses.

The value of Route 1B is that $\mathrm{H}_{\mathrm{collapse}}(j)$ is a weaker hypothesis than $\mathrm{H}_{\mathrm{tunnel}}(j)$: it requires only $\Delta_j \to 0$ (at any rate), not exponential decay at a specific rate.

---

## 6. Assessment and Forward Path

### 6.1 What has been accomplished

1. **The tail extinction problem has been reduced to two explicit, testable hypotheses** ($\mathrm{H}_{\mathrm{gap}}$ and $\mathrm{H}_{\mathrm{collapse}}$) via a rigorous conditional theorem (Proposition 8.35).
2. **The logical relationships between analytical routes have been clarified.** Route 1B is the analytically cleanest; Route 1A is the most powerful but analytically hardest.
3. **The remaining analytical gaps are precisely identified.** No hand-waving, no circular reasoning, no imported continuous-ODE theorems.

### 6.2 What remains open

| Open Problem | Difficulty | Connection |
| :--- | :--- | :--- |
| Prove $\mathrm{H}_{\mathrm{gap}}(11)$ | Medium–Hard | Related to Gate 2 (continuum spectral structure) |
| Prove $\mathrm{H}_{\mathrm{collapse}}(j)$ for $j \le 10$ | Medium | Implied by $\mathrm{H}_{\mathrm{tunnel}}$ or Gate 2 degeneracy |
| Prove $\mathrm{H}_{\mathrm{tunnel}}(j)$ | Hard | Open problem in discrete semiclassical analysis |
| Extend to $T(N) \to \infty$ | Hard | Requires Route 1A or 1D |

### 6.3 Recommended next steps

1. **Record Proposition 8.35 in Paper NR2** as a conditional theorem with explicit epistemic labeling.
2. **Continue `cell95` investigation** of $g_{11}(N)$ at high $N$ to strengthen the empirical foundation for $\mathrm{H}_{\mathrm{gap}}(11)$.
3. **Defer `cell98`** until the analytical landscape is fully mapped and a specific falsification/verification target is formulated.
4. **Investigate whether Gate 2 construction (Friedrichs extension, strong resolvent convergence) can directly yield $\mathrm{H}_{\mathrm{gap}}$ and $\mathrm{H}_{\mathrm{collapse}}$** — if so, the tail extinction theorem becomes an immediate corollary of the continuum limit.

---

*Working note for the Connes–CvS research programme. All conditional theorems are labelled explicitly. Empirical observations are distinguished from mathematical proofs.*
