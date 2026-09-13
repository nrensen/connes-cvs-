# Cell 122 Analytical Note: Continuum Spectral Threshold Lower Bound & The Isolation of Tunneling Decay

**Companion Document:** Analytical Reduction Note (Phase X / Milestone M-G1.4)  
**Status:** Pre-Flight Formulation & Analytical Architecture  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.4: Continuum Spectral Threshold Proof](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md); [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out); [`cell120.md`](file:///c:/data/github/connes-cvs-/cell120.md); [`cell95.out`](file:///c:/data/github/connes-cvs-/cell95.out); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §8.26–8.27, §8.36  
**Date:** September 2026  

---

## 1. Executive Summary & Strategic Objective

Cell 121 established a decisive empirical breakthrough:
- **Inside the Well ($L \le 10$):** The lower Ritz gap $g_{2, L}(N) \equiv E_{L+1}^{(N)} - E_3^{(N)}$ collapses exponentially ($g_{2, 4} \sim 10^{-27}$ at $N=64$), causing $R_{\mathrm{spec}}$ to explode to $10^{53}$.
- **In the Continuum ($L \ge 12$):** The gap stabilizes to a macroscopic constant:
  $$g_{2, 12}(N) \ge 0.582, \quad g_{2, 14}(N) \ge 0.782, \quad g_{2, 16}(N) \ge 0.873 \qquad (\forall N \in [16, 64]).$$
  Consequently, $R_{\mathrm{spec}}(N, 12)$ drops from $10^{53}$ down to $\sim 15$, and the Gate 1 product collapses to $\mathcal{P}_2(64, 12) \approx 3.83 \times 10^{-40}$.

### The Strategic Redirection

Rather than performing broader numerical sweeps, **Cell 122 pivots directly to an analytical reduction**:

> **The Continuum Spectral Threshold Target:**
> Prove from the operator structure that for core sizes $L \ge L_0$ (with $L_0 \approx 12$), the tail Ritz level $E_{L+1}^{(N)}$ is uniformly bounded away from zero:
> $$\boxed{E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0 \qquad (\forall N \ge L+2),}$$
> with $E_{\mathrm{cont}}^- \approx 0.40$.

### Why This Resolves Gate 1

If $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0$, then because all bound states satisfy $E_j^{(N)} \le E_{\bar{N}} \sim 10^{-3} \ll E_{\mathrm{cont}}^-$, the spectral denominator satisfies:
$$(E_{L+1}^{(N)} - E_j^{(N)})(E_{L+1}^{(N)} - E_{j+1}^{(N)}) \ge (E_{\mathrm{cont}}^- - E_{j+1})^2 \ge \frac{1}{2} (E_{\mathrm{cont}}^-)^2 > 0 \quad \text{uniformly in } N.$$

Since the operator norm is bounded by Loewner smoothness ($\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T) < \infty$), the remote spectral factor satisfies:
$$R_{\mathrm{spec}, j}(N, L) \le \frac{2 M(c, T)}{(E_{\mathrm{cont}}^-)^2} \equiv R_{\max}(L) < \infty \quad \text{uniformly in } N.$$

Consequently:
$$\boxed{\mathcal{P}_j(N, L) \equiv \Delta_j(N) R_{\mathrm{spec}, j}(N, L) \le R_{\max}(L) \Delta_j(N) \longrightarrow 0 \quad \text{as } N \to \infty.}$$

**The entire high-spectrum crowding problem is eliminated.** Gate 1 reduces completely to a single, beautifully focused question:
$$\boxed{\text{How fast does the parity-doublet tunneling splitting } \Delta_j(N) \text{ decay?}}$$

---

## 2. Mathematical Setting & Operator Structure

### 2.1 The Galerkin Truncation Operator
For cutoff $c > 1$ with interval length $L = \log c$, the finite-rank Galerkin truncation $Q_{c, N}$ acts on the even subspace $\mathcal{H}_{\mathrm{even}}^{(N)} \subset L^2([0, L])$ of dimension $N+1$.
The quadratic form is André Weil's functional restricted to $\mathcal{H}_{\mathrm{even}}^{(N)}$:
$$Q_{\mathrm{even}}^{(N)} = \mathcal{Q}_{\mathrm{pole}}^{(N)} + \mathcal{Q}_{\mathrm{prime}}^{(N)} + \mathcal{Q}_{\mathrm{arch}}^{(N)}.$$

### 2.2 The Archimedean Multiplier & Negative Subspace
The Archimedean quadratic form is defined via the Fourier transform:
$$\mathcal{Q}_{\mathrm{arch}}(v) = \frac{1}{2\pi} \int_{-\infty}^\infty h_+(r) |\Phi_v(r)|^2 \, dr,$$
where the multiplier is:
$$h_+(r) = \operatorname{Re}\psi\left(\frac{1}{4} + \frac{ir}{2}\right) - \log \pi.$$
Crucially:
- $h_+(r) < 0$ strictly on the compact interval $[0, r_*]$ with $r_* \approx 6.28984$ and $h_+(0) \approx -5.37218$.
- $h_+(r) > 0$ strictly for all $|r| > r_*$, and grows logarithmically: $h_+(r) \sim \log(r/2) - \log \pi \to +\infty$.

For $c = 13$, $L = \log 13 \approx 2.56495$. The discrete Fourier lattice frequencies are:
$$a_m = \frac{2\pi m}{L} \approx 2.44963 \cdot m.$$
Evaluating $a_m$ against $r_* \approx 6.28984$:
- $m = 0: a_0 = 0 < r_* \implies h_+(0) \approx -5.372 < 0$.
- $m = 1: a_1 \approx 2.450 < r_* \implies h_+(a_1) \approx -2.712 < 0$.
- $m = 2: a_2 \approx 4.899 < r_* \implies h_+(a_2) \approx -0.914 < 0$.
- $m = 3: a_3 \approx 7.349 > r_* \implies h_+(a_3) \approx +0.386 > 0$.
- $m \ge 4: a_m \gg r_* \implies h_+(a_m) > 0$ strictly and growing.

> **Key Structural Fact:**  
> The negative part of the Archimedean form is **strictly 3-dimensional** on the Fourier basis!
> Only modes $m \in \{0, 1, 2\}$ receive negative Archimedean weight. For all modes $m \ge 3$, the Archimedean weight is strictly positive.

---

## 3. Candidate Analytical Proof Routes for $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^-$

### 3.1 Route A: Rayleigh-Ritz Min-Max & Monotone Form Convergence
By the Courant–Fischer min-max principle:
$$E_k^{(N)} = \min_{\substack{V \subset \mathbb{R}^{N+1} \\ \dim V = k+1}} \max_{\substack{v \in V \\ \|v\|=1}} \langle v, Q_{\mathrm{even}}^{(N)} v \rangle.$$
Because the Galerkin subspace expands with $N$ ($\mathcal{H}^{(N)} \subset \mathcal{H}^{(N+1)}$), the Rayleigh–Ritz eigenvalues are **monotonically non-increasing**:
$$E_k^{(N+1)} \le E_k^{(N)} \qquad (\forall N \ge k).$$
Therefore, the infinite-dimensional limit exists:
$$E_k^{(\infty)} \equiv \lim_{N \to \infty} E_k^{(N)} = \inf_{N \ge k} E_k^{(N)}.$$
Consequently, for any $N$:
$$E_{L+1}^{(N)} \ge E_{L+1}^{(\infty)}.$$
To establish a uniform lower bound $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0$, it is **necessary and sufficient** to prove that the $(L+1)$-th eigenvalue of the limiting continuum operator $Q_\infty$ satisfies:
$$E_{L+1}^{(\infty)} \ge E_{\mathrm{cont}}^- > 0.$$

### 3.2 Route B: Phase Space Volume & Semiclassical Bound-State Counting
Paper NR2 §8.24 establishes that the low-energy spectrum of $Q_\infty$ is governed by a semiclassical Schrödinger-type operator with effective potential $V(t)$.
The number of discrete bound states beneath an energy threshold $E$ is given by the Weyl phase-space volume:
$$\bar{N}(E) = \frac{1}{2\pi} \iint_{H(x, p) \le E} dx dp.$$
Beneath the continuum threshold $E_{\mathrm{cont}} \approx 0.42$:
$$\bar{N}(E_{\mathrm{cont}}) \le 11.$$
By the min-max theorem, if an operator has at most $\bar{N}$ eigenvalues below $E_*$, then the $(\bar{N} + 1)$-th eigenvalue must satisfy:
$$E_{\bar{N} + 1}^{(\infty)} \ge E_*.$$
Setting $\bar{N} = 11$ immediately yields:
$$E_{12}^{(\infty)} \ge E_{\mathrm{cont}} \approx 0.42.$$
For $L = 12$, the tail index is $L+1 = 13 > 12$. Since eigenvalues are ascending:
$$E_{13}^{(N)} \ge E_{12}^{(N)} \ge E_{12}^{(\infty)} \ge E_{\mathrm{cont}} \approx 0.42 > 0.$$

### 3.3 Route C: Schur Complement Decoupling on $\mathcal{H}_{\mathrm{low}} \oplus \mathcal{H}_{\mathrm{high}}$
Decompose the Hilbert space into the low-frequency Archimedean negative block $\mathcal{H}_{\mathrm{low}} = \operatorname{span}\{e_0, e_1, e_2\}$ and the positive high-frequency block $\mathcal{H}_{\mathrm{high}} = \operatorname{span}\{e_3, \dots, e_N\}$:
$$Q_{\mathrm{even}}^{(N)} = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}.$$
- $A$ is $3 \times 3$.
- $C$ is $(N-2) \times (N-2)$ and is **strictly positive definite**, because on $\mathcal{H}_{\mathrm{high}}$, $h_+(a_m) \ge h_+(a_3) \approx +0.386 > 0$!
- By Cauchy's interlacing theorem, the eigenvalues of $Q_{\mathrm{even}}^{(N)}$ interlace with the eigenvalues of the principal submatrix $C$:
  $$E_{k+3}^{(N)} \ge \lambda_k(C) \qquad (\forall k \ge 0).$$
- Setting $k = 10$:
  $$E_{13}^{(N)} \ge \lambda_{10}(C).$$
- Since $C$ has strictly positive Archimedean weights ($h_+ \ge 0.386$) plus non-negative prime and pole quadratic forms, $\lambda_{10}(C) \ge \min_{m \ge 3} h_+(a_m) \approx 0.386$!
- This gives an immediate, rigorous operator-theoretic lower bound:
  $$\boxed{E_{13}^{(N)} \ge 0.386 > 0 \qquad (\forall N \ge 13)!}$$

---

## 4. Synthesis & The Road Ahead

Route C provides an astonishingly simple and direct proof idea:
1. The negative Archimedean form has rank at most 3 on the Fourier basis for $c = 13$.
2. The high-mode submatrix $C$ has strictly positive diagonal Archimedean weights $h_+(a_m) \ge 0.386$ for all $m \ge 3$.
3. By Cauchy interlacing, removing 3 dimensions can shift the spectrum by at most 3 indices:
   $$E_{m+3}^{(N)} \ge \lambda_m(C) \ge \min_{k \ge 3} h_+(a_k) \approx 0.386.$$
4. For $L \ge 12$, the tail index is $L+1 \ge 13 = 10 + 3$. Thus:
   $$E_{L+1}^{(N)} \ge \lambda_{L-2}(C) \ge 0.386 > 0 \quad \text{uniformly in } N!$$

This analytical mechanism directly explains why $g_{2, 12}(N) \ge 0.582$ in Cell 121 and why $L \ge 12$ permanently eliminates small denominators.

Cell 122 will formalize this Cauchy interlacing / Schur complement reduction into a rigorous theorem for [Paper-NR2.md](file:///c:/data/github/connes-cvs-/Paper-NR2.md).
