# Cell 125 Analytical Note: Variational Quasimode Codimension Bound — First-Principles Derivation

**Companion Document:** Analytical Research Note (Gate 1 / Milestone M-G1.5)  
**Status:** Working Analytical Note (Theoretical Derivation)  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.5: Variational Quasimode Codimension Bound](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell124.md`](file:///c:/data/github/connes-cvs-/cell124.md); [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md); [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)  
**Date:** September 2026  

---

## 1. Executive Summary & The Core Mandate

The review of [`cell124.md`](file:///c:/data/github/connes-cvs-/cell124.md) established the definitive strategic pivot for Gate 1:
$$\boxed{\text{We now need to prove a variational lower bound, not define a spectral projector.}}$$

1. **Retirement of Coordinate Blocks Ratified:**
   Coordinate truncation failed in Cell 123 because the physical bound-state eigenfunctions $u_j$ have non-compact Fourier support, bleeding tail mass $\mathcal{O}(\text{tail mass})$ into $\operatorname{span}\{e_M, \dots, e_N\}$ and producing artificial near-zero Ritz eigenvalues ($\lambda_{\min}(C_{12}) \approx 7.87 \times 10^{-6}$).
2. **The Tautology of $P_{\mathrm{cont}}$ Diagnosed:**
   Defining $P_{\mathrm{cont}} = I - \sum_{j=0}^{K-1} u_j u_j^T$ and observing $\lambda_{\min}(P_{\mathrm{cont}} Q P_{\mathrm{cont}}) = E_K$ is a restatement of the spectral theorem. It identifies the target subspace but does not prove continuum positivity independently.
3. **The Cell 125 Mandate:**
   To resolve the non-circularity problem, we must derive from first principles a rigorous, computable lower bound on the quadratic form:
   $$\lambda_{\min}(C) \equiv \inf_{\substack{v \perp \Phi \\ \|v\|=1}} \langle v, Q_{\mathrm{even}} v \rangle \ge c_* > 0,$$
   where $\Phi$ is an analytically constructed $K$-dimensional trial subspace ($K \approx 11$). By the Courant–Fischer min-max theorem:
   $$E_K \ge \lambda_{\min}(C) \ge c_* > 0$$
   then follows unconditionally, securing the continuum spectral denominator for Gate 1.

---

## 2. The Adapted Codimension Partitioning

Let $Q \equiv Q_{\mathrm{even}}^{(N)} \in \mathbb{R}^{(N+1) \times (N+1)}$ be the canonical even Galerkin matrix with ordered eigenvalues $E_0 \le E_1 \le \dots \le E_N$.

### Definition 2.1 (The Trial Subspace and Adapted Basis)
Let $\Phi \subset \mathbb{R}^{N+1}$ be a $K$-dimensional subspace ($K \ge 1$).
1. Choose an orthonormal basis for $\Phi$, represented as an isometry matrix:
   $$U_\Phi = [\phi_0, \dots, \phi_{K-1}] \in \mathbb{R}^{(N+1) \times K}, \qquad U_\Phi^T U_\Phi = I_K.$$
2. Let $U_{\Phi^\perp} \in \mathbb{R}^{(N+1) \times (N+1-K)}$ be an orthonormal basis for the orthogonal complement $V_\Phi \equiv \Phi^\perp$:
   $$U_{\Phi^\perp}^T U_{\Phi^\perp} = I_{N+1-K}, \qquad U_\Phi^T U_{\Phi^\perp} = 0.$$
3. The concatenated matrix $U \equiv [U_\Phi, U_{\Phi^\perp}]$ is an orthogonal matrix in $O(N+1)$.

### Proposition 2.1 (Block Partitioning & Exact Residual Identity)
In the adapted orthonormal basis $U$, the operator $Q$ is represented by the block symmetric matrix:
$$\widetilde{Q} \equiv U^T Q U = \begin{pmatrix} A & B^T \\ B & C \end{pmatrix},$$
where:
- $A \equiv U_\Phi^T Q U_\Phi \in \mathbb{R}^{K \times K}$ is the Rayleigh–Ritz compression to $\Phi$.
- $C \equiv U_{\Phi^\perp}^T Q U_{\Phi^\perp} \in \mathbb{R}^{(N+1-K) \times (N+1-K)}$ is the compression to $\Phi^\perp$.
- $B \equiv U_{\Phi^\perp}^T Q U_\Phi \in \mathbb{R}^{(N+1-K) \times K}$ is the cross-coupling block.

Furthermore, let $R \equiv Q U_\Phi - U_\Phi A$ be the residual matrix of the trial subspace $\Phi$. Then:
1. $U_\Phi^T R \equiv 0$, so $R$ has columns lying strictly in $\Phi^\perp$.
2. The coupling block $B$ satisfies the exact identity:
   $$B = U_{\Phi^\perp}^T R, \qquad R = U_{\Phi^\perp} B.$$
3. The operator 2-norm of the coupling block is identically the residual norm of the trial subspace:
   $$\|B\|_2 = \|R\|_2 \equiv \varepsilon.$$

*Proof.* By direct substitution: $U_\Phi^T R = U_\Phi^T Q U_\Phi - U_\Phi^T U_\Phi A = A - I_K A = 0$. Since $U U^T = U_\Phi U_\Phi^T + U_{\Phi^\perp} U_{\Phi^\perp}^T = I$, we have $R = U_{\Phi^\perp} U_{\Phi^\perp}^T R = U_{\Phi^\perp} B$. Because $U_{\Phi^\perp}$ is an isometry, $\|B\|_2 = \|U_{\Phi^\perp} B\|_2 = \|R\|_2 = \varepsilon$. $\blacksquare$

---

## 3. First-Principles Derivation of the Variational Lower Bound

We now derive the exact relation between the quadratic form on $\Phi^\perp$ and the exact eigenvalues of $Q$.

### Theorem 3.1 (Subspace Angle Rayleigh Bound)
Let $P_K \equiv \sum_{j=0}^{K-1} u_j u_j^T$ be the exact spectral projector for the lowest $K$ eigenvalues $\{E_0, \dots, E_{K-1}\}$ of $Q$, and let $P_K^\perp \equiv I - P_K = \sum_{j=K}^N u_j u_j^T$ be the spectral projector for the continuum sector $j \ge K$.
Let $\theta_{\max} \in [0, \pi/2]$ be the maximum principal angle between $\Phi$ and the exact bound-state invariant subspace $\mathcal{E}_K \equiv \operatorname{Ran}(P_K)$:
$$\sin\theta_{\max} \equiv \|P_K^\perp U_\Phi\|_2 = \|P_K U_{\Phi^\perp}\|_2 = \|P_K - P_\Phi\|_2.$$
Then for every unit vector $v \in \Phi^\perp$ ($\|v\|_2 = 1$):
$$\langle v, Q v \rangle \ge E_K \cos^2\theta_{\max} + E_0 \sin^2\theta_{\max}.$$
In particular, if $E_0 \ge 0$:
$$\lambda_{\min}(C) = \inf_{\substack{v \perp \Phi \\ \|v\|=1}} \langle v, Q v \rangle \ge E_K \cos^2\theta_{\max} = E_K (1 - \sin^2\theta_{\max}).$$

*Proof.* Decompose any unit vector $v \in \Phi^\perp$ into its exact spectral components:
$$v = v_b + v_c, \qquad v_b \equiv P_K v, \quad v_c \equiv P_K^\perp v.$$
Since $P_K$ and $P_K^\perp$ are mutually orthogonal projectors:
$$\|v_b\|_2^2 + \|v_c\|_2^2 = \|v\|_2^2 = 1, \qquad \langle v_b, Q v_c \rangle = 0.$$
Evaluating the quadratic form:
$$\langle v, Q v \rangle = \langle v_b, Q v_b \rangle + \langle v_c, Q v_c \rangle \ge E_0 \|v_b\|_2^2 + E_K \|v_c\|_2^2.$$
Substituting $\|v_c\|_2^2 = 1 - \|v_b\|_2^2$:
$$\langle v, Q v \rangle \ge E_K - (E_K - E_0) \|v_b\|_2^2.$$
Since $v \in \Phi^\perp$, $v = U_{\Phi^\perp} w$ for some $w \in \mathbb{R}^{N+1-K}$ with $\|w\|_2 = 1$. Thus:
$$\|v_b\|_2 = \|P_K U_{\Phi^\perp} w\|_2 \le \|P_K U_{\Phi^\perp}\|_2 = \sin\theta_{\max}.$$
Therefore $\|v_b\|_2^2 \le \sin^2\theta_{\max}$, which yields:
$$\langle v, Q v \rangle \ge E_K - (E_K - E_0) \sin^2\theta_{\max} = E_K (1 - \sin^2\theta_{\max}) + E_0 \sin^2\theta_{\max} = E_K \cos^2\theta_{\max} + E_0 \sin^2\theta_{\max}.$$
Taking the infimum over all unit vectors $v \in \Phi^\perp$ completes the proof. $\blacksquare$

---

## 4. Breaking the Non-Circularity Dilemma

Theorem 3.1 is an exact identity, but as noted in our review of Cell 124, bounding $\sin\theta_{\max}$ via the classical Davis–Kahan theorem:
$$\sin\theta_{\max} \le \frac{\|R\|_2}{\delta} = \frac{\varepsilon}{E_K - \lambda_{\max}(A)}$$
contains $E_K$ in the denominator $\delta$, which creates apparent circularity if one wishes to use the inequality to prove $E_K \ge c_* > 0$.

We now present two rigorous methods to break this circularity.

### Method 1: The Feshbach/Schur Determinantal Formulation & Epistemic Audit

Recall the block decomposition from Proposition 2.1:
$$\widetilde{Q} = \begin{pmatrix} A & B^T \\ B & C \end{pmatrix}, \qquad \|B\|_2 = \varepsilon.$$
Let $\lambda_{\max}(A) \le \bar{\mu} \ll 1$ be the maximum Rayleigh–Ritz energy on the trial subspace $\Phi$.

*Epistemic Audit of Resolvent Invertibility (Circularity Diagnosed):*  
One might be tempted to claim that if $Q - \lambda I$ is non-singular for all $\lambda \in (\bar{\mu} + \varepsilon, \gamma)$, then $E_K \ge \gamma$.  
However, as established in the review of Cell 125, this argument is **circular**:
1. Assuming that $(Q - \lambda I)^{-1}$ has no poles throughout an interval $(\bar{\mu} + \varepsilon, \gamma)$ assumes from the start that no eigenvalues exist in that interval—which is the very conclusion sought.
2. Furthermore, Kato–Temple / Bauer–Fike bounds only guarantee that each eigenvalue of $A$ has *some* eigenvalue of $Q$ within distance $\varepsilon$. It does not by itself prove that these eigenvalues correspond to the *first* $K$ eigenvalues of $Q$ without an independent multiplicity/order argument.
Consequently, Method 1 does not break the circularity and is retired as a proof of the gap.

---

### Method 2: Direct Operator Form Domination on $\Phi^\perp$ (The Open Analytical Target)

The true non-circular route is to prove that the compression $C = U_{\Phi^\perp}^T Q U_{\Phi^\perp}$ satisfies $\lambda_{\min}(C) \ge c_* > 0$ directly from the quadratic form:
$$\langle v, Q_{\mathrm{even}} v \rangle = \langle v, Q_{\mathrm{arch}} v \rangle + \langle v, Q_{\mathrm{prime}} v \rangle + \langle v, Q_{\mathrm{pole}} v \rangle \quad (\forall v \in \Phi^\perp).$$

Why did coordinate truncation $v \in \operatorname{span}\{e_{12}, \dots, e_N\}$ fail in Cell 123?
- In coordinate space, the basis vectors $e_m = \sqrt{2/L} \cos(2\pi m t / L)$ for $m \ge 12$ are unconstrained inside the potential well $[0, L]$.
- Specific linear combinations of high cosines can localize at the positions of the prime spikes $t_p = \log(p^k)$, exciting the negative eigenvalues of $Q_{\mathrm{prime}}$ ($\lambda_{\min}(C_{\mathrm{prime}}) \approx -2.373$) to cancel the Archimedean floor ($C_{\mathrm{arch}} \approx +1.553$).

### Open Proposition 4.2 (Subspace-Constrained Delocalization & Form Domination Target)
Suppose $\Phi$ contains:
1. An orthonormal basis for the low coordinate modes: $\operatorname{span}\{e_0, \dots, e_M\} \subset \Phi$ ($M \ge 3$, where $h_+(a_m) \ge 0.386$ for $m > M$).
2. The analytic well quasimodes $\{\phi_k^{\mathrm{well}}\}_{k=0}^{K-1}$ capturing the confined states of the potential well $V_{\mathrm{conf}}(t)$.

*The Central Open Analytical Question:*  
Does orthogonality to $\Phi$ ($v \in \Phi^\perp$) mathematically force $v$ to be a **delocalized scattering state** with bounded local amplitude at the prime coordinates:
$$\sup_{t \in [0, L]} |v(t)|^2 \stackrel{?}{\le} C_{\mathrm{deloc}} \|v\|^2?$$
*Epistemic Warning:*  
Orthogonality to a finite collection of localized well states does **not automatically imply pointwise delocalization**. High-energy scattering states could in principle still concentrate locally. Establishing an explicit, independent bound on $\langle v, Q_{\mathrm{prime}} v \rangle$ on $\Phi^\perp$ without assuming the continuum spectral gap is the exact analytical target for **Cell 126**.

---

## 5. The Approximate Invariant Subspace Separation Theorem

### Theorem 5.1 (Conditional Subspace Gap Separation Theorem)
Let $Q \in \mathbb{R}^{(N+1) \times (N+1)}$ be real symmetric.
Let $\Phi$ be a $K$-dimensional subspace with orthonormal basis $U_\Phi \in \mathbb{R}^{(N+1) \times K}$.
Let $A = U_\Phi^T Q U_\Phi$, and let $R = Q U_\Phi - U_\Phi A$ with $\|R\|_2 = \varepsilon$.
Suppose:
1. **Low-Energy Subspace Floor:** The Rayleigh–Ritz spectrum of $\Phi$ is strictly bounded:
   $$\lambda_{\max}(A) \le \bar{\mu}.$$
2. **Quasimode Residual:** $\varepsilon < \frac{c_* - \bar{\mu}}{2}$ for some target continuum threshold $c_* > \bar{\mu}$.
3. **Codimension Coercivity Hypothesis:** On the orthogonal complement $\Phi^\perp$, the operator satisfies:
   $$\lambda_{\min}(C) \equiv \inf_{v \perp \Phi, \|v\|=1} \langle v, Q v \rangle \ge c_*.$$

Then:
1. The spectrum of $Q$ splits into two disjoint components:
   $$\sigma(Q) = \sigma_{\mathrm{bound}}(Q) \cup \sigma_{\mathrm{cont}}(Q),$$
   where:
   $$\sigma_{\mathrm{bound}}(Q) \subset [\lambda_{\min}(A) - \varepsilon, \bar{\mu} + \varepsilon], \qquad \sigma_{\mathrm{cont}}(Q) \subset [c_* - \varepsilon^2/(c_* - \bar{\mu}), \infty).$$
2. The exact $(K+1)$-th eigenvalue satisfies:
   $$E_K \ge c_* - \frac{\varepsilon^2}{c_* - \bar{\mu}} > \bar{\mu} + \varepsilon \ge E_{K-1}.$$
3. The spectral gap separating the bound-state cluster from the continuum is at least:
   $$E_K - E_{K-1} \ge (c_* - \bar{\mu}) - \varepsilon - \frac{\varepsilon^2}{c_* - \bar{\mu}} > 0.$$

*Proof.* In the adapted basis $U = [U_\Phi, U_{\Phi^\perp}]$, write $\widetilde{Q} = \begin{pmatrix} A & B^T \\ B & C \end{pmatrix}$ with $\|B\|_2 = \varepsilon$.
For any $\lambda \in (\bar{\mu} + \varepsilon, c_* - \varepsilon^2/(c_* - \bar{\mu}))$, we examine the Schur complement:
$$S(\lambda) = C - \lambda I - B (A - \lambda I)^{-1} B^T.$$
Since $\lambda > \bar{\mu} + \varepsilon$, we have $\lambda I - A \succeq (\lambda - \bar{\mu}) I \succ 0$, so $(A - \lambda I)^{-1}$ is negative definite with:
$$\|(A - \lambda I)^{-1}\|_2 \le \frac{1}{\lambda - \bar{\mu}}.$$
Therefore:
$$- B (A - \lambda I)^{-1} B^T \succeq 0.$$
This means the second-order term is **positive semidefinite**!
Consequently, under Hypothesis 3 ($C \succeq c_* I$):
$$S(\lambda) \succeq C - \lambda I \succeq (c_* - \lambda) I \succ 0.$$
Since $S(\lambda)$ is strictly positive definite and $A - \lambda I$ is strictly negative definite, $\widetilde{Q} - \lambda I$ is non-singular for all $\lambda$ in this open interval.
Therefore, no eigenvalue of $Q$ can lie in $(\bar{\mu} + \varepsilon, c_* - \varepsilon^2/(c_* - \bar{\mu}))$.
Since the lowest $K$ eigenvalues of $Q$ lie below $\bar{\mu} + \varepsilon$, the $(K+1)$-th eigenvalue satisfies $E_K \ge c_* - \frac{\varepsilon^2}{c_* - \bar{\mu}}$. $\blacksquare$

*Epistemic Limitation of Theorem 5.1:*  
Theorem 5.1 is a **conditional transfer theorem**: it proves that *if* $C \succeq c_* I$, *then* the exact spectrum separates accordingly. It does **not** prove $C \succeq c_* I$ independently, because Hypothesis 3 assumes $\lambda_{\min}(C) \ge c_*$.

---

## 6. Significance & Strategic Mandate for Cell 126

Cell 125 has brought us to a clean, rigorous boundary:

1. **What is Established:**
   - The exact coupling identity $\|B\|_2 = \|R\|_2 = \varepsilon$ (Proposition 2.1).
   - The exact subspace angle lower bound $\langle v, Q v \rangle \ge E_K \cos^2\theta_{\max} + E_0 \sin^2\theta_{\max}$ on $\Phi^\perp$ (Theorem 3.1).
   - The conditional spectral-separation transfer theorem (Theorem 5.1).
2. **The Remaining Analytical Obstruction:**
   Proving that $C \equiv U_{\Phi^\perp}^T Q U_{\Phi^\perp} \succeq c_* I > 0$ directly from the quadratic form without presupposing $E_K$.
3. **The Core Mandate for Cell 126:**
   $$\boxed{\text{Gate 1 Denominator Problem} \quad \Longrightarrow \quad \text{Prove Delocalization / Control of } Q_{\mathrm{prime}} \text{ on } \Phi^\perp.}$$
   Cell 126 will investigate whether the indefinite prime form $|\langle v, Q_{\mathrm{prime}} v \rangle| \le C_{\mathrm{prime}} \|v\|^2$ can be explicitly bounded on $\Phi^\perp$ against the positive Archimedean background $\langle v, Q_{\mathrm{arch}} v \rangle \ge C_{\mathrm{arch}} \|v\|^2$ with $C_{\mathrm{arch}} - C_{\mathrm{prime}} > 0$, derived entirely independently of $E_K$.

---

## References
- [`cell124.md`](file:///c:/data/github/connes-cvs-/cell124.md) — Spectral-Subspace Projection & Min-Max Framework
- [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md) — Operator Decomposition Audit & Coordinate Falsification
- [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md) — Bound-State-to-Continuum Transition & Ritz Gaps
- [`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md) — Strategic Roadmap (Milestone M-G1.5)
- [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) — Semiclassical Barrier Mechanics
