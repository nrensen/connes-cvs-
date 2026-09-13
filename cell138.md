# CELL 138 — EXACT VARIATIONAL DECOMPOSITION OF THE COUPLING GAIN, OPERATOR INCOMPATIBILITY, AND FALSIFICATION OF LOW-DIMENSIONAL REDUCTION

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.6 / Variational Lower Bound & Coupled Operator Geometry)  
**Target Proposition:** The Exact Variational Decomposition of the Coupling Gain and Diagnostic Evaluation of Low-Dimensional Reductions:
1. **Theorem 138.1 (Exact Algebraic Decomposition of the Coupling Gain):**
   For any self-adjoint restoring operator $K_{\mathrm{rest}}$ and potential operator $W_\perp$ on $\mathcal{B}_{11}^\perp$, with split Weyl bound $\mu_0^{\mathrm{split}} \equiv \lambda_{\min}(K_{\mathrm{rest}}) - \lambda_{\max}(W_\perp)$ and coupled ground state $\mu_0 \equiv \lambda_{\min}(K_{\mathrm{rest}} - W_\perp)$ attained at normalized minimizer $v_{\mathrm{bad}}$, the coupling gain satisfies the exact two-term partition:
   $$\Delta_{\mathrm{coupling}} \equiv \mu_0 - \mu_0^{\mathrm{split}} \equiv \Delta K + \Delta W,$$
   where:
   $$\Delta K \equiv \langle v_{\mathrm{bad}}, K_{\mathrm{rest}} v_{\mathrm{bad}} \rangle - \lambda_{\min}(K_{\mathrm{rest}}) \ge 0 \quad (\text{restoring energy excitation above ground}),$$
   $$\Delta W \equiv \lambda_{\max}(W_\perp) - \langle v_{\mathrm{bad}}, W_\perp v_{\mathrm{bad}} \rangle \ge 0 \quad (\text{potential well harvest sacrifice below supremum}).$$
2. **Proposition / Diagnostic 138.2 (Extremal Spectral Misalignment):**
   Let $x_0$ be the normalized ground state of $K_{\mathrm{rest}}$ ($\lambda_{\min}(K_{\mathrm{rest}}) \approx 2.9315$), and $y_0$ the dominant state of $W_\perp$ ($\lambda_{\max}(W_\perp) = 4.2605$).
   The principal overlap is:
   $$\cos^2 \theta_0 \equiv |\langle x_0, y_0 \rangle|^2 = 0.000000 \qquad (\theta_0 = 90.00^\circ).$$
   The extremal restoring direction and extremal well-harvesting direction are nearly orthogonal, demonstrating strong incompatibility of the two individual variational optima. The quantitative coupling mechanism, however, involves the full cross-geometry of their low/high spectral sectors.
3. **Diagnostic 138.3 (Falsification of Low-Dimensional Subspace Reduction):**
   Define the $d$-dimensional coupled subspace ($d \le k + \ell$):
   $$V_{k, \ell} \equiv \operatorname{span}\{x_0, \dots, x_{k-1}\} \oplus \operatorname{span}\{y_0, \dots, y_{\ell-1}\} \subset \mathcal{B}_{11}^\perp.$$
   The proposed $4$-mode reduction $V_{2, 2}$ **fails as an asymptotically stable approximation**: its captured mass decays from $97.58\%$ at $N=28$ to $83.70\%$ at $N=64$, and its eigenvalue error rises from $0.0621$ to $0.3103$. The coupling mechanism is a **collective spectral geometry**, not a fixed finite-mode interaction.
4. **Diagnostic 138.4 (Commutator Scale & Variance Equipartition):**
   The commutator norm $\|[K_{\mathrm{rest}}, W_\perp]\|_{\mathrm{op}}$ remains macroscopic over $24 \le N \le 64$, increasing from $1.79$ to $2.73$; whether it has a finite asymptotic limit remains open.
   On the minimizer $v_{\mathrm{bad}}$, the variances satisfy $\operatorname{Var}(K_{\mathrm{rest}}) = \operatorname{Var}(W_\perp) = 0.4985$ identically due to the eigenvalue relation $(K_{\mathrm{rest}} - \mu_0 I) v_{\mathrm{bad}} = W_\perp v_{\mathrm{bad}}$. A non-zero commutator by itself is not sufficient to establish a lower bound; the coupling gain requires quantitative spectral tradeoff analysis.

**Companion Computational Script:** [`cell138.py`](file:///c:/data/github/connes-cvs-/cell138.py)  
**Certified Execution Output:** [`cell138.out`](file:///c:/data/github/connes-cvs-/cell138.out) (commit `189b8b9`)

---

## 1. Executive Context: Moving from Splitting Failure to Structural Dissection

Cell 137 certified the Three-Level Hierarchy on the continuum constraint subspace $\mathcal{B}_{11}^\perp$:
1. **Level 1 (Zero-Mode Enclosure):** $\varepsilon_0(64) = 4.00\%$ captures $96.00\%$ of the zero mode, securing the kinetic floor $\ge -0.0642$.
2. **Level 2 (Independent Operator Splitting):**
   $$\mu_0^{\mathrm{split}} \equiv \lambda_{\min}(K_{\mathrm{rest}}) - \|W_{\perp \mathcal{B}}\|_{\mathrm{op}} = 2.931526 - 4.260495 = -1.328969 < -0.50.$$
3. **Level 3 (Coupled Variational Minimum):**
   $$\mu_0 = -0.48697922 > -0.50 \qquad (\text{Finite-}N\text{ margin } +0.01302078 > 0).$$

Independent operator splitting fails because it discards the relative geometry between the restoring operator $K_{\mathrm{rest}}$ and the step well $W_\perp$. The true coupled ground state outperforms the split Weyl bound by a macroscopic gain:
$$\Delta_{\mathrm{coupling}} \equiv \mu_0 - \mu_0^{\mathrm{split}} = +0.841990.$$

The objective of **Cell 138** is to analyze the variational origin of this $+0.8420$ gap: decomposing it into exact constituent energy penalties, testing the alignment of extremal modes, and empirically evaluating whether the interaction can be truncated to a low-dimensional effective Hamiltonian.

---

## 2. Theorem 138.1: Exact Variational Decomposition of the Coupling Gain

Let $\mathcal{H}_q = \mathbb{R}^q$ with $q = N - 10$ represent the continuum subspace $(\mathcal{B}_{11}^{(N)})^\perp$.
Let $K \equiv \widehat{K}_{\mathrm{rest}}$ and $W \equiv \widehat{W}_\perp$ be real symmetric $q \times q$ matrices representing the restoring stiffness and projected step well, respectively.

Let $\omega_0 \equiv \lambda_{\min}(K)$ and $\nu_0 \equiv \lambda_{\max}(W)$.
The split Weyl lower bound is $\mu_0^{\mathrm{split}} = \omega_0 - \nu_0$.
The coupled competition operator is $Q = K - W$. Let $\mu_0 = \lambda_{\min}(Q)$, with corresponding normalized eigenvector $v_{\mathrm{bad}} \in \mathcal{H}_q$ ($\|v_{\mathrm{bad}}\|_2 = 1$):
$$\mu_0 = \langle v_{\mathrm{bad}}, K v_{\mathrm{bad}} \rangle - \langle v_{\mathrm{bad}}, W v_{\mathrm{bad}} \rangle.$$

### Proof of the Exact Two-Term Identity
Subtracting the split Weyl bound:
$$\mu_0 - \mu_0^{\mathrm{split}} = \left( \langle v_{\mathrm{bad}}, K v_{\mathrm{bad}} \rangle - \langle v_{\mathrm{bad}}, W v_{\mathrm{bad}} \rangle \right) - (\omega_0 - \nu_0) = \left( \langle v_{\mathrm{bad}}, K v_{\mathrm{bad}} \rangle - \omega_0 \right) + \left( \nu_0 - \langle v_{\mathrm{bad}}, W v_{\mathrm{bad}} \rangle \right).$$

Define:
$$\Delta K \equiv \langle v_{\mathrm{bad}}, K_{\mathrm{rest}} v_{\mathrm{bad}} \rangle - \lambda_{\min}(K_{\mathrm{rest}}) \ge 0,$$
$$\Delta W \equiv \lambda_{\max}(W_\perp) - \langle v_{\mathrm{bad}}, W_\perp v_{\mathrm{bad}} \rangle \ge 0.$$

Since $\omega_0 = \min_{\|v\|=1} \langle v, K v \rangle$ and $\nu_0 = \max_{\|v\|=1} \langle v, W v \rangle$, both terms are unconditionally non-negative:
$$\boxed{\Delta_{\mathrm{coupling}} \equiv \mu_0 - \mu_0^{\mathrm{split}} \equiv \Delta K + \Delta W \ge 0.}$$

### Conceptual Status: Decomposition vs Mechanism
This identity does not by itself explain the numerical magnitude of the coupling gain; rather, it **quantifies the exact variational tradeoff** enforced on the minimizer:
- **$\Delta K$ (Restoring Excess):** How much restoring stiffness $v_{\mathrm{bad}}$ accepts above the lowest possible state $\omega_0$.
  At $N=64$: $\langle v, K v \rangle = 3.222804 \implies \Delta K = 3.222804 - 2.931526 = +0.291278$ ($34.59\%$ of the gain).
- **$\Delta W$ (Well Harvest Deficit):** How much well depth $v_{\mathrm{bad}}$ forfeits below the maximal possible harvest $\nu_0 = 4.260495$.
  At $N=64$: $R_W(v_{\mathrm{bad}}) = 3.709783 \implies \Delta W = 4.260495 - 3.709783 = +0.550712$ ($65.41\%$ of the gain).

The minimizer pays roughly **one-third in excess kinetic excitation** and **two-thirds in forfeited well depth**.

---

## 3. Proposition / Diagnostic 138.2: Extremal Spectral Misalignment

Let $\{x_j\}_{j=0}^{q-1}$ be the orthonormal eigenbasis of $K_{\mathrm{rest}}$ ($\omega_0 < \omega_1 \le \dots$), and $\{y_j\}_{j=0}^{q-1}$ the orthonormal eigenbasis of $W_\perp$ ($\nu_0 > \nu_1 \ge \dots$).

The principal alignment between the ground restoring mode $x_0$ and dominant well mode $y_0$ is:
$$\cos^2 \theta_0 \equiv |\langle x_0, y_0 \rangle|^2.$$

The 50-dps calculation establishes:
$$\cos^2 \theta_0 = 0.000000 \qquad (\theta_0 = 90.00^\circ) \quad \forall N \in [24, 64].$$

### Interpretation & Scope
The extremal restoring direction and extremal well-harvesting direction are nearly orthogonal, demonstrating strong incompatibility of the two individual variational optima.
However, this numerical orthogonality does not by itself constitute an explanation of the quantitative coupling gain.
In particular, the relevant geometry is not merely $x_0 \perp y_0$. The actual minimizer $v_{\mathrm{bad}}$ involves the full cross-geometry of their low/high spectral sectors:
- In the $K_{\mathrm{rest}}$ eigenbasis, $v_{\mathrm{bad}}$ is relatively concentrated: $x_0$ holds $73.75\%$ mass, and the lowest 4 modes hold $89.38\%$.
- In the $W_\perp$ eigenbasis, $v_{\mathrm{bad}}$ is broadly distributed: $y_0$ holds only $21.06\%$, and the top 4 modes hold only $49.29\%$.

The quantitative coupling mechanism therefore involves the **collective cross-geometry of the entire low-$K$ and high-$W$ spectral sectors**, not merely the isolated pair $(x_0, y_0)$.

---

## 4. Diagnostic 138.3: Falsification of Low-Dimensional Subspace Reduction

To test whether the continuum competition operator can be captured by a low-dimensional effective Hamiltonian, Cell 138 evaluated the coupled trial subspaces:
$$V_{k, \ell} \equiv \operatorname{span}\{x_0, \dots, x_{k-1}\} \oplus \operatorname{span}\{y_0, \dots, y_{\ell-1}\} \subset \mathcal{B}_{11}^\perp.$$

We track the captured mass $\rho(k, \ell) \equiv \|P_{V_{k, \ell}} v_{\mathrm{bad}}\|^2$ and the effective ground-state eigenvalue $\mu_0^{(k, \ell)} \equiv \lambda_{\min}(P_{V_{k, \ell}} Q_{\mathrm{comp}} P_{V_{k, \ell}})$.

### Numerical Findings for the $4$-Mode Subspace $V_{2, 2}$:
| $N$ | Captured Mass $\|P_{V_{2, 2}} v_{\mathrm{bad}}\|^2$ | Eigenvalue Error $\mu_0^{(2, 2)} - \mu_0$ | Assessment |
| :---: | :---: | :---: | :---: |
| 24 | $96.99\%$ | $+0.0659$ | Moderate |
| 28 | $97.58\%$ | $+0.0621$ | Moderate |
| 32 | $97.40\%$ | $+0.0681$ | Moderate |
| 40 | $96.07\%$ | $+0.1035$ | Deteriorating |
| 48 | $88.86\%$ | $+0.2297$ | Poor |
| **64** | **$83.70\%$** | **$+0.3103$** | **FALSIFIED** |

### Mathematical Conclusion
The $V_{2, 2}$ reduction fails the pre-flight success criteria ($\ge 99\%$ capture and error $\le 0.05$).
Crucially, its performance **deteriorates steadily as $N$ increases**: at $N=64$, more than $16.3\%$ of the minimizer's mass escapes $V_{2, 2}$, producing a large eigenvalue error of $+0.3103$.

**Strategic Lesson:** The coupling mechanism is **not an interaction between a small, fixed set of modes**. Attempting to truncate the problem to a tiny effective Hamiltonian (e.g. $2 \times 2$ or $4 \times 4$) fails because the well deficit $\Delta W$ is distributed across many high-frequency potential modes.

---

## 5. Diagnostic 138.4: Commutator Scale & Variance Equipartition

### 5.1 Commutator Norm Growth
The operator commutator $C = [K_{\mathrm{rest}}, W_\perp] \in \operatorname{Skew}(\mathcal{B}_{11}^\perp)$ has spectral norm:
$$\|[K_{\mathrm{rest}}, W_\perp]\|_{\mathrm{op}} = \sqrt{\lambda_{\max}(-C^2)}.$$
Across discrete dimensions:
$$\begin{array}{c|cccccc}
N & 24 & 28 & 32 & 40 & 48 & 64 \\
\hline
\|[K, W]\|_{\mathrm{op}} & 1.7930 & 2.2322 & 2.4444 & 2.5246 & 2.5644 & 2.7345
\end{array}$$

The commutator norm remains macroscopic over $24 \le N \le 64$, increasing from $1.79$ to $2.73$; whether it has a finite asymptotic limit remains open. Furthermore, a non-zero commutator by itself is not sufficient to establish a positive lower bound above $-1/2$; an analytical connection between operator incompatibility, spectral dispersion, and the coupled minimum is required.

### 5.2 Variance Equipartition on $v_{\mathrm{bad}}$
The standard deviations of $K_{\mathrm{rest}}$ and $W_\perp$ evaluated on $v_{\mathrm{bad}}$ satisfy:
$$\sigma_K \equiv \sqrt{\langle v, K^2 v \rangle - \langle v, K v \rangle^2} = \sigma_W \equiv \sqrt{\langle v, W^2 v \rangle - \langle v, W v \rangle^2} = 0.706048 \quad (\text{at } N=64).$$

**Proof of Identity:**
The minimizer satisfies $(K - W) v_{\mathrm{bad}} = \mu_0 v_{\mathrm{bad}}$, which can be rewritten as:
$$(K - \mu_0 I) v_{\mathrm{bad}} = W v_{\mathrm{bad}}.$$
Taking the squared norm of both sides:
$$\|(K - \mu_0 I) v_{\mathrm{bad}}\|^2 = \|W v_{\mathrm{bad}}\|^2 \implies \langle v, K^2 v \rangle - 2 \mu_0 \langle v, K v \rangle + \mu_0^2 = \langle v, W^2 v \rangle.$$
Substituting $\mu_0 = \langle v, K v \rangle - \langle v, W v \rangle$:
$$\langle v, K^2 v \rangle - \langle v, K v \rangle^2 = \langle v, W^2 v \rangle - \langle v, W v \rangle^2 \quad \Longrightarrow \quad \operatorname{Var}(K) \equiv \operatorname{Var}(W).$$
This is an exact algebraic consequence of the eigenvalue equation.
Furthermore, because $C = [K, W]$ is real skew-symmetric, $\langle v, C v \rangle = 0$ identically for any real vector. Table 5 provides a diagnostic of joint variance structure, not an uncertainty bound driven by a non-zero commutator expectation.

---

## 6. Certified 50-DPS Numerical Results

Executed via external compute node (`cell138.out`, commit `189b8b9`), runtime $170.10\text{ s}$ across $N \in [24, 64]$:

### Table 1: Exact Algebraic Decomposition of the Coupling Gain
$$\Delta_{\mathrm{coupling}} \equiv \mu_0 - \mu_0^{\mathrm{split}} \equiv \Delta K + \Delta W$$

| $N$ | $\mu_0^{\mathrm{split}}$ | $\mu_0$ | $\Delta_{\mathrm{coupling}}$ | $\Delta K$ | $\Delta W$ | $\% \Delta K$ | $\% \Delta W$ | Closure Error |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $-0.762703$ | $-0.00807113$ | $+0.754632$ | $0.220956$ | $0.533676$ | $29.28\%$ | $70.72\%$ | $1.34 \times 10^{-51}$ |
| 28 | $-1.199009$ | $-0.28445278$ | $+0.914556$ | $0.330996$ | $0.583561$ | $36.19\%$ | $63.81\%$ | $1.34 \times 10^{-51}$ |
| 32 | $-1.268002$ | $-0.38662701$ | $+0.881375$ | $0.310568$ | $0.570807$ | $35.24\%$ | $64.76\%$ | $1.07 \times 10^{-50}$ |
| 40 | $-1.312160$ | $-0.45755948$ | $+0.854601$ | $0.294572$ | $0.560029$ | $34.47\%$ | $65.53\%$ | $2.67 \times 10^{-51}$ |
| 48 | $-1.324539$ | $-0.47803442$ | $+0.846505$ | $0.292687$ | $0.553818$ | $34.58\%$ | $65.42\%$ | $1.34 \times 10^{-51}$ |
| **64** | **$-1.328969$** | **$-0.48697922$** | **$+0.841990$** | **$0.291278$** | **$0.550712$** | **$34.59\%$** | **$65.41\%$** | **$6.68 \times 10^{-51}$** |

---

### Table 2: Extremal Eigenvector Misalignment & Angular Separation

| $N$ | $\cos^2 \theta_0$ | Angle $\theta_0$ | $|\langle w_{\mathrm{bad}}, x_0 \rangle|^2$ | $|\langle w_{\mathrm{bad}}, y_0 \rangle|^2$ | Alignment Status |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $0.000000$ | $90.00^\circ$ | $85.96\%$ | $0.00\%$ | ORTHOGONAL |
| 28 | $0.000000$ | $90.00^\circ$ | $62.24\%$ | $0.00\%$ | ORTHOGONAL |
| 32 | $0.000000$ | $90.00^\circ$ | $69.26\%$ | $0.00\%$ | ORTHOGONAL |
| 40 | $0.000000$ | $90.00^\circ$ | $73.36\%$ | $0.00\%$ | ORTHOGONAL |
| 48 | $0.000000$ | $90.00^\circ$ | $73.79\%$ | $0.00\%$ | ORTHOGONAL |
| **64** | **$0.000000$** | **$90.00^\circ$** | **$73.75\%$** | **$0.00\%$** | **ORTHOGONAL** |

---

### Table 3: Minimizer Spectral Mass Distribution in $K_{\mathrm{rest}}$ and $W_\perp$ Bases

| $N$ | $K$-Mode 1 | $K$-Modes 1–2 | $K$-Modes 1–4 | $W$-Mode 1 | $W$-Modes 1–2 | $W$-Modes 1–4 |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $85.96\%$ | $87.04\%$ | $93.17\%$ | $58.34\%$ | $58.37\%$ | $61.60\%$ |
| 28 | $62.24\%$ | $88.86\%$ | $90.37\%$ | $45.99\%$ | $57.51\%$ | $66.33\%$ |
| 32 | $69.26\%$ | $88.77\%$ | $90.22\%$ | $38.53\%$ | $58.41\%$ | $58.50\%$ |
| 40 | $73.36\%$ | $88.94\%$ | $89.20\%$ | $33.05\%$ | $52.37\%$ | $58.36\%$ |
| 48 | $73.79\%$ | $77.79\%$ | $89.22\%$ | $29.07\%$ | $45.58\%$ | $58.88\%$ |
| **64** | **$73.75\%$** | **$74.38\%$** | **$89.38\%$** | **$21.06\%$** | **$31.46\%$** | **$49.29\%$** |

---

### Table 4: Low-Dimensional Subspace Reduction Models ($V_{1, 1}$ vs $V_{2, 2}$)

| $N$ | $V_{1, 1}$ Capture | $\mu_0^{(1, 1)}$ | $\mathrm{Err}(1, 1)$ | $V_{2, 2}$ Capture | $\mu_0^{(2, 2)}$ | $\mathrm{Err}(2, 2)$ | Asymptotic Model Status |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $96.70\%$ | $+0.061116$ | $0.069187$ | $96.99\%$ | $+0.057835$ | $0.065906$ | Sub-critical |
| 28 | $81.33\%$ | $+0.010924$ | $0.295376$ | $97.58\%$ | $-0.222362$ | $0.062091$ | Sub-critical |
| 32 | $83.32\%$ | $-0.088364$ | $0.298263$ | $97.40\%$ | $-0.318576$ | $0.068051$ | Sub-critical |
| 40 | $84.19\%$ | $-0.153904$ | $0.303655$ | $96.07\%$ | $-0.354066$ | $0.103494$ | Deteriorating |
| 48 | $83.29\%$ | $-0.157323$ | $0.320711$ | $88.86\%$ | $-0.248347$ | $0.229688$ | Failing |
| **64** | **$81.01\%$** | **$-0.136009$** | **$0.350971$** | **$83.70\%$** | **$-0.176728$** | **$0.310251$** | **FALSIFIED** |

---

### Table 5: Commutator Norm & Variance Diagnostics

| $N$ | $\|[K_{\mathrm{rest}}, W_\perp]\|_{\mathrm{op}}$ | $\sigma_K$ on $v_{\mathrm{bad}}$ | $\sigma_W$ on $v_{\mathrm{bad}}$ | Variance $\sigma^2$ |
| :---: | :---: | :---: | :---: | :---: |
| 24 | $1.793025$ | $0.650032$ | $0.650032$ | $0.422541$ |
| 28 | $2.232241$ | $0.696441$ | $0.696441$ | $0.485030$ |
| 32 | $2.444361$ | $0.709715$ | $0.709715$ | $0.503696$ |
| 40 | $2.524613$ | $0.708392$ | $0.708392$ | $0.501820$ |
| 48 | $2.564443$ | $0.707260$ | $0.707260$ | $0.500217$ |
| **64** | **$2.734524$** | **$0.706048$** | **$0.706048$** | **$0.498504$** |

---

## 7. Strategic Synthesis & Forward Path to Cell 139

Cell 138 establishes three foundational results:

1. **Exact Variational Tradeoff:**
   The $+0.841990$ coupling gain is partitioned into $\Delta K = 0.291278$ ($34.6\%$) restoring excess and $\Delta W = 0.550712$ ($65.4\%$) forfeited well harvest.
2. **Extremal Incompatibility:**
   The ground restoring mode $x_0$ and dominant well mode $y_0$ are mutually orthogonal ($\theta_0 = 90.00^\circ$).
3. **Decisive Falsification of Low-Dimensional Truncation:**
   The 4-mode subspace $V_{2, 2}$ fails asymptotically ($83.70\%$ capture, error $+0.3103$). The coupling is fundamentally an **infinite-dimensional collective spectral geometry**.

### Strategic Pivot for Cell 139
Instead of attempting larger fixed-dimensional trial spaces ($V_{3, 3}, V_{4, 4}$), the research programme should pivot to proving a **functional spectral tradeoff inequality**:
$$\Delta K \ge F(\Delta W) \qquad \forall v \in \mathcal{B}_{11}^\perp, \ \|v\|=1,$$
where:
$$\Delta K(v) = \sum_j (\omega_j - \omega_0) |\langle x_j, v \rangle|^2, \qquad \Delta W(v) = \sum_k (\nu_0 - \nu_k) |\langle y_k, v \rangle|^2.$$
The cross-Gram matrix $O_{j, k} \equiv |\langle x_j, y_k \rangle|^2$ bridges the two spectral expansions, converting the coupling gain into a quantitative incompatibility theorem.
