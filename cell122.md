# Cell 122 Analytical Note: Continuum Spectral Threshold Lower Bound & The Three-Lemma Architecture

**Companion Document:** Analytical Reduction Note (Gate 1 / Milestone M-G1.4)  
**Status:** Calibrated Research Note / Theoretical Framework & Proof Obligations  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.4: Continuum Spectral Threshold Proof](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md); [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out); [`cell120.md`](file:///c:/data/github/connes-cvs-/cell120.md); [`cell95.out`](file:///c:/data/github/connes-cvs-/cell95.out); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §8.26–8.27, §8.36  
**Date:** September 2026  

---

## 1. Executive Summary & Strategic Objective

Cell 121 established a decisive empirical breakthrough in the Gate 1 architecture:
- **Inside the Well ($L \le 10$):** The lower Ritz gap $g_{2, L}(N) \equiv E_{L+1}^{(N)} - E_3^{(N)}$ collapses exponentially ($g_{2, 4} \sim 10^{-27}$ at $N=64$), causing $R_{\mathrm{spec}}$ to explode to $10^{53}$.
- **In the Continuum ($L \ge 12$):** The gap stabilizes to macroscopic values ($g_{2, 12} \ge 0.582, g_{2, 14} \ge 0.782, g_{2, 16} \ge 0.873$ across $N \in [16, 64]$), dropping $R_{\mathrm{spec}}$ to $\sim 15$ and collapsing the Gate 1 product to $\mathcal{P}_2(64, 12) \approx 3.83 \times 10^{-40}$.

### The Strategic Redirection

Rather than performing broader numerical sweeps, **Cell 122 pivots directly to the analytical mechanism**:
Can we prove from the operator structure that for core sizes $L \ge L_0$ (with $L_0 \approx 12$), the tail Ritz level $E_{L+1}^{(N)}$ is uniformly bounded away from zero:
$$\boxed{E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0 \qquad (\forall N \ge L+2),}$$
ideally with $E_{\mathrm{cont}}^- \approx 0.40$?

### Separation of Denominator and Numerator Theorems

Controlling the remote spectral factor $R_{\mathrm{spec}}(N, L)$ requires two distinct, independent mathematical components:
1. **The Denominator Theorem:**
   $$E_{L+1}^{(N)} - E_{j+1}^{(N)} \ge \delta > 0 \qquad (\forall N \ge L+2, \; L \ge L_0).$$
2. **The Numerator Theorem:**
   $$\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \le M(c, T) < \infty \qquad (\forall N \ge 1).$$

Together, these two theorems guarantee:
$$R_{\mathrm{spec}, j}(N, L) \le \frac{M(c, T)}{\delta^2} \equiv R_{\max}(L) < \infty \quad \text{uniformly in } N,$$
which completely eliminates high-spectrum crowding and reduces Gate 1 to the asymptotic decay rate of the parity-doublet tunneling splitting $\Delta_j(N)$.

---

## 2. The Archimedean Multiplier Sign Change & Representation Hygiene

### 2.1 The Three-Mode Negative Multiplier Window
For cutoff $c = 13$, the interval length is $L = \log 13 \approx 2.56495$. The discrete Fourier lattice frequencies are:
$$a_m = \frac{2\pi m}{L} \approx 2.44963 \cdot m.$$
Evaluating $a_m$ against the Archimedean sign-change point $r_* \approx 6.28984$ (where $h_+(r_*) = 0$):
- $m = 0: a_0 = 0 < r_* \implies h_+(0) \approx -5.37218 < 0$.
- $m = 1: a_1 \approx 2.450 < r_* \implies h_+(a_1) \approx -2.71184 < 0$.
- $m = 2: a_2 \approx 4.899 < r_* \implies h_+(a_2) \approx -0.91421 < 0$.
- $m = 3: a_3 \approx 7.349 > r_* \implies h_+(a_3) \approx +0.38604 > 0$.
- $m \ge 4: a_m \gg r_* \implies h_+(a_m) > 0$ strictly and growing logarithmically.

The Archimedean multiplier $h_+(r)$ changes sign cleanly between $m = 2$ and $m = 3$.

### 2.2 Critical Representation Distinction (Dense Operator vs Multiplier)
It is tempting to write $Q_{\mathrm{arch}} = \operatorname{diag}(h_+(a_0), h_+(a_1), \dots)$ and deduce that the high-mode block $m \ge 3$ is automatically coercive with lower bound $h_+(a_3) \approx 0.386$.
**This is an illegitimate conflation of continuous Fourier multipliers with finite-interval Galerkin matrices:**
- In the continuum on $\mathbb{R}$, $h_+(r)$ acts as a diagonal Fourier multiplier.
- In the finite-rank Galerkin model on $[0, L]$ (`connes_cvs/operator.py`), $Q_{\mathrm{arch}}$ is a **dense divided-difference matrix**:
  $$Q_{\mathrm{arch}}[m, n] = \frac{\psi_{\mathrm{arch}}(m) - \psi_{\mathrm{arch}}(n)}{m - n} \quad (m \ne n), \qquad Q_{\mathrm{arch}}[n, n] = \psi_{\mathrm{arch}}'(n).$$
- The basis functions $e_k(t) = e^{2\pi i k t / L}$ on $[0, L]$ do **not** diagonalize the singular integral operator on $\mathbb{R}$; finite boundary truncation introduces non-zero off-diagonal coupling.
- Therefore, $h_+(a_m) > 0$ for $m \ge 3$ does **not** automatically imply that the principal submatrix $Q_{\mathrm{arch}}|_{\mathcal{H}_{\mathrm{high}}}$ has eigenvalues bounded below by $0.386$. This must be proved as an operator inequality.

---

## 3. The Three-Lemma Proof Architecture for the Continuum Threshold

To establish $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0$ rigorously, the analytical argument is factored into three precise lemmas:

### Lemma A (Archimedean High-Frequency Coercivity)
Let $\mathcal{H}_{\mathrm{high}} = \operatorname{span}\{e_3, \dots, e_N\}$ be the subspace spanned by modes $m \ge 3$. Prove that the dense Galerkin submatrix $C_{\mathrm{arch}} = Q_{\mathrm{arch}}|_{\mathcal{H}_{\mathrm{high}}}$ satisfies:
$$\langle v, Q_{\mathrm{arch}} v \rangle \ge c_{\mathrm{arch}} \|v\|^2 \qquad (\forall v \in \mathcal{H}_{\mathrm{high}}),$$
for some strictly positive constant $c_{\mathrm{arch}} > 0$ (candidate: $c_{\mathrm{arch}} \approx h_+(a_3) \approx 0.386$).

### Lemma B (Positivity / Lower Bound of Arithmetic and Pole Sectors)
In `connes_cvs/operator.py`, the full Galerkin operator is $Q = Q_{\mathrm{arch}} + Q_{\mathrm{prime}} + Q_{\mathrm{pole}}$, where:
$$Q_{\mathrm{prime}}[m, n] = \frac{\psi_{\mathrm{prime}}(m) - \psi_{\mathrm{prime}}(n)}{m - n}, \qquad Q_{\mathrm{pole}}[m, n] = \frac{\psi_{\mathrm{pole}}(m) - \psi_{\mathrm{pole}}(n)}{m - n}.$$
Prove that the arithmetic and pole contributions on $\mathcal{H}_{\mathrm{high}}$ are positive semidefinite, or at least bounded below:
$$\langle v, (Q_{\mathrm{prime}} + Q_{\mathrm{pole}}) v \rangle \ge -c_{\mathrm{pert}} \|v\|^2 \qquad (\forall v \in \mathcal{H}_{\mathrm{high}}),$$
with $c_{\mathrm{pert}} < c_{\mathrm{arch}}$, so that the total high-mode block $C = Q_{\mathrm{even}}|_{\mathcal{H}_{\mathrm{high}}}$ satisfies:
$$C \succeq c_* I \qquad (c_* = c_{\mathrm{arch}} - c_{\mathrm{pert}} > 0).$$

### Lemma C (Finite-Rank Cauchy Interlacing & The Bound-State Obstruction)

Let $Q_{\mathrm{even}}^{(N)}$ be the $(N+1) \times (N+1)$ even Galerkin matrix, with eigenvalues indexed in standard 0-based ascending order:
$$E_0^{(N)} \le E_1^{(N)} \le \dots \le E_N^{(N)}.$$

#### The 3-Mode Submatrix ($M = 3$)
Partition $Q_{\mathrm{even}}^{(N)}$ by removing the 3 lowest modes $m \in \{0, 1, 2\}$, leaving the $(N-2) \times (N-2)$ principal submatrix $C_3 = Q_{\mathrm{even}}|_{\operatorname{span}\{e_3, \dots, e_N\}}$, with ascending eigenvalues $\lambda_0(C_3) \le \lambda_1(C_3) \le \dots \le \lambda_{N-3}(C_3)$.

By Cauchy's interlacing theorem for deleting 3 rows and columns:
$$E_k^{(N)} \le \lambda_k(C_3) \le E_{k+3}^{(N)} \qquad (\forall k \in \{0, \dots, N-3\}).$$

**1. The Fatal Upper Bound on $\lambda_{\min}(C_3)$:**
Setting $k = 0$:
$$\boxed{\lambda_{\min}(C_3) = \lambda_0(C_3) \le E_3^{(N)} \approx 2.98 \times 10^{-38} \ll 0.386.}$$
Because $E_3$ is the 4th bound state of the potential well, its energy is exponentially close to zero. Consequently, **$C_3$ cannot be coercive with constant $0.386$**. Deleting only 3 basis functions leaves the remaining $\sim 8$ bound states ($E_3, \dots, E_{10}$) supported inside $C_3$, forcing $\lambda_0(C_3), \dots, \lambda_7(C_3)$ to cluster near zero.

**2. The Interior Lower Bound on $E_{13}^{(N)}$:**
Setting $k = 10$:
$$E_{13}^{(N)} \ge \lambda_{10}(C_3) \ge E_{10}^{(N)}.$$
This bounds $E_{13}$ from below by the **11th eigenvalue** $\lambda_{10}(C_3)$ (which is near the continuum threshold), **not** by the spectral bottom $\lambda_{\min}(C_3)$.

#### The Core-Cut Submatrix ($M = 12$)
To obtain a truly coercive submatrix whose minimum eigenvalue is macroscopic, one must project out the **entire bound-state core**:
$$C_{12} = Q_{\mathrm{even}}|_{\operatorname{span}\{e_{12}, \dots, e_N\}}.$$
On this subspace, all $\sim 11$ bound states have been eliminated. Cauchy interlacing for deleting $M = 12$ rows/columns gives:
$$E_{12+k}^{(N)} \ge \lambda_k(C_{12}) \ge \lambda_{\min}(C_{12}) \qquad (\forall k \ge 0).$$
Setting $k = 1$ (the first tail mode for $L = 12$):
$$\boxed{E_{13}^{(N)} \ge \lambda_1(C_{12}) \ge \lambda_{\min}(C_{12}).}$$
If $C_{12} \succeq c_{12} I > 0$, then $E_{13}^{(N)} \ge c_{12} > 0$ uniformly in $N$.

---

## 4. Epistemic Audit of Proof Obligations

| Proof Obligation | Mathematical Question | Status |
| :--- | :--- | :---: |
| **1. Coercivity of $C_3$ vs Bound-State Obstruction** | Does $\lambda_{\min}(C_3) \le E_3 \approx 10^{-38}$ hold, confirming that $C_3 \succeq 0.386 I$ is mathematically impossible? | **RESOLVED NEGATIVELY (Cauchy Bound)** |
| **2. Coercivity of $C_{12}$ on $\operatorname{span}\{e_{12}, \dots, e_N\}$** | Does the core-cut submatrix satisfy $C_{12} \succeq c_{12} I > 0$ with $c_{12} \approx 0.50$ uniformly across $N$? | **OPEN (Cell 123 Target)** |
| **3. Archimedean High-Mode Coercivity** | Does the dense divided-difference matrix $Q_{\mathrm{arch}}|_{\mathcal{H}_{\mathrm{high}}}$ satisfy $C_{\mathrm{arch}} \succeq c_{\mathrm{arch}} I > 0$? | **OPEN (Cell 123 Target)** |
| **4. Prime Form Sign Structure** | In `connes_cvs/operator.py`, $\psi_{\mathrm{prime}}$ carries an explicit minus sign. What is the signature of $Q_{\mathrm{prime}}$ on $\mathcal{H}_{\mathrm{high}}$? | **OPEN (Cell 123 Target)** |
| **5. Pole Form Sign Structure** | Is the zeta-pole divided-difference matrix $Q_{\mathrm{pole}}$ positive semidefinite on $\mathcal{H}_{\mathrm{high}}$? | **OPEN (Cell 123 Target)** |

---

## 5. Forward Mandate for Cell 123

To resolve these proof obligations directly from the code, **Cell 123** will execute an **operator decomposition audit**:
1. Extract the individual component matrices $Q_{\mathrm{arch}}, Q_{\mathrm{prime}}, Q_{\mathrm{pole}}$ using the exact definitions in `connes_cvs/operator.py`.
2. Compute the spectrum of $C_3$ on $\operatorname{span}\{e_3, \dots, e_N\}$ to verify $\lambda_{\min}(C_3) \le E_3$ and determine the individual component bounds.
3. Compute the spectrum of $C_{12}$ on $\operatorname{span}\{e_{12}, \dots, e_N\}$ to test the Continuum Coercivity Conjecture $C_{12} \succeq c_{12} I > 0$.
4. Sweep $M \in \{3, 4, 6, 8, 10, 12, 14, 16\}$ to observe the exact bound-state-to-continuum transition in $\lambda_{\min}(C_M)$.

