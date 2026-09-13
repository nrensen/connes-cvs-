# CELL 137 — THE PROJECTED STEP-POTENTIAL OPERATOR, FULL 11-DIMENSIONAL CONSTANT-MODE ENCLOSURE, AND CONSTRAINED OPERATOR SUPREMUM

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.6 / Variational Lower Bound & Projected Step Potential)  
**Target Proposition:** Operator-Theoretic Enclosures on the Continuum Constraint Subspace $\mathcal{B}_\infty^\perp$:
1. **Theorem 137.1 (Subspace-Invariant Constant-Mode Enclosure):**
   For the full 11-dimensional bound-state space $\mathcal{B}_{11} = \operatorname{span}\{u_0, \dots, u_{10}\}$ of $Q_{\mathrm{even}}^{(N)}$, the maximal constant-mode amplitude on $\mathcal{B}_{11}^\perp$ is an invariant of the orthogonal projection:
   $$\sup_{\substack{T \in \mathcal{B}_{11}^\perp \\ \|T\|_{L^2} = 1}} |\langle T, e_0 \rangle|^2 = \|P_{\mathcal{B}_{11}^\perp} e_0\|^2 = 1 - \langle e_0, P_{11} e_0 \rangle = 1 - \sum_{k=0}^{10} |(u_k)_0|^2 \equiv \varepsilon_0(N).$$
   This eliminates sensitivity to individual eigenvector branch-crossings (such as $u_0 \leftrightarrow u_1$ near $N \in [48, 64]$) and establishes a rigorous zero-frequency floor:
   $$\mathcal{A}[T] + \mathcal{D}^{\mathrm{per}}[T] \ge \varepsilon_0(N) h_+(0) + (1 - \varepsilon_0(N)) \inf_{m \ge 1} \Omega(m) \qquad \forall T \in \mathcal{B}_{11}^\perp, \ \|T\|_{L^2} = 1.$$
2. **Theorem 137.2 (The Projected Step-Potential Operator & Supremum):**
   Define the projected step-potential operator on $\mathcal{B}_{11}^\perp$:
   $$W_{\perp \mathcal{B}} \equiv P_{\mathcal{B}_{11}^\perp} \widetilde{W} P_{\mathcal{B}_{11}^\perp}, \qquad \|W_{\perp \mathcal{B}}\|_{\mathrm{op}} = \sup_{\substack{T \in \mathcal{B}_{11}^\perp \\ \|T\|_{L^2} = 1}} \langle T, W T \rangle = \lambda_{\max}(U_{\mathrm{cont}}^T \widetilde{W} U_{\mathrm{cont}}).$$
   Evaluating $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$ resolves whether the observed $37.31\%$ well harvest:
   $$\frac{R_W(v_{\mathrm{bad}})}{W(L)} = \frac{3.709783}{9.943769} \approx 37.31\%$$
   is a universal geometric constraint enforced by $\mathcal{B}_{11}^\perp$ ($\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} \approx 3.71$), or an energy-minimization feature of the competition ground state $v_{\mathrm{bad}}$ balancing potential depth against kinetic cost.
3. **Proposition 137.3 (Operator-Splitting Lower Bound vs Coupled Ground State):**
   By Weyl's inequality, the competition deficit satisfies:
   $$\mu_0 \ge \lambda_{\min}\left( P_{\mathcal{B}_{11}^\perp} (\widehat{\Omega} + \Delta\widetilde{\mathcal{D}}) P_{\mathcal{B}_{11}^\perp} \right) - \|W_{\perp \mathcal{B}}\|_{\mathrm{op}}.$$
   Quantify the gap between this decoupled operator-splitting bound and the true coupled eigenvalue $\mu_0 \approx -0.4870$.

**Verification / Falsification Criteria:**
1. **Constant-Mode Leakage Floor:** Verify that $\varepsilon_0(N) = 1 - \sum_{k=0}^{10} |(u_k)_0|^2$ stabilizes monotonically as $N \to 64$.
2. **Step Potential Compression Audit:** Compute $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$ and determine whether $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} < \|\widetilde{W}\|_{\mathrm{op}} \le W(L)$.
3. **Minimizer Harvest Gap:** Quantify $\Delta_W \equiv \|W_{\perp \mathcal{B}}\|_{\mathrm{op}} - R_W(v_{\mathrm{bad}})$ to determine whether the $37\%$ phenomenon is a subspace constraint or an energy-driven feature.
4. **Finite-$N$ Margin:** Re-verify that the coupled competition ground state satisfies $\mu_0(N) > -1/2$ with margin $\Delta_{\mathrm{margin}} \approx +0.0130 > 0$.

**Companion Computational Script:** [`cell137.py`](file:///c:/data/github/connes-cvs-/cell137.py)  

---

## 1. Executive Context & Epistemic Motivation

In Cell 136, the operator-theoretic foundation of the continuum competition functional was rigorously established:
- $\mathcal{D}$ is an unconditionally bounded quadratic form on $L^2([0, L])$ with $\|\mathcal{D}\|_{\mathrm{op}} \le 2 W(L) \approx 19.8875$, meaning $\mathcal{Q}(\mathcal{D}) = L^2([0, L])$ (Theorem 136.1).
- $\mathcal{W}$ is an $L^\infty$ multiplication form with $\|\mathcal{W}\|_{\mathrm{op}} \le W(L) \approx 9.9438$ (Theorem 136.3).
- The exact quadratic form domain is the logarithmic Sobolev space $\mathcal{Q}(\mathcal{E}) = H^{\log}([0, L])$, semibounded from below by $(h_+(0) - W(L)) \approx -15.316$, uniquely generating a self-adjoint Friedrichs Hamiltonian $H_{\mathrm{comp}}$ with form domain $\mathcal{Q}(H_{\mathrm{comp}}) = H^{\log}([0, L])$ (Theorem 136.2).
- The bound-state projector sequence stabilizes rapidly: $\|P_{11}^{(48)} - P_{11}^{(64)}\|_{\mathrm{op}} = 1.47 \times 10^{-3}$.

However, as highlighted in the reviewer critique, a serious logical gap was exposed in the preliminary lower bound argument:
> **The fatal gap is the claimed $38\%$ bound on $\mathcal{W}[T]$.**  
> We observed $\mathcal{W}[v_{\mathrm{bad}}^{(64)}] / W(L) \approx 37.31\%$ for one specific numerical minimizer. This is not a universal inequality over all $T \in \mathcal{B}_\infty^\perp \cap H^{\log}$. Riemann–Lebesgue says Fourier coefficients of a fixed $L^1$ function decay; it does not prevent arbitrary unit vectors in an infinite-dimensional orthogonal complement from concentrating mass on $[\log 2, L]$.

To turn the empirical finite-$N$ margin ($\mu_0^{(64)} + 0.50 \approx +0.0130$) into an authentic operator-theoretic theorem, we must abandon heuristic assumptions about arbitrary unit vectors and investigate the actual operators on $\mathcal{B}_\infty^\perp$.

Cell 137 attacks two concrete operator-theoretic objects:
1. The **exact constant-mode leakage** $\|P_{\mathcal{B}_{11}^\perp} e_0\|^2$ under the full 11-dimensional bound-state projector.
2. The **projected step-potential operator** $W_{\perp \mathcal{B}} \equiv P_{\mathcal{B}_{11}^\perp} \widetilde{W} P_{\mathcal{B}_{11}^\perp}$ and its operator norm $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$.

---

## 2. The Full 11-Dimensional Constant-Mode Enclosure

Let $\mathcal{H}_N = \mathbb{R}^{N+1}$ denote the even Galerkin truncation space, with orthonormal basis $e_m$ ($m = 0, \dots, N$), where $e_0 = (1, 0, \dots, 0)^T$ represents the flat spatial zero mode $T_{e_0}(t) \equiv 1$.

Let $u_0, u_1, \dots, u_{10}$ be the lowest 11 orthonormal eigenvectors of the full Friedrichs operator $Q_{\mathrm{even}}^{(N)}$, spanning the bound-state subspace:
$$\mathcal{B}_{11}^{(N)} \equiv \operatorname{span}\{u_0, u_1, \dots, u_{10}\} \subset \mathcal{H}_N.$$

The orthogonal projector onto $\mathcal{B}_{11}^{(N)}$ is:
$$P_{11} = \sum_{k=0}^{10} u_k u_k^T.$$

### 2.1 Theorem 137.1 (Subspace-Invariant Constant-Mode Enclosure)
For any unit vector $T \in (\mathcal{B}_{11}^{(N)})^\perp$ with $\|T\|_2 = 1$, the zero-frequency Fourier coefficient $c_0 = \langle T, e_0 \rangle$ satisfies the exact upper bound:
$$\boxed{|c_0|^2 \le \varepsilon_0(N) \equiv 1 - \langle e_0, P_{11} e_0 \rangle = 1 - \sum_{k=0}^{10} |(u_k)_0|^2.}$$
Furthermore, the bound is sharp and achieved by the normalized vector $T_* = P_{\mathcal{B}_{11}^\perp} e_0 / \|P_{\mathcal{B}_{11}^\perp} e_0\|_2$:
$$\sup_{\substack{T \in (\mathcal{B}_{11}^{(N)})^\perp \\ \|T\|_2 = 1}} |\langle T, e_0 \rangle|^2 = \varepsilon_0(N).$$

*Proof.*  
Every unit vector $T \in \mathcal{H}_N$ with $T \perp \mathcal{B}_{11}^{(N)}$ belongs to the range of the complementary projector $P_\perp = I - P_{11}$.  
By Cauchy–Schwarz in $\mathcal{H}_N$:
$$|\langle T, e_0 \rangle|^2 = |\langle P_\perp T, e_0 \rangle|^2 = |\langle T, P_\perp e_0 \rangle|^2 \le \|T\|_2^2 \|P_\perp e_0\|_2^2 = \|P_\perp e_0\|_2^2.$$
Since $P_\perp$ is an orthogonal projector ($P_\perp^2 = P_\perp = P_\perp^T$):
$$\|P_\perp e_0\|_2^2 = \langle P_\perp e_0, P_\perp e_0 \rangle = \langle e_0, P_\perp e_0 \rangle = \langle e_0, (I - P_{11}) e_0 \rangle = 1 - \langle e_0, P_{11} e_0 \rangle.$$
Expanding $P_{11} = \sum_{k=0}^{10} u_k u_k^T$:
$$\langle e_0, P_{11} e_0 \rangle = \sum_{k=0}^{10} |\langle u_k, e_0 \rangle|^2 = \sum_{k=0}^{10} |(u_k)_0|^2.$$
Hence $\|P_\perp e_0\|_2^2 = 1 - \sum_{k=0}^{10} |(u_k)_0|^2 \equiv \varepsilon_0(N)$.  
Equality is achieved when $T = P_\perp e_0 / \|P_\perp e_0\|_2$, completing the proof. $\blacksquare$

### 2.2 Resolution of the Branch-Crossing Forensics
In `cell136.out`, the single-vector quantity $|u_0(0)|^2$ exhibited an apparent discontinuity between $N=40$ and $N=64$:
$$N=24: 29.90\%, \quad N=40: 29.27\%, \quad N=48: 21.33\%, \quad N=64: 1.56\%.$$
This occurred because near $N \in [48, 64]$, the nearly degenerate lowest bound states undergo a rotation/crossing, shifting zero-mode mass from the nominal index $k=0$ into excited indices $k \in \{1, \dots, 10\}$ (as documented in Cells 112a–114).

Theorem 137.1 demonstrates why analyzing $u_0$ in isolation is mathematically defective: the true invariant of the physical system is the **entire 11-dimensional bound-state subspace $\mathcal{B}_{11}$**.  
The total captured mass $\sum_{k=0}^{10} |(u_k)_0|^2$ is invariant under any orthogonal change of basis among the bound states $u_k \mapsto \sum_j O_{kj} u_j$.

### 2.3 The Rigorous Kinetic/Periodic Multiplier Floor
Recall the Two-Regime Multiplier Theorem (Cell 130):
$$\Omega(m) \equiv h_+\left(\frac{2\pi m}{L}\right) + 4 M(m) \ge 0.156708 \qquad \forall m \ge 1,$$
with only the $m=0$ mode negative: $\Omega(0) = h_+(0) \approx -5.372183$.

For any unit vector $T \in \mathcal{B}_{11}^\perp$:
$$\sum_{m=0}^N |c_m|^2 \Omega(m) = |c_0|^2 \Omega(0) + \sum_{m=1}^N |c_m|^2 \Omega(m) \ge |c_0|^2 h_+(0) + (1 - |c_0|^2) \inf_{m \ge 1} \Omega(m).$$
Because $h_+(0) < \inf_{m \ge 1} \Omega(m)$, the right-hand side is a strictly decreasing function of $|c_0|^2$. Substituting the maximum possible leakage $|c_0|^2 \le \varepsilon_0(N)$ from Theorem 137.1 yields the unconditional lower bound:
$$\boxed{\mathcal{A}[T] + \mathcal{D}^{\mathrm{per}}[T] \ge \varepsilon_0(N) h_+(0) + (1 - \varepsilon_0(N)) \times 0.156708 \qquad \forall T \in \mathcal{B}_{11}^\perp, \ \|T\|_2 = 1.}$$

---

## 3. The Projected Step-Potential Operator $W_{\perp \mathcal{B}}$

Let $U_{\mathrm{cont}} \in \mathbb{R}^{(N+1) \times (N-10)}$ be the matrix whose columns are the remaining orthonormal eigenvectors $u_{11}, \dots, u_N$ of $Q_{\mathrm{even}}^{(N)}$. The columns of $U_{\mathrm{cont}}$ form an orthonormal basis for $(\mathcal{B}_{11}^{(N)})^\perp$:
$$U_{\mathrm{cont}}^T U_{\mathrm{cont}} = I_{N-10}, \qquad U_{\mathrm{cont}} U_{\mathrm{cont}}^T = P_\perp = I - P_{11}.$$

### 3.1 Definition of the Projected Step Potential
The projected step-potential operator on $(\mathcal{B}_{11}^{(N)})^\perp$ is represented in this orthonormal basis by the $(N-10) \times (N-10)$ symmetric matrix:
$$\widehat{W}_{\perp} \equiv U_{\mathrm{cont}}^T \widetilde{W} U_{\mathrm{cont}}.$$

Its operator norm on $(\mathcal{B}_{11}^{(N)})^\perp$ is:
$$\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} \equiv \sup_{\substack{T \in \mathcal{B}_{11}^\perp \\ \|T\|_2 = 1}} \langle T, \widetilde{W} T \rangle = \lambda_{\max}(\widehat{W}_{\perp}).$$

### 3.2 The Fundamental Question: Universal Geometric Bound vs Minimizer Shape
The observed minimizer $v_{\mathrm{bad}} \in \mathcal{B}_{11}^\perp$ with $\|v_{\mathrm{bad}}\|_2 = 1$ satisfies $v_{\mathrm{bad}} = U_{\mathrm{cont}} w_{\mathrm{bad}}$ for some unit vector $w_{\mathrm{bad}} \in \mathbb{R}^{N-10}$. Its well harvest is the Rayleigh quotient:
$$R_W(v_{\mathrm{bad}}) = v_{\mathrm{bad}}^T \widetilde{W} v_{\mathrm{bad}} = w_{\mathrm{bad}}^T \widehat{W}_{\perp} w_{\mathrm{bad}}.$$

By the Rayleigh–Ritz theorem:
$$R_W(v_{\mathrm{bad}}) \le \lambda_{\max}(\widehat{W}_{\perp}) = \|W_{\perp \mathcal{B}}\|_{\mathrm{op}} \le \|\widetilde{W}\|_{\mathrm{op}} \le W(L) \approx 9.9438.$$

We now formulate the diagnostic dichotomy:
- **Case A (Universal Subspace Compression):**  
  If $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} \approx R_W(v_{\mathrm{bad}}) \approx 3.71$, then the $37.31\%$ well harvest is a **universal geometric property of $\mathcal{B}_{11}^\perp$**. Orthogonality to the 11 bound states alone prevents *any* wavepacket from harvesting more than $\approx 3.71$ from the step well.
- **Case B (Variational Selection):**  
  If $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} \gg R_W(v_{\mathrm{bad}})$ (e.g., $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} \approx 4.26 \approx \|\widetilde{W}\|_{\mathrm{op}}$), then the subspace $\mathcal{B}_{11}^\perp$ still contains states that harvest a larger fraction of the well. In this case, the $37.31\%$ harvest is **not** a universal geometric bound, but an **energy-selected shape**: the minimizer $v_{\mathrm{bad}}$ chooses to harvest only $3.71$ because deeper well harvesting would require steep gradients that incur excessive kinetic ($\mathcal{A}$) and translation ($\mathcal{D}^{\mathrm{true}}$) penalties.

Cell 137 computes $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$ to decide definitively between Case A and Case B.

---

## 4. Operator Splitting vs Coupled Eigenvalue Lower Bound

The discrete competition operator on the subspace $\mathcal{B}_{11}^\perp$ is:
$$\widehat{Q}_{\mathrm{comp}} = U_{\mathrm{cont}}^T \mathcal{Q}_{\mathrm{comp}} U_{\mathrm{cont}} = U_{\mathrm{cont}}^T (\Omega_{\mathrm{diag}} + \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}} - \widehat{W}_{\perp}.$$

Define the restoring operator on $\mathcal{B}_{11}^\perp$:
$$\widehat{K}_{\mathrm{rest}} \equiv U_{\mathrm{cont}}^T (\Omega_{\mathrm{diag}} + \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}} = U_{\mathrm{cont}}^T (D_{\mathrm{mult}} + D_{\mathrm{per}} + \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}} = U_{\mathrm{cont}}^T (D_{\mathrm{mult}} + D_{\mathrm{true}}) U_{\mathrm{cont}}.$$

Since $\widehat{Q}_{\mathrm{comp}} = \widehat{K}_{\mathrm{rest}} - \widehat{W}_{\perp}$, Weyl's lower bound inequality gives:
$$\mu_0 = \lambda_{\min}(\widehat{Q}_{\mathrm{comp}}) \ge \lambda_{\min}(\widehat{K}_{\mathrm{rest}}) - \lambda_{\max}(\widehat{W}_{\perp}) = \lambda_{\min}(\widehat{K}_{\mathrm{rest}}) - \|W_{\perp \mathcal{B}}\|_{\mathrm{op}}.$$

We define the **Operator-Splitting Lower Bound**:
$$\mu_0^{\mathrm{split}} \equiv \lambda_{\min}(\widehat{K}_{\mathrm{rest}}) - \|W_{\perp \mathcal{B}}\|_{\mathrm{op}}.$$

The gap $\Delta_{\mathrm{coupling}} \equiv \mu_0 - \mu_0^{\mathrm{split}} \ge 0$ represents the **coupling gain over the split Weyl bound**, measuring the extent to which the joint variational minimization of $K_{\mathrm{rest}} - W_\perp$ overcomes the pessimistic independent operator bound.

---

## 5. Mathematical Discovery: Case B (Variational Selection)

The high-precision execution of [`cell137.py`](file:///c:/data/github/connes-cvs-/cell137.py) certified two foundational facts:

1. **Full 11-Dimensional Constant-Mode Enclosure:**
   $$\varepsilon_0(64) = 4.00\%, \qquad \langle e_0, P_{11} e_0 \rangle = 96.00\%.$$
   The full 11-mode constraint subspace $\mathcal{B}_{11}$ captures $96.00\%$ of the zero mode, leaving at most $4.00\%$ leakage on $\mathcal{B}_{11}^\perp$. This yields the solid kinetic floor:
   $$\mathcal{A}[T] + \mathcal{D}^{\mathrm{per}}[T] \ge \varepsilon_0 h_+(0) + (1 - \varepsilon_0) \times 0.156708 \approx -0.064197 \qquad \forall T \in \mathcal{B}_{11}^\perp, \ \|T\| = 1.$$
2. **Projected Step-Potential Norm Decisively Evaluated:**
   $$\|P_{\mathcal{B}_{11}^\perp} \widetilde{W} P_{\mathcal{B}_{11}^\perp}\|_{\mathrm{op}} = 4.260495 \equiv \|\widetilde{W}\|_{\mathrm{op}}.$$
   The 11-mode constraint **does not materially compress the step potential at all** (identical to 6 decimal places).  
   *Conclusion (Decisive finite-$N$ evidence against Case A):* The observed $\approx 37.31\%$ well harvest ($R_W \approx 3.7098$) is **not** an artifact of the constraint projector $\mathcal{B}_{11}^\perp$. The subspace still contains admissible directions harvesting significantly more ($\|W_{\perp \mathcal{B}}\|_{\mathrm{op}} - R_W(v_{\mathrm{bad}}) = 0.550712$, a $13\%$ relative gap). It is an **energy-selection phenomenon (Case B)**: deeper well harvesting would require steep spatial gradients that incur excessive kinetic ($\mathcal{A}$) and translation ($\mathcal{D}^{\mathrm{true}}$) penalties.

---

## 6. Implementation Forensics & Regression Certification Protocol

### 6.1 Diagnosis of the Sign Error in `build_Delta_D_tilde_closed`
In the initial draft of `cell137.py`, a minus sign was omitted in the boundary truncation form:
```python
# BROKEN:
integrand_val = (mp.mpf("8") / L_PARAM) * sin_prod * J_val
entry += w_q * integrand_val

# REPAIRED (per Theorem 3.3 and Cell 135):
term = -(mp.mpf("8") / L_PARAM) * w_q * sin_prod * J_val
entry += term
```
Because $\Delta\mathcal{D}[v_{\mathrm{bad}}] < 0$, omitting this minus sign inverted $\Delta D \mapsto -\Delta D$, which corrupted $K_{\mathrm{neg}} = \widetilde{W} - \Delta D$ by $+2|\Delta D| \approx +15$ and shifted $Q_{\mathrm{comp}}$, causing $\mu_0$ to report $+6.998$ instead of the true $-0.48697922$.

### 6.2 Hard Regression Test Against Cell 135
To guarantee absolute consistency across the repository, `cell137.py` incorporates an automated regression check on $v_{\mathrm{bad}}^{(64)}$ against the certified Cell 135 baseline:

| Quantity | Evaluated $v_{\mathrm{bad}}^T (\cdot) v_{\mathrm{bad}}$ | Cell 135 Certified Reference | Status |
| :---: | :---: | :---: | :---: |
| $\mathcal{A}[T]$ | $1.173761$ | $1.173761$ | **MATCH** |
| $\mathcal{D}_{\mathrm{per}}[T]$ | $7.770307$ | $7.770307$ | **MATCH** |
| $\Delta\mathcal{D}[T]$ | $-5.721264$ | $-5.721264$ | **MATCH** |
| $\mathcal{D}_{\mathrm{true}}[T]$ | $2.049043$ | $2.049043$ | **MATCH** |
| $\mathcal{W}[T]$ | $3.709783$ | $3.709783$ | **MATCH** |
| $\mathcal{E}[T]$ ($\mu_0$) | $-0.48697922$ | $-0.48697922$ | **CERTIFIED** |

The algebraic identity $1.173761 - 3.709783 + 2.049043 = -0.48697922$ closes to $< 10^{-45}$, proving exact mathematical consistency with Cells 135–136.

---

## 7. Certified 50-DPS Numerical Results

Executed via external compute node (`cell137.out`, commit `28f8015`), runtime $185.04\text{ s}$ across $N \in [24, 64]$:

### Table 1: 11-Dimensional Zero-Mode Mass Enclosure & Kinetic Floor
$$\varepsilon_0(N) = 1 - \sum_{k=0}^{10} |(u_k)_0|^2 = \max_{T \perp \mathcal{B}_{11}} |\langle T, e_0 \rangle|^2$$

| $N$ | $|u_0(0)|^2$ | Captured Mass in $\mathcal{B}_{11}$ | Max Leakage $\varepsilon_0(N)$ | $|v_{\mathrm{bad}}(0)|^2$ | Kinetic Floor |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $29.90\%$ | $97.08\%$ | $2.92\%$ | $1.74\%$ | $-0.004651$ |
| 28 | $29.64\%$ | $96.68\%$ | $3.32\%$ | $2.05\%$ | $-0.027093$ |
| 32 | $29.48\%$ | $96.36\%$ | $3.64\%$ | $2.39\%$ | $-0.044584$ |
| 40 | $29.27\%$ | $96.06\%$ | $3.94\%$ | $2.76\%$ | $-0.060947$ |
| 48 | $21.33\%$ | $96.01\%$ | $3.99\%$ | $2.84\%$ | $-0.063803$ |
| **64** | **$1.56\%$** | **$96.00\%$** | **$4.00\%$** | **$2.86\%$** | **$-0.064197$** |

*Key finding:* Despite individual eigenvector rotation ($|u_0(0)|^2$ drops from $29.90\%$ to $1.56\%$), the total captured zero-mode mass remains invariant at $96.00\%$, bounding leakage to $\le 4.00\%$ and providing a rigorous kinetic floor of $\ge -0.0642$.

---

### Table 2: Projected Step-Potential Operator Norm
$$W_{\perp \mathcal{B}} = P_{\mathcal{B}_{11}^\perp} \widetilde{W} P_{\mathcal{B}_{11}^\perp}, \qquad \|W_{\perp \mathcal{B}}\|_{\mathrm{op}} = \lambda_{\max}(\widehat{W}_{\perp})$$

| $N$ | $\|\widetilde{W}\|_{\mathrm{op}}$ (Unprojected) | $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$ (Projected) | $\lambda_{\max}$ | Ratio to $W(L)$ | Compression |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $4.260495$ | $4.244045$ | $4.244045$ | $42.68\%$ | $99.61\%$ |
| 28 | $4.260495$ | $4.258124$ | $4.258124$ | $42.82\%$ | $99.94\%$ |
| 32 | $4.260495$ | $4.260232$ | $4.260232$ | $42.84\%$ | $99.99\%$ |
| 40 | $4.260495$ | $4.260493$ | $4.260493$ | $42.85\%$ | $100.00\%$ |
| 48 | $4.260495$ | $4.260495$ | $4.260495$ | $42.85\%$ | $100.00\%$ |
| **64** | **$4.260495$** | **$4.260495$** | **$4.260495$** | **$42.85\%$** | **$100.00\%$** |

---

### Table 3: Term-by-Term Energy Decomposition & Regression Audit

| $N$ | $\mathcal{A}[T]$ | $\mathcal{D}_{\mathrm{per}}[T]$ | $\Delta\mathcal{D}[T]$ | $\mathcal{D}_{\mathrm{true}}[T]$ | $\mathcal{W}[T]$ | $\mathcal{E}[T]$ ($\mu_0$) | Status |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $1.436458$ | $7.643216$ | $-5.377377$ | $2.265839$ | $3.710369$ | $-0.00807113$ | OK |
| 28 | $1.392251$ | $7.665071$ | $-5.667212$ | $1.997859$ | $3.674564$ | $-0.28445278$ | OK |
| 32 | $1.307048$ | $7.704757$ | $-5.709007$ | $1.995750$ | $3.689425$ | $-0.38662701$ | OK |
| 40 | $1.202064$ | $7.752913$ | $-5.712072$ | $2.040841$ | $3.700465$ | $-0.45755948$ | OK |
| 48 | $1.178727$ | $7.764826$ | $-5.714910$ | $2.049916$ | $3.706678$ | $-0.47803442$ | OK |
| **64** | **$1.173761$** | **$7.770307$** | **$-5.721264$** | **$2.049043$** | **$3.709783$** | **$-0.48697922$** | **CERTIFIED** |

---

### Table 4: Universal Supremum vs Variational Minimizer Harvest

| $N$ | $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$ | $R_W(v_{\mathrm{bad}})$ | Harvest $\%$ | Harvest Gap $\Delta_W$ | Dichotomy Verdict |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $4.244045$ | $3.710369$ | $37.31\%$ | $0.533676$ | Case B (Variational) |
| 28 | $4.258124$ | $3.674564$ | $36.95\%$ | $0.583561$ | Case B (Variational) |
| 32 | $4.260232$ | $3.689425$ | $37.10\%$ | $0.570807$ | Case B (Variational) |
| 40 | $4.260493$ | $3.700465$ | $37.21\%$ | $0.560029$ | Case B (Variational) |
| 48 | $4.260495$ | $3.706678$ | $37.28\%$ | $0.553818$ | Case B (Variational) |
| **64** | **$4.260495$** | **$3.709783$** | **$37.31\%$** | **$0.550712$** | **Case B (Variational)** |

---

### Table 5: Operator Splitting vs Coupled Competition Ground State
$$\mu_0 \ge \mu_0^{\mathrm{split}} \equiv \lambda_{\min}(K_{\mathrm{rest}}) - \|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$$

| $N$ | $\lambda_{\min}(K_{\mathrm{rest}})$ | $\|W_{\perp \mathcal{B}}\|_{\mathrm{op}}$ | $\mu_0^{\mathrm{split}}$ | $\mu_0$ (Coupled) | Coupling Gain $\Delta_{\mathrm{coupling}}$ | Margin to $-1/2$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 24 | $3.481342$ | $4.244045$ | $-0.762703$ | $-0.00807113$ | $+0.754632$ | $+0.49192887$ |
| 28 | $3.059115$ | $4.258124$ | $-1.199009$ | $-0.28445278$ | $+0.914556$ | $+0.21554722$ |
| 32 | $2.992230$ | $4.260232$ | $-1.268002$ | $-0.38662701$ | $+0.881375$ | $+0.11337299$ |
| 40 | $2.948333$ | $4.260493$ | $-1.312160$ | $-0.45755948$ | $+0.854601$ | $+0.04244052$ |
| 48 | $2.935956$ | $4.260495$ | $-1.324539$ | $-0.47803442$ | $+0.846505$ | $+0.02196558$ |
| **64** | **$2.931526$** | **$4.260495$** | **$-1.328969$** | **$-0.48697922$** | **$+0.841990$** | **$+0.01302078$** |

---

## 8. Epistemic Synthesis: The Three-Level Hierarchy

The certified results of Cell 137 establish a clear three-level conceptual hierarchy:

1. **Level 1 — Subspace Zero-Mode Enclosure:**
   $$\mathcal{A}[T] + \mathcal{D}^{\mathrm{per}}[T] \ge \varepsilon_0(64) h_+(0) + (1 - \varepsilon_0(64)) \times 0.156708 \ge -0.064197.$$
   The full 11-dimensional constraint eliminates $96.00\%$ of the dangerous negative Archimedean zero mode ($h_+(0) \approx -5.3722$), proving that low-frequency Archimedean leakage cannot derail positivity on $\mathcal{B}_{11}^\perp$.

2. **Level 2 — Decoupled Operator Splitting (Weyl Bound):**
   $$\mu_0 \ge \mu_0^{\mathrm{split}} = \lambda_{\min}(K_{\mathrm{rest}}) - \|W_{\perp \mathcal{B}}\|_{\mathrm{op}} = 2.931526 - 4.260495 = -1.328969.$$
   This proves that **separate coercivity is mathematically insufficient**. Treating the restoring operator $K_{\mathrm{rest}}$ and the step well $W_\perp$ independently yields a bound well below $-1/2$, because the unrestricted supremum of the step potential on $\mathcal{B}_{11}^\perp$ remains large ($4.2605$).

3. **Level 3 — Coupled Variational Minimization:**
   $$\mu_0 = -0.48697922 > -1/2 \qquad (\text{Finite-}N\text{ margin } +0.01302078).$$
   The true coupled operator outperforms the split Weyl bound by:
   $$\Delta_{\mathrm{coupling}} \equiv \mu_0 - \mu_0^{\mathrm{split}} \approx +0.841990.$$
   At $N=64$, the positive margin above $-1/2$ arises entirely from the coupled variational problem and is lost under independent operator splitting. Because the finite-$N$ margin decreases monotonically across tested dimensions ($+0.4919 \to +0.2155 \to +0.1134 \to +0.0424 \to +0.0220 \to +0.0130$), this is an empirical finite-$N$ structural fact, not a theorem of continuum positivity. Because $K_{\mathrm{rest}}$ and $W_\perp$ do not commute, any wavepacket attempting to harvest the full well depth ($4.2605$) incurs severe restoring penalties, forcing the true minimizer into a compromised state harvesting only $3.7098$.

This shifts the scientific frontier to **Cell 138**: investigating the exact variational decomposition, operator incompatibility, and the collective spectral geometry of $K_{\mathrm{rest}}$ and $W_\perp$.

