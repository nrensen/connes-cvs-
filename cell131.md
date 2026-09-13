# CELL 131 — ARCHIMEDEAN OFF-DIAGONAL CONTROL & PROJECTED NEGATIVE-POTENTIAL COMPRESSION ON $\Phi^\perp$

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.5)  
**Target Proposition:** Subspace Coercivity & Positivity on the Continuum Subspace:
$$\widehat{Q}_{\mathrm{even}} \equiv U_{\mathrm{cont}}^T Q_{\mathrm{even}} U_{\mathrm{cont}} = \widehat{\Omega} + \widehat{\Delta Q}_{\mathrm{arch}} - \widehat{\mathcal{K}}_{\mathrm{neg}} + \widehat{Q}_{\mathrm{pole}} \succeq E_{11}(N) I > 0$$
on the codimension-11 continuum spectral subspace $\Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}}) \subset \mathcal{H}_{\mathrm{even}}(N)$, with the asymptotic behavior $\liminf_{N \to \infty} E_{11}(N) \stackrel{?}{>} 0$ formulated as an open question.  
**Verification / Falsification Criteria:**  
1. **Probe 1 (Archimedean Off-Diagonal Control):** Measure the compressed operator norm $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} = \|U_{\mathrm{cont}}^T [Q_{\mathrm{arch}}^{\mathrm{even}} - \operatorname{diag}(h_+(a_m))] U_{\mathrm{cont}}\|_{\mathrm{op}}$ across $N \in \{16, 20, 24, 28, 32, 40, 48, 64\}$. If $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} \gg c_0 \approx 0.1567$, off-diagonal divided differences dominate the diagonal multiplier and the scalar multiplier cannot serve as the autonomous lower bound.
2. **Probe 2 (Projected Negative Potential Control):** Measure $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) = \lambda_{\max}(U_{\mathrm{cont}}^T (\widetilde{W} - \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}})$. If $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) < \lambda_{\min}(\widehat{\Omega})$, direct operator form domination $\widehat{\Omega} \succ \widehat{\mathcal{K}}_{\mathrm{neg}}$ holds on $\Phi^\perp$. If $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) > \lambda_{\min}(\widehat{\Omega})$, coercivity requires localized phase cancellation within the coupled Friedrichs form rather than separate form domination.  
**Companion Computational Script:** [`cell131.py`](file:///c:/data/github/connes-cvs-/cell131.py)  

---

## 1. Executive Summary & Epistemic Progression

### 1.1 The Context of the Investigation
The progression through Cells 127–130 successfully opened the black box of the prime operator:
1. **Cell 127:** Discovered the exact physical-space shifted autocorrelation representation:
   $$\langle v, Q_{\mathrm{prime}}^{\mathrm{even}} v \rangle = -\frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt.$$
2. **Cell 128:** Proved the Identical Vanishing Theorem $\mathcal{B}_q[v] \equiv 0$ for boundary mass mismatch via midpoint reflection symmetry $T_v(L - t) \equiv T_v(t)$.
3. **Cell 129:** Proved the **Two-Regime Multiplier Theorem**, establishing that the effective diagonal multiplier $\Omega(m) \equiv h_+(a_m) + 4 M(m)$ is unconditionally strictly positive across all modes:
   $$\inf_{m \ge 1} \Omega(m) \ge c_0 \equiv h_+(a_3) \approx +0.1567 > 0,$$
   completely overcoming the Diophantine incommensurability obstruction $\inf_{m \ge 1} M(m) = 0$.
4. **Cell 130:** Proved and certified to 50-digit precision the **Exact Component Decomposition Identity**:
   $$Q_{\mathrm{prime}}^{\mathrm{even}} \equiv -\widetilde{W} + \widetilde{\mathcal{D}}^{\mathrm{per}} + \Delta\widetilde{\mathcal{D}},$$
   with $\widetilde{W} \succeq 0$, $\widetilde{\mathcal{D}}^{\mathrm{per}} \succeq 0$, and $\Delta\widetilde{\mathcal{D}} \preceq 0$.

### 1.2 The Analytical Mission of Cell 131
With the exact tripartite decomposition certified, the central question of Gate 1 has been distilled down to its fundamental core:

$$\textbf{Can the negative operator } \mathcal{K}_{\mathrm{neg}} \equiv \widetilde{W} - \Delta\widetilde{\mathcal{D}} \textbf{ be controlled relative to the positive backbone on } \Phi^\perp?$$

To answer this, Cell 131 investigates two decoupled, precise analytical questions:
1. **Archimedean Off-Diagonal Control:** How large is the off-diagonal divided-difference remainder $\Delta Q_{\mathrm{arch}} = Q_{\mathrm{arch}}^{\mathrm{even}} - \operatorname{diag}(h_+(a_m))$ when projected onto $\Phi^\perp$?
2. **Projected Negative-Potential Control:** Does orthogonal projection onto the continuum spectral subspace $\Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}})$ suppress the negative potential $\mathcal{K}_{\mathrm{neg}} = \widetilde{W} - \Delta\widetilde{\mathcal{D}}$ sufficiently to permit positive form domination?

---

## 2. Mathematical Formulation of the Continuum Compression

### 2.1 The Codimension-11 Spectral Subspace $\Phi^\perp$
On the discrete even Galerkin space $\mathbb{R}^{N+1}$ ($N \ge 16$), the full Friedrichs matrix is:
$$Q_{\mathrm{even}} = Q_{\mathrm{arch}}^{\mathrm{even}} + Q_{\mathrm{prime}}^{\mathrm{even}} + Q_{\mathrm{pole}}^{\mathrm{even}}.$$
Let $Q_{\mathrm{even}} u_k = E_k u_k$ be its complete orthonormal eigendecomposition, sorted in ascending energy order:
$$E_0 \le E_1 \le \dots \le E_{10} < E_{11} \le \dots \le E_N.$$
For $c = 13$, the potential well carries $\bar{N}_{\mathrm{bound}} \approx 11$ bound states ($E_0, \dots, E_{10} \ll 1$), separated from the continuum by the finite-$N$ spectral gap:
$$g_{11}(N) \equiv E_{11}(N) - E_{10}(N) > 0.$$
*(Note on Epistemic Calibration: Systematic sweeps across $N \in [16, 64]$ demonstrate that $E_{11}(N)$ decreases from $1.976$ down to $0.006025$, while $E_{10} \downarrow 6.99 \times 10^{-6}$. The bound-state cluster and the continuum remain separated by three orders of magnitude ($g_{11} \approx 0.006018$, $E_{11}/E_{10} \approx 862$), but the absolute floor decreases with $N$, refuting the low-$N$ hypothesis of a fixed macroscopic continuum threshold at $\approx 0.58$.)*

We define the bound-state projector:
$$P_{\mathrm{bound}} \equiv \sum_{k=0}^{10} u_k u_k^T \in \mathbb{R}^{(N+1) \times (N+1)},$$
and the continuum spectral projector:
$$P_{\mathrm{cont}} \equiv I - P_{\mathrm{bound}} = \sum_{k=11}^N u_k u_k^T \in \mathbb{R}^{(N+1) \times (N+1)}.$$
The continuum subspace $\Phi^\perp = \operatorname{Ran}(P_{\mathrm{cont}})$ has dimension $q \equiv N - 10$.  
Let $U_{\mathrm{cont}} \in \mathbb{R}^{(N+1) \times q}$ be the orthonormal column matrix:
$$U_{\mathrm{cont}} \equiv \big[ u_{11}, u_{12}, \dots, u_N \big].$$
Because $U_{\mathrm{cont}}^T U_{\mathrm{cont}} = I_q$, $U_{\mathrm{cont}}$ is an exact isometry from $\mathbb{R}^q$ into $\mathbb{R}^{N+1}$.

### 2.2 Compressed Operators on $\Phi^\perp$
For any symmetric operator $A$ on $\mathcal{H}_{\mathrm{even}}(N)$, its compression to $\Phi^\perp$ is the $q \times q$ symmetric matrix:
$$\widehat{A} \equiv U_{\mathrm{cont}}^T A U_{\mathrm{cont}}.$$
By construction, the full Friedrichs operator compresses to a strictly diagonal matrix:
$$\widehat{Q}_{\mathrm{even}} = U_{\mathrm{cont}}^T Q_{\mathrm{even}} U_{\mathrm{cont}} = \operatorname{diag}\big(E_{11}, E_{12}, \dots, E_N\big).$$
Consequently:
$$\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) = E_{11}(N) > 0.$$
The continuum spectral floor is an exact, unconditional identity on $\Phi^\perp$ for each finite $N$. The central question is how its constituent physical pieces balance and how $E_{11}(N)$ behaves asymptotically as $N \to \infty$.

---

## 3. The Constituent Operators and the Exact Form Balance

From the Cell 130 exact decomposition, $Q_{\mathrm{even}}$ expands into:
$$Q_{\mathrm{even}} = Q_{\mathrm{arch}}^{\mathrm{even}} + \widetilde{\mathcal{D}}^{\mathrm{per}} - \widetilde{W} + \Delta\widetilde{\mathcal{D}} + Q_{\mathrm{pole}}^{\mathrm{even}}.$$
We define the constituent operators as follows:

### 3.1 The Diagonal Positive Multiplier Backbone ($\Omega_{\mathrm{diag}}$)
The continuous Fourier multiplier symbol induces the diagonal matrix:
$$D_{\mathrm{mult}} \equiv \operatorname{diag}\big(h_+(0), h_+(a_1), \dots, h_+(a_N)\big).$$
Adding the diagonal periodic translation-defect matrix $\widetilde{\mathcal{D}}^{\mathrm{per}} = \operatorname{diag}(0, 4M(1), \dots, 4M(N))$:
$$\Omega_{\mathrm{diag}} \equiv D_{\mathrm{mult}} + \widetilde{\mathcal{D}}^{\mathrm{per}} = \operatorname{diag}\big(h_+(0), \Omega(1), \Omega(2), \dots, \Omega(N)\big),$$
where $\Omega(m) = h_+(a_m) + 4M(m)$.  
- For $m \ge 1$: $\Omega(m) \ge c_0 \approx +0.1567 > 0$ by the Two-Regime Multiplier Theorem.
- For $m = 0$: $\Omega(0) = h_+(0) \approx -5.3722 < 0$.

#### Proposition 3.1 (Zero-Mode Isolation on $\Phi^\perp$)
For any unit vector $w \in \mathbb{R}^q$ ($\|w\|_2 = 1$), the corresponding physical state $v = U_{\mathrm{cont}} w \in \Phi^\perp$ satisfies:
$$\langle w, \widehat{\Omega} w \rangle = h_+(0) |v_0|^2 + \sum_{m=1}^N \Omega(m) v_m^2 \ge h_+(0) |v_0|^2 + c_0 (1 - |v_0|^2).$$
Because the bound states $\{u_0, \dots, u_{10}\}$ capture the spatial average of the potential well, the zero-mode leakage into $\Phi^\perp$:
$$\kappa_0(N) \equiv \|P_{\mathrm{cont}} e_0\|_2^2 = \sum_{k=11}^N |(u_k)_0|^2 = 1 - \sum_{k=0}^{10} |(u_k)_0|^2$$
governs the extent to which the negative $h_+(0)$ penetrates the continuum subspace. If $\kappa_0(N) \ll 1$, $\Phi^\perp$ is insulated from $h_+(0)$.

### 3.2 The Archimedean Off-Diagonal Defect ($\Delta Q_{\mathrm{arch}}$)
The divided-difference matrix $Q_{\mathrm{arch}}^{\mathrm{even}}$ deviates from its diagonal multiplier by:
$$\Delta Q_{\mathrm{arch}} \equiv Q_{\mathrm{arch}}^{\mathrm{even}} - D_{\mathrm{mult}}.$$
Its entries are:
$$(\Delta Q_{\mathrm{arch}})_{mn} = \begin{cases} (Q_{\mathrm{arch}}^{\mathrm{even}})_{mn}, & m \ne n, \\ \psi_{\mathrm{arch}}'(m) - h_+(a_m), & m = n. \end{cases}$$
Its compression onto $\Phi^\perp$ is:
$$\widehat{\Delta Q}_{\mathrm{arch}} \equiv U_{\mathrm{cont}}^T \Delta Q_{\mathrm{arch}} U_{\mathrm{cont}}.$$

### 3.3 The Total Negative Potential ($\mathcal{K}_{\mathrm{neg}}$)
The negative contributions to $Q_{\mathrm{prime}}^{\mathrm{even}}$ combine into:
$$\mathcal{K}_{\mathrm{neg}} \equiv \widetilde{W} - \Delta\widetilde{\mathcal{D}} \succeq 0.$$
Because $\widetilde{W} \succeq 0$ and $\Delta\widetilde{\mathcal{D}} \preceq 0$, $\mathcal{K}_{\mathrm{neg}}$ is strictly **positive semi-definite**. Its compression onto $\Phi^\perp$ is:
$$\widehat{\mathcal{K}}_{\mathrm{neg}} \equiv U_{\mathrm{cont}}^T (\widetilde{W} - \Delta\widetilde{\mathcal{D}}) U_{\mathrm{cont}} \succeq 0.$$

### 3.4 The Exact Subspace Decomposition Identity
Projecting all components onto $\Phi^\perp$ yields the exact algebraic balance:
$$\boxed{\widehat{Q}_{\mathrm{even}} \equiv \widehat{\Omega} + \widehat{\Delta Q}_{\mathrm{arch}} - \widehat{\mathcal{K}}_{\mathrm{neg}} + \widehat{Q}_{\mathrm{pole}},}$$
where $\widehat{Q}_{\mathrm{pole}} \equiv U_{\mathrm{cont}}^T Q_{\mathrm{pole}}^{\mathrm{even}} U_{\mathrm{cont}}$.  
Taking the trace:
$$\operatorname{Tr}(\widehat{Q}_{\mathrm{even}}) \equiv \operatorname{Tr}(\widehat{\Omega}) + \operatorname{Tr}(\widehat{\Delta Q}_{\mathrm{arch}}) - \operatorname{Tr}(\widehat{\mathcal{K}}_{\mathrm{neg}}) + \operatorname{Tr}(\widehat{Q}_{\mathrm{pole}}).$$

---

## 4. Analytical Hypotheses Under Test

### Hypothesis H-OffDiag (Archimedean Off-Diagonal Decoupling)
The off-diagonal divided-difference coupling decays like $\mathcal{O}(|m - n|^{-1})$ away from the diagonal. On the high-frequency subspace $\Phi^\perp$, destructive phase interference suppresses the off-diagonal operator norm:
$$\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} < \delta_{\mathrm{arch}} < \infty.$$
*Test:* Measure $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}}$ across discrete dimensions $N \in [16, 64]$ and compare with the diagonal floor $c_0 \approx 0.1567$.

### Hypothesis H-Supp (Negative-Potential Subspace Suppression)
On the unprojected Galerkin space, the step potential has maximum eigenvalue $\lambda_{\max}(\widetilde{W}) \approx 9.945$, which is large.  
However, the 11 bound states $u_0, \dots, u_{10}$ localize strongly inside the potential well. Orthogonal projection onto $\Phi^\perp$ forces test functions to be orthogonal to the potential well eigenstates:
$$\int_0^L T_v(t) T_{u_k}(t) \, dt = 0 \qquad (k = 0, \dots, 10).$$
*Hypothesis:* The compressed operator $\widehat{\mathcal{K}}_{\mathrm{neg}}$ satisfies:
$$\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}}) \ll \lambda_{\max}(\mathcal{K}_{\mathrm{neg}}) \approx 9.945.$$
*Quantitative Metric:* The suppression factor:
$$\mathcal{S}_{\mathrm{supp}}(N) \equiv \frac{\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}})}{\lambda_{\max}(\mathcal{K}_{\mathrm{neg}})} \ll 1.$$

### Hypothesis H-Dom (Coupled Backbone Form Domination)
Define the coupled competition operator:
$$\widehat{\mathcal{Q}}_{\mathrm{comp}} \equiv \widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}.$$
*Hypothesis:* The positive translation lift $4M(m)$ and Archimedean growth $h_+(a_m)$ dominate the projected negative potential:
$$\lambda_{\min}\big(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}\big) > 0.$$

---

## 5. Pre-Flight Computational Protocol (`cell131.py`)

The companion Python script [`cell131.py`](file:///c:/data/github/connes-cvs-/cell131.py) implements this audit under strict numerical standards:

1. **Parameters & Precision:**
   - Precision: `mp.mp.dps = 50`
   - Primary Cutoff: $c = 13$, $L = \log 13 \approx 2.56494935746$
   - Parameter $T$: $T = 600$
   - Dimension Grid: $N \in [16, 20, 24, 28, 32, 40, 48, 64]$
2. **Persistent Matrix Retrieval:**
   - Uses `get_galerkin_matrix` from [`cell.py`](file:///c:/data/github/connes-cvs-/cell.py) with persistent JSON caching to assemble $Q_{\mathrm{even}}$ without redundant quadrature evaluations.
3. **Subspace Isometry Construction:**
   - Computes eigendecomposition $Q_{\mathrm{even}} = V \operatorname{diag}(E) V^T$.
   - Sorts eigenvalues via `symmetric_eigenvalues` to guarantee exact ordering.
   - Extracts $U_{\mathrm{cont}} = V[:, 11:N+1]$.
4. **Audited Quantities per Dimension $N$:**
   - Continuum floor check: $\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) \stackrel{?}{=} E_{11}$.
   - Zero-mode leakage: $\kappa_0(N) = \|P_{\mathrm{cont}} e_0\|^2$.
   - Backbone compression: $\lambda_{\min}(\widehat{\Omega})$, $\lambda_{\max}(\widehat{\Omega})$, $\|\widehat{\Omega}\|_F$.
   - Archimedean off-diagonal defect: $\lambda_{\min}(\widehat{\Delta Q}_{\mathrm{arch}})$, $\lambda_{\max}(\widehat{\Delta Q}_{\mathrm{arch}})$, $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}}$.
   - Negative potential compression: $\lambda_{\min}(\widehat{\mathcal{K}}_{\mathrm{neg}})$, $\lambda_{\max}(\widehat{\mathcal{K}}_{\mathrm{neg}})$, $\|\widehat{\mathcal{K}}_{\mathrm{neg}}\|_{\mathrm{op}}$, and suppression ratio $\mathcal{S}_{\mathrm{supp}}$.
   - Coupled competition spectrum: $\lambda_{\min}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}})$, $\lambda_{\max}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}})$.
   - Pole contribution: $\lambda_{\min}(\widehat{Q}_{\mathrm{pole}})$, $\lambda_{\max}(\widehat{Q}_{\mathrm{pole}})$.
   - Trace identity residual: $|\operatorname{Tr}(\widehat{Q}_{\mathrm{even}}) - [\operatorname{Tr}(\widehat{\Omega}) + \operatorname{Tr}(\widehat{\Delta Q}_{\mathrm{arch}}) - \operatorname{Tr}(\widehat{\mathcal{K}}_{\mathrm{neg}}) + \operatorname{Tr}(\widehat{Q}_{\mathrm{pole}})]|$.

---

## 6. Audited Computational Results (`cell131.out` at 50 dps)

The high-precision computational suite [`cell131.py`](file:///c:/data/github/connes-cvs-/cell131.py) was executed to 50 decimal digits across $N \in [16, 64]$ with $c = 13$, $L = \log 13 \approx 2.56494935746$, and $T = 600$. The raw certified outputs are summarized below:

### Table 1: Continuum Subspace Dimension, Bound-State Gap & Zero-Mode Leakage
| $N$ | dim | $q$ | $E_{10}$ | $E_{11}$ (Floor) | Gap $g_{11}$ | $\kappa_0(e_0)$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 16 | 17 | 6 | $1.04298$ | $1.97648$ | $0.933505$ | $0.015361$ |
| 20 | 21 | 10 | $0.128165$ | $1.21128$ | $1.08312$ | $0.0221508$ |
| 24 | 25 | 14 | $0.00535019$ | $0.411765$ | $0.406415$ | $0.0291847$ |
| 28 | 29 | 18 | $0.000747055$ | $0.133716$ | $0.132969$ | $0.0332437$ |
| 32 | 33 | 22 | $0.000137141$ | $0.0415847$ | $0.0414475$ | $0.0364073$ |
| 40 | 41 | 30 | $1.68603 \times 10^{-5}$ | $0.0124547$ | $0.0124378$ | $0.0393668$ |
| 48 | 49 | 38 | $8.10097 \times 10^{-6}$ | $0.00924271$ | $0.00923461$ | $0.0398835$ |
| 64 | 65 | 54 | $6.98939 \times 10^{-6}$ | $0.00602464$ | $0.00601765$ | $0.0399547$ |

*Verification:* $\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) \equiv E_{11}$ exact to 50 decimal digits across all $N$.

### Table 2: Archimedean Off-Diagonal Defect $\Delta Q_{\mathrm{arch}}$ on Full Space vs $\Phi^\perp$
| $N$ | $\|\Delta Q\|_{\mathrm{full}}$ | $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}}$ | $[\lambda_{\min}, \lambda_{\max}]$ on $\Phi^\perp$ |
| :---: | :---: | :---: | :---: |
| 16 | $2.57983$ | $0.309411$ | $[1.6594 \times 10^{-9}, 0.30941]$ |
| 20 | $2.58090$ | $0.503063$ | $[-1.0122 \times 10^{-7}, 0.50306]$ |
| 24 | $2.58160$ | $0.646195$ | $[-2.3812 \times 10^{-6}, 0.64619]$ |
| 28 | $2.58208$ | $0.699990$ | $[-1.4037 \times 10^{-5}, 0.69999]$ |
| 32 | $2.58241$ | $0.734244$ | $[-6.8591 \times 10^{-5}, 0.73424]$ |
| 40 | $2.58282$ | $0.765904$ | $[-0.0006209, 0.76590]$ |
| 48 | $2.58304$ | $0.774924$ | $[-0.0024522, 0.77492]$ |
| 64 | $2.58323$ | $0.779593$ | $[-0.0136630, 0.77959]$ |

### Table 3: Negative Potential $\mathcal{K}_{\mathrm{neg}} = \widetilde{W} - \Delta\widetilde{\mathcal{D}}$ on Full Space vs $\Phi^\perp$
| $N$ | $\lambda_{\max}(\mathrm{full})$ | $\lambda_{\max}(\Phi^\perp)$ | $\lambda_{\min}(\Phi^\perp)$ | Suppression Ratio $\mathcal{S}_{\mathrm{supp}}$ |
| :---: | :---: | :---: | :---: | :---: |
| 16 | $11.7582$ | $10.3544$ | $6.39211$ | $0.880610$ |
| 20 | $11.7638$ | $10.3617$ | $5.93294$ | $0.880808$ |
| 24 | $12.3829$ | $12.3755$ | $5.88162$ | $0.999404$ |
| 28 | $12.3903$ | $12.3798$ | $5.89335$ | $0.999154$ |
| 32 | $12.3906$ | $12.3772$ | $5.89636$ | $0.998918$ |
| 40 | $12.7732$ | $12.7730$ | $5.89139$ | $0.999989$ |
| 48 | $12.8125$ | $12.8123$ | $5.21223$ | $0.999991$ |
| 64 | $12.8165$ | $12.8164$ | $5.15710$ | $0.999991$ |

### Table 4: Coupled Backbone $\widehat{\Omega}$, Competition Operator & Net Spectrum
| $N$ | $\lambda_{\min}(\widehat{\Omega})$ | $\lambda_{\min}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}})$ | $\lambda_{\min}(\widehat{Q}_{\mathrm{even}})$ | $|\Delta \operatorname{Tr}|$ Residual |
| :---: | :---: | :---: | :---: | :---: |
| 16 | $8.81772$ | $+1.725600$ | $1.97648$ | $3.635 \times 10^{-49}$ |
| 20 | $7.38634$ | $+0.878704$ | $1.21128$ | $4.704 \times 10^{-49}$ |
| 24 | $7.22742$ | $-0.008071$ | $0.411765$ | $9.408 \times 10^{-49}$ |
| 28 | $7.22795$ | $-0.284453$ | $0.133716$ | $5.132 \times 10^{-49}$ |
| 32 | $7.23745$ | $-0.386627$ | $0.041585$ | $3.421 \times 10^{-49}$ |
| 40 | $7.19916$ | $-0.457559$ | $0.012455$ | $3.421 \times 10^{-49}$ |
| 48 | $7.19915$ | $-0.478034$ | $0.009243$ | $2.737 \times 10^{-48}$ |
| 64 | $7.19915$ | $-0.486979$ | $0.006025$ | $1.026 \times 10^{-48}$ |

---

## 7. Epistemic Verdicts & Falsification Summary

The numerical evidence delivered by Cell 131 is exceptionally decisive. It resolves the three analytical hypotheses formulated in Section 4:

### 7.1 Hypothesis H-OffDiag: FALSIFIED
- The off-diagonal Archimedean operator norm does not decay to zero on $\Phi^\perp$. Instead, it increases with dimension and plateaus near $\|\widehat{\Delta Q}_{\mathrm{arch}}\|_{\mathrm{op}} \approx 0.780$ at $N = 64$.
- Because $0.780 \gg c_0 = 0.1567$, off-diagonal divided-difference coupling is an $\mathcal{O}(1)$ operator on the continuum subspace.
- **Epistemic Conclusion:** The diagonal multiplier $\operatorname{diag}(h_+(a_m) + 4M(m))$ cannot serve as an autonomous lower bound for the continuum operator; the off-diagonal terms cannot be discarded as minor perturbations.

### 7.2 Hypothesis H-Supp: DECISIVELY FALSIFIED
- Projection onto the continuum subspace $\Phi^\perp$ produces essentially zero suppression of the negative potential:
  $$\mathcal{S}_{\mathrm{supp}}(64) = \frac{12.816423}{12.816543} = 0.99999065 \approx 1.$$
- **Epistemic Conclusion:** The heuristic intuition that "bound states live in the well, so orthogonality to bound states shields continuum states from the well" is mathematically false in this discrete Galerkin setting. Coordinate localization and spectral projection cannot be identified.

### 7.3 Hypothesis H-Dom: FALSIFIED
- The separate competition operator $\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}$ fails to be positive definite for $N \ge 24$, falling to $\lambda_{\min} = -0.486979$ at $N = 64$.
- **Epistemic Conclusion:** The scalar backbone $\widehat{\Omega}$ does not autonomously dominate the negative potential $\widehat{\mathcal{K}}_{\mathrm{neg}}$ on $\Phi^\perp$.

### 7.4 The Emergent Reality: Coupled Cancellation Mechanism STRONGLY SUPPORTED
Despite the failure of all three separate domination hypotheses:
$$\boxed{\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) = E_{11}(N) = 0.00602464 > 0 \quad \text{strictly at } N = 64.}$$
The constituent operators balance as:
$$\widehat{Q}_{\mathrm{even}} = \underbrace{(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}})}_{\min = -0.48698} + \underbrace{\widehat{\Delta Q}_{\mathrm{arch}}}_{\|\cdot\| = 0.77959} + \underbrace{\widehat{Q}_{\mathrm{pole}}}_{\max = 0.29714}.$$
The full operator is lifted into strict positivity by the cooperative coupling between the Archimedean off-diagonal operator, the pole projector, and the competition Hamiltonian. Positivity is an **emergent property of the coupled Friedrichs form**, not a product of separate term-by-term domination.

### 7.5 Subspace Geometry of the Zero Mode
While approximately $4\%$ of the zero-mode direction penetrates $\Phi^\perp$ ($\kappa_0(64) \approx 0.03995$), the compressed backbone remains strictly bounded below:
$$\lambda_{\min}(\widehat{\Omega}) = 7.19915 \gg 0.$$
A crude unconstrained vector with $4\%$ zero-mode mass would yield $h_+(0) \kappa_0 + c_0(1 - \kappa_0) \approx -0.06 < 0$. The fact that the actual compressed operator is $+7.199$ proves that the geometry of the spectral subspace $\Phi^\perp$ strongly couples zero-mode leakage to high-frequency components, preventing alignment with the negative direction.

---

## 8. Strategic Pivot: Formulation of Cell 132

Having decisively falsified independent operator domination, the research programme must pivot from asking *"Does $\Omega$ dominate $\mathcal{K}_{\mathrm{neg}}$?"* to analyzing the **fine geometric structure of the coupled cancellation**:

1. **Eigenvector Dissection of the Competition Minimum:**
   Let $w_{\mathrm{bad}} \in \mathbb{R}^q$ be the normalized eigenvector corresponding to the negative eigenvalue $\lambda_{\min}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}) \approx -0.487$.
2. **Rayleigh Quotient Budget on the Negative Direction:**
   Evaluate the exact quadratic form components along $w_{\mathrm{bad}}$:
   $$R_{\mathrm{comp}}(w_{\mathrm{bad}}) \equiv \langle w_{\mathrm{bad}}, (\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}) w_{\mathrm{bad}} \rangle \approx -0.487,$$
   $$R_{\mathrm{arch}}(w_{\mathrm{bad}}) \equiv \langle w_{\mathrm{bad}}, \widehat{\Delta Q}_{\mathrm{arch}} w_{\mathrm{bad}} \rangle,$$
   $$R_{\mathrm{pole}}(w_{\mathrm{bad}}) \equiv \langle w_{\mathrm{bad}}, \widehat{Q}_{\mathrm{pole}} w_{\mathrm{bad}} \rangle.$$
   If $R_{\mathrm{arch}}(w_{\mathrm{bad}}) + R_{\mathrm{pole}}(w_{\mathrm{bad}}) \approx +0.493 + \epsilon$, it proves that the Archimedean off-diagonal operator acts as a **targeted restoring force** precisely aligned with the vulnerable directions of the negative potential well.
3. **Full Operator Norm Residual Audit:**
   Audit the full matrix residual:
   $$\|\widehat{Q}_{\mathrm{even}} - (\widehat{\Omega} + \widehat{\Delta Q}_{\mathrm{arch}} - \widehat{\mathcal{K}}_{\mathrm{neg}} + \widehat{Q}_{\mathrm{pole}})\|_{\mathrm{op}} < 10^{-45}$$
   to certify entrywise operator balance beyond the scalar trace identity.

---

## References
- [`cell130.md`](file:///c:/data/github/connes-cvs-/cell130.md) — Exact Component Decomposition of $Q_{\mathrm{prime}}^{\mathrm{even}}$
- [`cell129.md`](file:///c:/data/github/connes-cvs-/cell129.md) — Two-Regime Multiplier Theorem & Monotonicity of $h_+(r)$
- [`cell128.md`](file:///c:/data/github/connes-cvs-/cell128.md) — Identical Vanishing of Boundary Mismatch & Step Potential
- [`cell127.md`](file:///c:/data/github/connes-cvs-/cell127.md) — Shifted Autocorrelation Representation of $Q_{\mathrm{prime}}$
- [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md) — Core Submatrix Spectra & Cauchy Interlacing Bound
- [`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md) — Strategic Roadmap (Gate 1 Pipeline)

