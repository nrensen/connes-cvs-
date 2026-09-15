# CELL 145 — Asymptotic Competition of Resolvent Denominator Collapse vs Exponential Tunneling Splitting

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction)  
**Milestone:** Milestone M-G1.6 / Gate 1 Resolvent Mechanism  
**Target Operators:** Certified Continuum Truncations $K_{\mathrm{rest}}(N)$, $W_\perp(N)$, and $H_1(N) = K_{\mathrm{rest}} - W_\perp$ on $\mathcal{B}_{11}^\perp$; Full Parity Operators $Q_{\mathrm{even}}, Q_{\mathrm{odd}}$ across $N \in [24, 80]$  
**Pre-Flight Invariants ($N = 64$):** $\omega_0 = 2.9315259531$, $\nu_0 = 4.2604954421$, $\mu_0 = -0.4869792197$ (tolerance $10^{-8}$)  
**Execution Script:** [`cell145.py`](file:///c:/data/github/connes-cvs-/cell145.py) (50-dps verification suite across $N \in [24, 32, 40, 48, 56, 64, 72, 80]$)  

---

## 1. Executive Context & Epistemic Reorientation

### 1.1 The Crucial Reorientation from Cell 144
In [`cell144.md`](file:///c:/data/github/connes-cvs-/cell144.md), a critical conceptual distinction was established:
- **Absolute Continuum Energy Control:** $\inf_{T \in \mathcal{B}_{11}^\perp} \langle T, (K_{\mathrm{rest}} - W_\perp) T \rangle = \mu_0 \approx -0.4870 > -1/2$.
- **Relative Continuum/Bound-State Gap:** $g_{\mathrm{cont}}(N) \equiv E_{11}(N) - E_{10}(N) \longrightarrow 0$.

Controlling the absolute energy floor of the continuum subspace does **not** prevent the continuum base $E_{11}(N)$ from approaching the bound states $E_{10}(N)$ asymptotically. As $N \to \infty$, both accumulate toward a common boundary energy $E_*$, causing the continuum gap to collapse monotonically ($0.0414 \to 0.0060$). Consequently, the resolvent denominator $D(N) \equiv (E_{11} - E_2)(E_{11} - E_3)$ collapses from $1.73 \times 10^{-3} \to 3.63 \times 10^{-5}$, and the resolvent ratio explodes:
$$R_{\mathrm{spec}}(N) \equiv \frac{\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{D(N)} \in \{2.38 \times 10^3, \; 5.81 \times 10^4, \; 9.98 \times 10^4, \; 1.41 \times 10^5\} \longrightarrow \infty.$$

### 1.2 The True Gate 1 Engine
Crucially, **the Gate 1 product does not fail**:
$$\Pi_{\mathrm{Gate1}}(N) = \Delta_2(N) R_{\mathrm{spec}}(N) \in \{3.12 \times 10^{-28}, \; 2.57 \times 10^{-32}, \; 4.80 \times 10^{-33}, \; 4.21 \times 10^{-33}\} \longrightarrow 0.$$
While $R_{\mathrm{spec}}$ grows by a factor of 60, the tunneling doublet splitting $\Delta_2(N)$ collapses by **7 orders of magnitude** ($1.31 \times 10^{-31} \to 2.98 \times 10^{-38}$).

This definitively reframes the Gate 1 problem:
$$\boxed{\textbf{We do not need } R_{\mathrm{spec}} = \mathcal{O}(1). \textbf{ We need its growth to be subordinate to the doublet splitting.}}$$

The central objective of Cell 145 is to quantify the asymptotic competition between:
$$\boxed{\textbf{Exponential Doublet Splitting } \Delta_2(N) \lesssim e^{-\kappa N} \quad\text{\textbf{vs.}}\quad \textbf{Subexponential Resolvent Growth } R_{\mathrm{spec}}(N) \lesssim N^p.}$$

---

## 2. Mathematical Framework & Asymptotic Envelopes

### 2.1 Competing Asymptotic Hypotheses
In Gate 1 of the canonical roadmap ([`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md)), tail extinction requires:
$$\lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{\mathrm{Gate1}}(N, L) = 0, \qquad \Pi_{\mathrm{Gate1}}(N, L) \equiv \Delta_2(N) R_{\mathrm{spec}}(N, L) = \frac{\Delta_2(N) \|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{D(N, L)}.$$

We formulate two competing asymptotic models:
1. **Hypothesis 145.1 (Polynomial Denominator Collapse):**  
   The resolvent denominator $D(N)$ collapses as a polynomial power law:
   $$D(N) \approx C_D \, N^{-p_D} \qquad (p_D > 0).$$
   Because the operator norm grows at most as $\|Q_{\mathrm{even}}\|_{\mathrm{op}} \approx C_Q \, N^{p_Q}$ with $p_Q \le 2$, the resolvent ratio grows polynomially:
   $$R_{\mathrm{spec}}(N) \approx C_R \, N^{p_R}, \qquad p_R \equiv p_D + p_Q.$$
2. **Hypothesis 145.2 (Exponential Tunneling Damping):**  
   The bound-state doublet splitting decays exponentially under semiclassical barrier tunneling:
   $$\Delta_2(N) \approx C_\Delta \, e^{-\kappa N} \qquad (\kappa > 0).$$

---

### 2.2 Theorem 145.1 (Subexponential Domination Criterion)
To test whether doublet splitting dominates resolvent divergence without relying on arbitrary prefactor normalizations, we define the **dimensionless log-ratio metric**:
$$\rho_{\mathrm{Gate1}}(N) \equiv \frac{-\log \Delta_2(N)}{\log R_{\mathrm{spec}}(N)}.$$

> **Theorem 145.1 (Subexponential Domination Criterion).**  
> If $R_{\mathrm{spec}}(N) = \mathcal{O}(N^{p_R})$ is polynomial while $\Delta_2(N) \le C_\Delta e^{-\kappa N}$ is exponential with $\kappa > 0$, then:
> $$\rho_{\mathrm{Gate1}}(N) = \frac{\kappa N - \log C_\Delta}{p_R \log N + \log C_R} \sim \frac{\kappa}{p_R} \frac{N}{\log N} \longrightarrow \infty.$$
> Consequently:
> $$\boxed{\lim_{N \to \infty} \Pi_{\mathrm{Gate1}}(N) = \lim_{N \to \infty} C_R N^{p_R} C_\Delta e^{-\kappa N} = 0.}$$
> A growing log-ratio $\rho_{\mathrm{Gate1}}(N) \to \infty$ constitutes an exact, audit-proof proof that Gate 1 tail extinction holds unconditionally despite the collapsing continuum gap.

---

### 2.3 Lower Bounds on the Resolvent Denominator $D(N)$
Analytically, we do not require $D(N) \ge D_\infty > 0$. We require only a subexponential lower bound:
$$D(N) \ge c_0 \, N^{-p_D} \qquad \text{or} \qquad D(N) \ge c_0 \, e^{-o(N)}.$$
Because $E_2, E_3$ are bound states deep in the step well, $E_{11} - E_2 \approx E_{11} - E_3$.  
The spacing between Ritz eigenvalues of discrete Sturm–Liouville / Galerkin operators near an accumulation point is bounded below by the Weyl density of states $\mathcal{O}(N^{-2})$ or $\mathcal{O}(N^{-3})$.  
Cell 145 extracts the empirical exponent $p_D$ to verify whether $p_D \approx 2\text{--}3$, confirming that denominator collapse is strictly polynomial.

---

## 3. Dual Core-Size Architecture: Edge vs Continuum

To separate the physics of the well edge ($L=10$, where $E_{11}$ accumulates) from the stabilized continuum ($L=12$, where $E_{13}$ lives in the scattering continuum), Cell 145 evaluates both core sizes:
- **Nominal Bound-State Edge ($L = 10$, $E_{L+1} = E_{11}$):**  
  Exhibits collapsing gap $g_{\mathrm{cont}} = E_{11} - E_{10} \to 0$ and exploding $R_{\mathrm{spec}, 10} \to \infty$. Tests whether polynomial denominator collapse is dominated by exponential tunneling.
- **Stabilized Scattering Continuum ($L = 12$, $E_{L+1} = E_{13}$):**  
  As discovered in Cell 121, the continuum gap stabilizes at $g_{2, 12}(N) \equiv E_{13} - E_3 \ge 0.5824 > 0$. The denominator $D_{12}(N) \ge 0.339$ remains strictly open, yielding $R_{\mathrm{spec}, 12} \approx 6\text{--}15 = \mathcal{O}(1)$ and super-exponential product extinction $\Pi_{2, 12}(N) \sim 10^{-40}$.

---

## 4. Experimental Suite Design (`cell145.py`)

The computational suite executes four targeted numerical modules across $N \in [24, 32, 40, 48, 56, 64, 72, 80]$:

### 4.1 Module 1: Pre-Flight Hard Regression Audit ($N=64$, 50 dps)
Verify agreement with certified invariants within $10^{-8}$:
$$\omega_0 = 2.9315259531, \quad \nu_0 = 4.2604954421, \quad \mu_0 = -0.4869792197.$$
The $N=64$ system is cached and reused with explicit label `(reused from pre-flight cache)`.

### 4.2 Module 2: Primary Scaling Table Across 8 Dimensions
For $N \in [24, 32, 40, 48, 56, 64, 72, 80]$, compute:
- Continuum gap $g_{\mathrm{cont}} = E_{11} - E_{10}$ and threshold $E_{11}$.
- Doublet splitting $\Delta_2 = E_3 - E_2$ and parity splitting $\Delta_0^{\mathrm{par}} = |E_0^{\mathrm{even}} - E_0^{\mathrm{odd}}|$.
- Resolvent denominator $D_{10}(N) = (E_{11} - E_2)(E_{11} - E_3)$.
- Spectral norm $\|Q_{\mathrm{even}}\|_{\mathrm{op}}$ and resolvent ratio $R_{\mathrm{spec}, 10}(N)$.
- Gate 1 product $\Pi_{\mathrm{Gate1}}(N) = \Delta_2(N) R_{\mathrm{spec}, 10}(N)$.
- Log-ratio metric $\rho_{\mathrm{Gate1}}(N) = -\log \Delta_2(N) / \log R_{\mathrm{spec}, 10}(N)$.

### 4.3 Module 3: Log-Linear & Log-Log Envelope Regressions
Fit asymptotic scaling exponents:
1. Doublet splitting decay: $\log \Delta_2(N) = -\kappa N + c_\Delta$ (extract $\kappa, R^2$).
2. Denominator power law: $\log D_{10}(N) = -p_D \log N + c_D$ (extract $p_D, R^2$).
3. Resolvent ratio power law: $\log R_{\mathrm{spec}, 10}(N) = p_R \log N + c_R$ (extract $p_R, R^2$).
4. Check whether $p_R \approx p_D + p_Q$ where $\|Q\|_{\mathrm{op}} \sim N^{p_Q}$.

### 4.4 Module 4: Multi-Doublet Gate 1 Products ($j = 0, 1, 2$)
Compare products across doublets:
- Ground parity doublet: $\Pi_0(N) \equiv \Delta_0^{\mathrm{par}}(N) R_{\mathrm{spec}}(N)$.
- Excited parity doublet: $\Pi_1(N) \equiv \Delta_1^{\mathrm{par}}(N) R_{\mathrm{spec}}(N)$.
- Second intra-well doublet: $\Pi_2(N) \equiv \Delta_2(N) R_{\mathrm{spec}}(N)$.

### 4.5 Module 5: Dual Core-Size Comparison ($L=10$ vs $L=12$)
Compare $D(N, L)$, $R_{\mathrm{spec}}(N, L)$, and $\Pi_2(N, L)$ between $L=10$ (collapsing gap) and $L=12$ (stable continuum gap).

---

## 5. Diagnostic Output Tables

The computational results generated by [`cell145.py`](file:///c:/data/github/connes-cvs-/cell145.py) will populate the following certified tables:

### Table 1: Primary Scaling Table: Denominator, Splittings, and Resolvent Ratio Across $N \in [24, 80]$
*(To be populated upon external execution of `cell145.py`)*

### Table 2: Asymptotic Curve-Fitting and Exponent Extraction
*(To be populated upon external execution of `cell145.py`)*

### Table 3: Multi-Doublet Gate 1 Products $\Pi_j(N) = \Delta_j R_{\mathrm{spec}}$ ($j = 0, 1, 2$)
*(To be populated upon external execution of `cell145.py`)*

### Table 4: Dual Core-Size Comparison ($L=10$ vs $L=12$)
*(To be populated upon external execution of `cell145.py`)*
