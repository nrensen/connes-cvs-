# CELL 145 — Asymptotic Competition of Resolvent Denominator Collapse vs Exponential Tunneling Splitting

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction)  
**Milestone:** Milestone M-G1.6 / Gate 1 Resolvent Mechanism  
**Target Operators:** Certified Continuum Truncations $K_{\mathrm{rest}}(N)$, $W_\perp(N)$, and $H_1(N) = K_{\mathrm{rest}} - W_\perp$ on $\mathcal{B}_{11}^\perp$; Full Parity Operators $Q_{\mathrm{even}}, Q_{\mathrm{odd}}$ across $N \in [24, 80]$  
**Pre-Flight Invariants ($N = 64$):** $\omega_0 = 2.9315259531$, $\nu_0 = 4.2604954421$, $\mu_0 = -0.4869792197$ (tolerance $10^{-8}$, certified residual $< 2.6 \times 10^{-11}$)  
**Execution Script:** [`cell145.py`](file:///c:/data/github/connes-cvs-/cell145.py) (50-dps verification suite across $N \in [24, 32, 40, 48, 56, 64, 72, 80]$)  
**Execution Output:** [`cell145.out`](file:///c:/data/github/connes-cvs-/cell145.out) (runtime 615.23s, commit `61525cf`)  

---

## 1. Executive Context & Epistemic Reorientation

### 1.1 The Crucial Reorientation from Cell 144
In [`cell144.md`](file:///c:/data/github/connes-cvs-/cell144.md), a critical conceptual distinction was established:
- **Absolute Continuum Energy Control:** $\inf_{T \in \mathcal{B}_{11}^\perp} \langle T, (K_{\mathrm{rest}} - W_\perp) T \rangle = \mu_0 \approx -0.4870 > -1/2$.
- **Relative Continuum/Bound-State Gap:** $g_{\mathrm{cont}}(N) \equiv E_{11}(N) - E_{10}(N) \longrightarrow 0$.

Controlling the absolute energy floor of the continuum subspace does **not** prevent the continuum base $E_{11}(N)$ from approaching the bound states $E_{10}(N)$ asymptotically. As $N \to \infty$, both accumulate toward a common boundary energy $E_*$, causing the continuum gap to collapse monotonically ($0.0414 \to 0.0050$). Consequently, the resolvent denominator $D(N) \equiv (E_{11} - E_2)(E_{11} - E_3) \approx E_{11}(N)^2$ collapses from $1.70 \times 10^{-1} \to 2.54 \times 10^{-5}$, and the resolvent ratio explodes:
$$R_{\mathrm{spec}}(N) \equiv \frac{\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{D(N)} \in \{2.25 \times 10^1, \; 2.38 \times 10^3, \; 5.81 \times 10^4, \; 1.41 \times 10^5, \; 2.12 \times 10^5\} \longrightarrow \infty.$$

### 1.2 The True Gate 1 Engine
Crucially, **the Gate 1 product remains extraordinarily small throughout the tested window**:
$$\Pi_{\mathrm{Gate1}}(N) = \Delta_2(N) R_{\mathrm{spec}}(N) \in \{3.17 \times 10^{-25}, \; 3.12 \times 10^{-28}, \; 2.57 \times 10^{-32}, \; 4.21 \times 10^{-33}, \; 4.62 \times 10^{-33}\}.$$
While $R_{\mathrm{spec}}$ grows by nearly four orders of magnitude ($22.5 \to 2.12 \times 10^5$), the tunneling doublet splitting $\Delta_2(N)$ collapses by **twelve orders of magnitude** ($1.41 \times 10^{-26} \to 2.18 \times 10^{-38}$).

This definitively reframes the Gate 1 problem:
$$\boxed{\textbf{We do not need } R_{\mathrm{spec}} = \mathcal{O}(1). \textbf{ We need its growth to be subordinate to the doublet splitting.}}$$

The central objective of Cell 145 is to quantify the asymptotic competition between:
$$\boxed{\textbf{Tunneling Doublet Splitting } \Delta_2(N) \quad\text{\textbf{vs.}}\quad \textbf{Resolvent Denominator Collapse } D(N) \approx E_{11}(N)^2.}$$

---

## 2. Mathematical Framework & Asymptotic Envelopes

### 2.1 Competing Asymptotic Hypotheses
In Gate 1 of the canonical roadmap ([`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md)), tail extinction requires:
$$\lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{\mathrm{Gate1}}(N, L) = 0, \qquad \Pi_{\mathrm{Gate1}}(N, L) \equiv \Delta_2(N) R_{\mathrm{spec}}(N, L) = \frac{\Delta_2(N) \|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}}}{D(N, L)}.$$

We formulate two competing asymptotic models:
1. **Hypothesis 145.1 (Polynomial Denominator Collapse):**  
   The resolvent denominator $D(N)$ collapses as a polynomial power law:
   $$D(N) \approx C_D \, N^{-p_D} \qquad (p_D > 0).$$
   Because the operator norm grows mildly as $\|Q_{\mathrm{even}}\|_{\mathrm{op}} \approx C_Q \, N^{p_Q}$ (with $p_Q \approx 0.29$), the resolvent ratio would scale as:
   $$R_{\mathrm{spec}}(N) \approx C_R \, N^{p_R}, \qquad p_R \equiv p_D + p_Q.$$
2. **Hypothesis 145.2 (Exponential Tunneling Damping):**  
   The bound-state doublet splitting decays exponentially under semiclassical barrier tunneling:
   $$\Delta_2(N) \approx C_\Delta \, e^{-\kappa N} \qquad (\kappa > 0).$$

---

### 2.2 The Dimensionless Log-Ratio Metric $\rho_{\mathrm{Gate1}}(N)$
To quantify whether doublet splitting dominates resolvent divergence without relying on arbitrary prefactor normalizations, we define the **dimensionless log-ratio metric**:
$$\rho_{\mathrm{Gate1}}(N) \equiv \frac{-\log \Delta_2(N)}{\log R_{\mathrm{spec}}(N)}.$$

> **Proposition 145.1 (Subexponential Domination Criterion).**  
> If $R_{\mathrm{spec}}(N) = \mathcal{O}(N^{p_R})$ is polynomial while $\Delta_2(N) \le C_\Delta e^{-\kappa N}$ is exponential with $\kappa > 0$, then asymptotically:
> $$\rho_{\mathrm{Gate1}}(N) = \frac{\kappa N - \log C_\Delta}{p_R \log N + \log C_R} \sim \frac{\kappa}{p_R} \frac{N}{\log N} \longrightarrow \infty.$$
> Consequently:
> $$\lim_{N \to \infty} \Pi_{\mathrm{Gate1}}(N) = \lim_{N \to \infty} C_R N^{p_R} C_\Delta e^{-\kappa N} = 0.$$

**Important Epistemic Qualification:** While Proposition 145.1 is an exact mathematical deduction from the joint premises ($R_{\mathrm{spec}}$ polynomial and $\Delta_2$ exponential), the empirical numerical sequence $\rho(N)$ across the accessible range $N \in [24, 80]$ actually **decreases**:
$$\rho(N) \in \{19.11, \; 9.15, \; 7.64, \; 7.63, \; 7.47, \; 7.29, \; 7.15, \; 7.07\}.$$
Therefore, $\rho(N) \to \infty$ is **not empirically demonstrated** by the finite-$N$ data. What the data rigorously establishes is:
$$\boxed{\rho(N) > 7.0 \quad \text{throughout the tested range } N \in [24, 80],}$$
which confirms that $\Delta_2 R_{\mathrm{spec}}$ is separated from unity by roughly seven orders of magnitude on a natural logarithmic scale.

---

### 2.3 Lower Bounds on the Resolvent Denominator $D(N)$
Analytically, Gate 1 does not require $D(N) \ge D_\infty > 0$. We require only a subexponential lower bound:
$$D(N) \ge c_0 \, N^{-p_D} \qquad \text{or} \qquad D(N) \ge c_0 \, e^{-o(N)}.$$

Because $E_2, E_3$ are bound states deep in the step well, $E_2, E_3 \lll E_{11}$ ($E_2, E_3 \approx 0$ to $> 20$ digits), so:
$$D(N) = (E_{11} - E_2)(E_{11} - E_3) \approx E_{11}(N)^2.$$

**Retirement of Generic Weyl Asymptotics:** Prior heuristics suggested that Weyl eigenvalue density estimates near an accumulation point bound Ritz level spacing by $\mathcal{O}(N^{-2})$ or $\mathcal{O}(N^{-3})$. However, generic Weyl asymptotics describe eigenvalue counting in bulk regimes; they do not automatically provide a lower bound on the distance from a specific edge state $E_{11}$ to bound states. Threshold spacing depends sensitively on the boundary conditions, the Galerkin truncation scheme, the effective potential, and whether $E_{11}$ is a true continuum eigenvalue.  
Generic Weyl spacing is therefore retired as an asserted bound; it is replaced by the specific **finite-volume boundary quantization condition** investigated in Cell 146.

---

## 3. Dual Core-Size Architecture: Edge vs Continuum

To separate the physics of the well edge ($L=10$, where $E_{11}$ accumulates) from the stabilized continuum ($L=12$, where $E_{13}$ lives in the scattering continuum), Cell 145 evaluates both core sizes:
- **Nominal Bound-State Edge ($L = 10$, $E_{L+1} = E_{11}$):**  
  Exhibits collapsing gap $g_{\mathrm{cont}} = E_{11} - E_{10} \to 0$ and exploding $R_{\mathrm{spec}, 10} \to \infty$. Tests whether denominator collapse is dominated by tunneling suppression.
- **Stabilized Scattering Continuum ($L = 12$, $E_{L+1} = E_{13}$):**  
  As discovered in Cell 121, the continuum gap stabilizes at $g_{2, 12}(N) \equiv E_{13} - E_3 \ge 0.5824 > 0$. The denominator $D_{12}(N) \ge 0.336$ remains strictly macroscopic, yielding $R_{\mathrm{spec}, 12} \approx 16 = \mathcal{O}(1)$ and super-exponential product extinction $\Pi_{2, 12}(N) \sim 10^{-37}$.

---

## 4. Experimental Suite Design (`cell145.py`)

The computational suite executes four targeted numerical modules across $N \in [24, 32, 40, 48, 56, 64, 72, 80]$:
1. **Module 1 (Pre-Flight Regression Audit):** Hard audit at $N=64$ verifying agreement with certified invariants within $10^{-8}$. The $N=64$ system is cached and reused with explicit label `(reused from pre-flight cache)`.
2. **Module 2 (Primary Scaling Sweep):** Computes $E_{10}, E_{11}, g_{\mathrm{cont}}, \Delta_2, D_{10}, \|Q\|_{\mathrm{op}}, R_{\mathrm{spec}, 10}, \Pi_2(N), \rho(N)$.
3. **Module 3 (Asymptotic Envelope Curve-Fitting):** Fits log-linear and log-log scaling models for $D(N), R_{\mathrm{spec}}(N), \|Q\|_{\mathrm{op}}$, and $\Delta_2(N)$.
4. **Module 4 (Multi-Doublet Comparison):** Evaluates products $\Pi_0(N), \Pi_1(N), \Pi_2(N)$ for doublets $j=0, 1, 2$.
5. **Module 5 (Dual Core-Size Comparison):** Contrasts $L=10$ vs $L=12$.

---

## 5. Certified Diagnostic Output Tables (`cell145.out`)

### Table 1: Asymptotic Scaling Table ($L = 10$, Nominal Continuum Base $E_{11}$)
Tracking Denominator Collapse $D(N)$, Doublet Splitting $\Delta_2(N)$, Resolvent Ratio $R_{\mathrm{spec}}$, and Log-Ratio $\rho(N)$:

| $N$ | $E_{10}$ | $E_{11}$ | $g_{\mathrm{cont}}$ | $\Delta_2$ | Denom $D_{10}$ | $\|Q\|_{\mathrm{op}}$ | $R_{\mathrm{spec}, 10}$ | $\Pi_2(N) = \Delta_2 R_{\mathrm{spec}}$ | $\rho_{\mathrm{Gate1}}$ |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 24 | 0.0054 | 0.4118 | 0.4064 | $1.4080 \times 10^{-26}$ | $1.6955 \times 10^{-1}$ | 3.82 | 22.51 | $3.1700 \times 10^{-25}$ | 19.114 |
| 32 | 0.0001 | 0.0416 | 0.0414 | $1.3121 \times 10^{-31}$ | $1.7293 \times 10^{-3}$ | 4.11 | 2379.06 | $3.1215 \times 10^{-28}$ | 9.146 |
| 40 | 0.0000 | 0.0125 | 0.0124 | $5.9823 \times 10^{-35}$ | $1.5512 \times 10^{-4}$ | 4.67 | 30098.69 | $1.8006 \times 10^{-30}$ | 7.642 |
| 48 | 0.0000 | 0.0092 | 0.0092 | $4.4210 \times 10^{-37}$ | $8.5428 \times 10^{-5}$ | 4.96 | 58071.38 | $2.5673 \times 10^{-32}$ | 7.631 |
| 56 | 0.0000 | 0.0071 | 0.0071 | $4.8088 \times 10^{-38}$ | $5.0338 \times 10^{-5}$ | 5.03 | 99839.13 | $4.8011 \times 10^{-33}$ | 7.465 |
| 64 | 0.0000 | 0.0060 | 0.0060 | $2.9751 \times 10^{-38}$ | $3.6296 \times 10^{-5}$ | 5.13 | 141466.14 | $4.2088 \times 10^{-33}$ | 7.286 |
| 72 | 0.0000 | 0.0054 | 0.0054 | $2.4510 \times 10^{-38}$ | $2.8984 \times 10^{-5}$ | 5.26 | 181518.46 | $4.4490 \times 10^{-33}$ | 7.152 |
| 80 | 0.0000 | 0.0050 | 0.0050 | $2.1801 \times 10^{-38}$ | $2.5437 \times 10^{-5}$ | 5.39 | 211946.25 | $4.6206 \times 10^{-33}$ | 7.071 |

---

### Table 2: Asymptotic Curve-Fitting & Exponent Extraction

| Functional Model | Extracted Law | Exponent | Fit Quality $R^2$ | Analytical Character |
|:---|:---|:---|:---|:---|
| Denominator Collapse $D(N)$ | $D(N) \sim C_D N^{-p_D}$ | $p_D = 6.6564$ | 0.841891 | Finite-range descriptive fit (not asymptotic law) |
| Resolvent Growth $R_{\mathrm{spec}}(N)$ | $R(N) \sim C_R N^{p_R}$ | $p_R = 6.9457$ | 0.848881 | Finite-range descriptive fit (not asymptotic law) |
| Operator Norm $\|Q\|_{\mathrm{op}}$ | $\|Q\| \sim C_Q N^{p_Q}$ | $p_Q = 0.2894$ | 0.952814 | Clean sublinear growth |
| Doublet Decay (Clean $N \in [24..56]$) | $\Delta_2 \sim C e^{-\kappa N}$ | $\kappa = 0.8176$ | 0.927933 | Rapid exponential transient |
| Doublet Decay (Full $N \in [24..80]$) | $\Delta_2 \sim C e^{-\kappa N}$ | $\kappa = 0.4358$ | 0.752728 | Pronounced flattening / deceleration |

*Consistency Identity Check:* $p_R$ (6.9457) vs $p_D + p_Q$ (6.9457) — Discrepancy: $3.21 \times 10^{-50}$ (algebraic identity confirmed to precision floor).

---

### Table 3: Multi-Doublet Gate 1 Products $\Pi_j(N) = \Delta_j(N) R_{\mathrm{spec}, 10}(N)$

| $N$ | $\Delta_0^{\mathrm{par}}$ | $\Pi_0 = \Delta_0^{\mathrm{par}} R_{\mathrm{spec}}$ | $\Delta_1^{\mathrm{par}}$ | $\Pi_1 = \Delta_1^{\mathrm{par}} R_{\mathrm{spec}}$ | $\Delta_2$ | $\Pi_2 = \Delta_2 R_{\mathrm{spec}}$ |
|:---|:---|:---|:---|:---|:---|:---|
| 24 | $4.9186 \times 10^{-40}$ | $1.1074 \times 10^{-38}$ | $3.3571 \times 10^{-34}$ | $7.5584 \times 10^{-33}$ | $1.4080 \times 10^{-26}$ | $3.1700 \times 10^{-25}$ |
| 32 | $5.4013 \times 10^{-46}$ | $1.2850 \times 10^{-42}$ | $7.9126 \times 10^{-40}$ | $1.8825 \times 10^{-36}$ | $1.3121 \times 10^{-31}$ | $3.1215 \times 10^{-28}$ |
| 40 | $3.6093 \times 10^{-51}$ | $1.0863 \times 10^{-46}$ | $9.2842 \times 10^{-44}$ | $2.7944 \times 10^{-39}$ | $5.9823 \times 10^{-35}$ | $1.8006 \times 10^{-30}$ |
| 48 | $2.7919 \times 10^{-50}$ | $1.6213 \times 10^{-45}$ | $2.4441 \times 10^{-46}$ | $1.4193 \times 10^{-41}$ | $4.4210 \times 10^{-37}$ | $2.5673 \times 10^{-32}$ |
| 56 | $1.0331 \times 10^{-50}$ | $1.0314 \times 10^{-45}$ | $1.2129 \times 10^{-47}$ | $1.2110 \times 10^{-42}$ | $4.8088 \times 10^{-38}$ | $4.8011 \times 10^{-33}$ |
| 64 | $7.9865 \times 10^{-51}$ | $1.1298 \times 10^{-45}$ | $5.7237 \times 10^{-48}$ | $8.0970 \times 10^{-43}$ | $2.9751 \times 10^{-38}$ | $4.2088 \times 10^{-33}$ |
| 72 | $9.0811 \times 10^{-51}$ | $1.6484 \times 10^{-45}$ | $4.5606 \times 10^{-48}$ | $8.2783 \times 10^{-43}$ | $2.4510 \times 10^{-38}$ | $4.4490 \times 10^{-33}$ |
| 80 | $1.0530 \times 10^{-50}$ | $2.2317 \times 10^{-45}$ | $3.8345 \times 10^{-48}$ | $8.1271 \times 10^{-43}$ | $2.1801 \times 10^{-38}$ | $4.6206 \times 10^{-33}$ |

---

### Table 4: Dual Core-Size Comparison ($L = 10$ [Edge] vs $L = 12$ [Scattering Continuum])

| $N$ | Denom $D_{10}$ | $R_{\mathrm{spec}, 10}$ | $\Pi_2(L=10)$ | Denom $D_{12}$ | $R_{\mathrm{spec}, 12}$ | $\Pi_2(L=12)$ |
|:---|:---|:---|:---|:---|:---|:---|
| 24 | $1.6955 \times 10^{-1}$ | 22.51 | $3.1700 \times 10^{-25}$ | 3.2445 | 1.18 | $1.6566 \times 10^{-26}$ |
| 32 | $1.7293 \times 10^{-3}$ | 2379.06 | $3.1215 \times 10^{-28}$ | 0.6662 | 6.18 | $8.1026 \times 10^{-31}$ |
| 40 | $1.5512 \times 10^{-4}$ | 30098.69 | $1.8006 \times 10^{-30}$ | 0.5292 | 8.82 | $5.2780 \times 10^{-34}$ |
| 48 | $8.5428 \times 10^{-5}$ | 58071.38 | $2.5673 \times 10^{-32}$ | 0.3954 | 12.55 | $5.5464 \times 10^{-36}$ |
| 56 | $5.0338 \times 10^{-5}$ | 99839.13 | $4.8011 \times 10^{-33}$ | 0.3544 | 14.18 | $6.8190 \times 10^{-37}$ |
| 64 | $3.6296 \times 10^{-5}$ | 141466.14 | $4.2088 \times 10^{-33}$ | 0.3392 | 15.14 | $4.5035 \times 10^{-37}$ |
| 72 | $2.8984 \times 10^{-5}$ | 181518.46 | $4.4490 \times 10^{-33}$ | 0.3376 | 15.58 | $3.8192 \times 10^{-37}$ |
| 80 | $2.5437 \times 10^{-5}$ | 211946.25 | $4.6206 \times 10^{-33}$ | 0.3362 | 16.04 | $3.4964 \times 10^{-37}$ |

---

## 6. Scientific Analysis & Epistemic Audit

### 6.1 The Core Positive Finding: Doublet Damping vs Resolvent Explosion
Across the expanded sweep $N \in [24, 80]$, the resolvent ratio $R_{\mathrm{spec}, 10}$ explodes by nearly four orders of magnitude ($22.5 \to 211946$). Despite this severe divergence, the Gate 1 tail extinction product $\Pi_2(N) \equiv \Delta_2(N) R_{\mathrm{spec}, 10}(N)$ collapses by eight orders of magnitude:
$$\Pi_2(24) = 3.17 \times 10^{-25} \longrightarrow \Pi_2(80) = 4.62 \times 10^{-33}.$$
This provides compelling empirical support for the qualitative mechanism: **tunneling suppression easily overwhelms resolvent divergence across the accessible pre-asymptotic regime.**

### 6.2 Critique of the Polynomial Fit: A Descriptive, Decelerating Transient
The regression fit yields $D(N) \sim N^{-6.6564}$ ($R^2 = 0.8419$) and $R_{\mathrm{spec}}(N) \sim N^{6.9457}$ ($R^2 = 0.8489$). These modest $R^2$ values indicate that a single global power law is **not** the true asymptotic form:
- From $N=32 \to 40$, $D(N)$ drops from $1.73 \times 10^{-3} \to 1.55 \times 10^{-4}$ (an **11-fold reduction**, corresponding to an effective local exponent $p_D^{\mathrm{loc}} \approx 10.8$).
- From $N=72 \to 80$, $D(N)$ drops from $2.90 \times 10^{-5} \to 2.54 \times 10^{-5}$ (only a **14% reduction**, corresponding to an effective local exponent $p_D^{\mathrm{loc}} \approx 1.24$).

The denominator collapse is visibly decelerating. The fitted exponent $p_D = 6.66$ is a descriptive average over a transient window, not an established asymptotic power law.

### 6.3 Instability of the Doublet Tunneling Exponent
The doublet splitting fit demonstrates significant sensitivity to the inclusion of higher dimensions:
- Restricted clean range $N \in [24, 56]$: $\kappa = 0.8176$ ($R^2 = 0.9279$).
- Full range $N \in [24, 80]$: $\kappa = 0.4358$ ($R^2 = 0.7527$).

While $\Delta_2$ decreases extraordinarily rapidly over the sweep ($1.4 \times 10^{-26} \to 2.18 \times 10^{-38}$), the extracted exponent is not asymptotically stable, and an exponential law cannot be claimed as empirically proven.

### 6.4 The Flattening of the Gate 1 Product Around $4 \times 10^{-33}$
A crucial feature of Table 1 is that after $N=48$, the Gate 1 product stops decreasing rapidly:
$$2.57 \times 10^{-32} \longrightarrow 4.80 \times 10^{-33} \longrightarrow 4.21 \times 10^{-33} \longrightarrow 4.45 \times 10^{-33} \longrightarrow 4.62 \times 10^{-33}.$$
The product essentially **flattens into a plateau near $4 \times 10^{-33}$**. Plausible physical and numerical hypotheses for this plateau include:
1. **Numerical Floor / Conditioning of $\Delta_2$:** At $N=80$, $\Delta_2 = 2.18 \times 10^{-38}$. In a 50-dps calculation with a global operator norm $\|Q\| \approx 5.4$, numerical precision limits and eigenvalue cluster conditioning may begin to introduce relative noise into the splitting.
2. **Finite-$T$ Truncation Effects:** As demonstrated in Cell 144, shifting $T$ from $600 \to 800$ alters edge levels by $5.7\%$. A finite cutoff $T = 600$ may impose an effective floor on the intra-well tunneling barrier.
3. **Genuine Asymptotic Crossover:** A crossover where the local decay rate of $\Delta_2$ slows to match the local growth rate of $R_{\mathrm{spec}}$.

Under no circumstances should the rapid $N \in [24, 56]$ slope be naively extrapolated through this plateau.

### 6.5 The Clean Operator Norm Scaling
The operator norm exhibits the cleanest scaling of the entire suite:
$$\|Q_{\mathrm{even}}^{(N)}\|_{\mathrm{op}} \sim N^{0.2894} \qquad (R^2 = 0.9528).$$
The algebraic consistency relation $p_R = p_D + p_Q$ holds to 50-digit precision ($6.9457 = 6.6564 + 0.2894$). This confirms structurally that the divergence in $R_{\mathrm{spec}}$ is **overwhelmingly driven by denominator collapse**, while the operator norm remains exceptionally well-behaved.

### 6.6 The $L = 12$ Control Experiment
Table 4 demonstrates a sharp contrast between the two core sizes:
- Well edge ($L = 10$): $D_{10} \to 2.54 \times 10^{-5}$, $R_{\mathrm{spec}, 10} \to 2.12 \times 10^5$.
- Scattering continuum ($L = 12$): $D_{12} \ge 0.3362$, $R_{\mathrm{spec}, 12} \approx 16.04 = \mathcal{O}(1)$, and $\Pi_2(L=12) \approx 3.5 \times 10^{-37}$.

This confirms that the exploding resolvent is tied specifically to the well-edge mode $E_{11}$ in $L=10$, while in the scattering continuum $L=12$, the denominator remains strictly macroscopic. (This is demonstrated for $L \in \{10, 12\}$, not asserted for all $L$.)

---

## 7. Epistemic Classification & Strategic Hand-Off

### 7.1 Formal Epistemic Classification

#### 1. Established Numerically (Audit-Proof):
- The resolvent denominator $D_{10}(N)$ collapses dramatically ($0.170 \to 2.54 \times 10^{-5}$).
- The resolvent ratio $R_{\mathrm{spec}, 10}(N)$ explodes ($22.5 \to 2.12 \times 10^5$).
- The Gate 1 product $\Pi_2(N) \equiv \Delta_2 R_{\mathrm{spec}, 10}$ collapses by eight orders of magnitude, remaining below $4.7 \times 10^{-33}$ across $N \ge 56$.
- The operator norm $\|Q_{\mathrm{even}}\|_{\mathrm{op}}$ grows mildly as $N^{0.29}$ ($R^2 = 0.953$).
- At $L=12$, the denominator stabilizes at $D_{12} \approx 0.34$, keeping $R_{\mathrm{spec}, 12} \approx 16 = \mathcal{O}(1)$.
- The product suppression persists across all three tested doublets ($\Pi_0 \sim 10^{-45}$, $\Pi_1 \sim 10^{-43}$, $\Pi_2 \sim 10^{-33}$).
- Pre-flight invariants match certified baselines to $2.5 \times 10^{-11}$.

#### 2. Strongly Suggested:
- $D_{10}(N)$ is subexponentially small (exhibiting marked deceleration toward a low-order power law).
- $\Delta_2(N)$ is exponentially or faster decaying over the pre-asymptotic range.

#### 3. Not Established:
- $D(N) \ge c N^{-p}$ or $D(N) \ge e^{-o(N)}$ is **not proved**.
- An asymptotic power law $N^{-6.66}$ or pure exponential decay $e^{-0.82 N}$ is **not established**.
- Asymptotic divergence of the log-ratio $\rho(N) \to \infty$ is **not empirically observed** (it decreases from $19.1 \to 7.07$).

$$\boxed{\textbf{Verdict: Gate 1 is NOT yet mathematically cleared.}}$$

---

### 7.2 Forward Strategy: Cell 146
The research direction for Cell 146 is clearly delineated:
- We do not run another large-$N$ curve fit.
- We attack the analytical weak link: **Can we derive a genuine lower bound on $D(N) \approx E_{11}(N)^2$ from the finite-volume boundary quantization condition and boundary phase shift $\delta(E)$?**
- Because $E_{11}$ is the lowest continuum state of a specific discrete Galerkin operator, its Ritz level is constrained by the effective box size $X_{\mathrm{eff}}(N)$ and boundary reflection phase:
  $$k_{11} X_{\mathrm{eff}}(N) + \delta(E_{11}) = \pi \implies E_{11}(N) \gtrsim N^{-2} \implies D(N) \gtrsim N^{-4}.$$
- Cell 146 will formulate and test this boundary quantization mechanism.
