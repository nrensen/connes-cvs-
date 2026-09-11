# Cell 105 Analytical Note: The Quantitative Localization Programme for the Connes–CvS Ground State

**Companion Computational Script:** [`cell105.py`](file:///c:/data/github/connes-cvs-/cell105.py) | **Verification Log:** [`cell105.out`](file:///c:/data/github/connes-cvs-/cell105.out)  
**Status:** Working Research Note (Gate 1 / Route 1D / Phase VIII)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell104.md`](file:///c:/data/github/connes-cvs-/cell104.md); [`cell103.md`](file:///c:/data/github/connes-cvs-/cell103.md); Paper NR2 Section 8; Cells 98–104  
**Date:** September 2026  

---

## Executive Summary

Cell 104 established the exact structural reduction of the Feshbach decoupling problem:
$$\boxed{\Delta E_{11}^{\mathrm{Fesh}}(M) \asymp \|v^{(Q)}\|^2, \qquad \|w_M\|^2 \asymp \|v^{(Q)}\|^2.}$$
By proving the **Exact Eigenvector Complementarity Identity** ($w_M = -(C_M - E_{11} I) v^{(Q)}$) and the **Exact Rayleigh Shift Energy Identity** ($\Delta E_{11}^{\mathrm{Fesh}} = \langle v^{(Q)}, (C_M - E_{11} I) v^{(Q)} \rangle$), Cell 104 completely eliminated the Feshbach resolvent inverse $(C_M - E_{11} I)^{-1}$ as an analytical obstacle. The back-reaction onto the ground state does not suffer from small-denominator singularities.

However, as the reviewer's evaluation of Cell 104 precisely emphasized:
$$\boxed{\textbf{Exact complementarity does not by itself prove that the tail is exponentially small.}}$$
The logical chain is strictly:
$$\boxed{\text{exact finite-}N\text{ complementarity} \;+\; \text{ground-state localization } \|v^{(Q)}\|^2 \to 0 \;\Longrightarrow\; \text{Feshbach decoupling } \Delta E_{11}^{\mathrm{Fesh}} \to 0.}$$

The entire Gate 1 Feshbach decoupling programme has therefore collapsed onto a single, universal analytical question:
$$\boxed{\textbf{Why is the ground-state eigenvector } v_N \textbf{ so extraordinarily localized, and can we prove it quantitatively?}}$$

This note launches the **Quantitative Localization Programme for $v_N$**, systematically examining four analytical routes:
1. **Bernstein–Paley–Wiener Continuum Embedding (Route A):** Boundary flatness $T_\infty^{(k)}(0) = T_\infty^{(k)}(L) = 0$ in physical space rigorously yields super-polynomial decay $|v_m| = \mathcal{O}(m^{-k})$ for all $k \ge 1$; analyticity in a complex Bernstein strip yields exponential decay $|v_m| \le C e^{-\sigma m}$.
2. **Discrete Combes–Thomas Resolvent Inversion (Route B):** Inverting the complementarity identity gives $v^{(Q)} = -(C_M - E_{11} I)^{-1} w_M$; off-diagonal decay in $C_M$ combined with the spectral gap $\delta_M \ge \delta_\infty > 0$ generates exponential spatial decay in the resolvent kernel.
3. **The Coordinate Inversion Recurrence (Route C):** Exact diagonal balancing $v_k = -\frac{1}{H_{kk} - E_{11}} \sum_{m \neq k} H_{mk} v_m$ governs how modal amplitudes propagate from the low-frequency core into the high-frequency tail.
4. **Variational Sobolev Kinetic Energy Bounds (Route D):** Unconditional, deterministic polynomial tail enclosures $\|v^{(Q)}\|^2 \le M^{-2s} \|T_v\|_{H^s}^2$ derived directly from the near-zero Rayleigh quotient $\langle v_N, H v_N \rangle \approx 0$.

---

## 1. Epistemic Calibration & The Exact Problem Formulation

### 1.1 The Epistemic Boundary
Following the audit of Cell 104, we maintain strict epistemic discipline:
- **Established Mathematical Fact (Tier 1):** The two-sided tail-mass sandwich (Theorem 3 of Cell 104) is an exact algebraic theorem valid for all finite $N \ge 1$:
  $$\delta_M \|v^{(Q)}\|^2 \le \Delta E_{11}^{\mathrm{Fesh}}(M) \le (\|H\|_{\mathrm{op}} - E_{11}) \|v^{(Q)}\|^2.$$
- **Empirical Numerical Fact:** Across tested cutoffs $M \in \{24, 32, 48, 64\}$ at $N=192$, the tail mass plummets across 26 decimal orders of magnitude:
  $$\|v^{(Q)}(24)\|^2 \approx 9.22 \times 10^{-26} \;\longrightarrow\; \|v^{(Q)}(64)\|^2 \approx 2.70 \times 10^{-52}.$$
- **Analytical Target (Active Problem):** Prove an unconditional upper bound:
  $$\mathcal{T}(M) \equiv \|v^{(Q)}(M)\|^2 = \sum_{m=M+1}^N v_m^2 \le C f(M),$$
  where $f(M) \to 0$ as $M \to \infty$ (at least polynomially $f(M) = M^{-2s}$, and ideally exponentially $f(M) = e^{-2\sigma M}$).

### 1.2 The Conceptual Resolution: Cancellation and Tail Localization are Duals
As the reviewer precisely observed, Cell 104 resolved the paradox of the tiny projected coupling $w_M = B_M^T v^{(P)}$:
- Previously, we looked at $w_M = B_M^T v^{(P)}$ and wondered how a vector formed from algebraically decaying matrix entries ($H_{mk} \sim 1/k$) could have norm $\|w_M\| \sim 10^{-26}$ (norm squared $10^{-51}$).
- We thought we had to determine whether this was driven by kernel decay, ground-state localization, or miraculous destructive interference.
- Theorem 1 reveals that:
  $$w_M = -(C_M - E_{11} I) v^{(Q)} \quad \Longleftrightarrow \quad \sum_{m=0}^M H_{mk} v_m = -(H_{kk} - E_{11}) v_k - \sum_{\substack{m=M+1 \\ m \neq k}}^N H_{mk} v_m.$$
  At the level of the exact finite matrix, **"cancellation" and "tail localization" are not competing explanations; they are two descriptions of the exact same high-sector eigenvector equation.**
- The low-mode terms really are individually substantial (at $M=64$, $\sum_{m \le M} |v_m H_{m, M+1}| \approx 3.03 \times 10^{-2}$), while their signed sum is $1.26 \times 10^{-26}$, yielding a cancellation ratio of $2.4 \times 10^{24}$. But that cancellation is *forced* by the eigenvector equation because the tail $v^{(Q)}$ is of order $10^{-26}$.

### 1.3 The Observed Tail Rate Drift
In `cell104.out`, the empirical exponential decay rate $\sigma_M \equiv -\frac{1}{2M} \log(\|v^{(Q)}\|^2)$ exhibited a slow, monotonic decline:
$$\sigma_{24} \approx 1.201 \;\longrightarrow\; \sigma_{32} \approx 1.172 \;\longrightarrow\; \sigma_{48} \approx 1.112 \;\longrightarrow\; \sigma_{64} \approx 0.928.$$
This decline indicates that the asymptotic decay is not a pure single-exponent exponential $e^{-\sigma m}$, but possesses sub-leading algebraic or logarithmic modulation:
$$|v_m| \approx C m^{-\gamma} e^{-\sigma m} \quad \implies \quad \sigma_m = \sigma + \gamma \frac{\log m}{m}.$$
As $m$ increases, the positive correction $\gamma \frac{\log m}{m}$ decreases toward zero, causing $\sigma_m$ to drift downward toward the true asymptotic rate $\sigma$.
Distinguishing between a modulated exponential ($|v_m| \sim m^{-\gamma} e^{-\sigma m}$), a stretched exponential ($|v_m| \sim e^{-\sigma m^\beta}$), and a pure power law ($|v_m| \sim m^{-p}$) is a central numerical objective of `cell105.py`.

---

## 2. Route A: Bernstein–Paley–Wiener Continuum Embedding

### 2.1 Trigonometric Duality
The normalized eigenvector $v_N = (v_0, v_1, \dots, v_N)^T \in \mathbb{R}^{N+1}$ defines a unique real-valued trigonometric polynomial on the interval $[0, L]$ (where $L = \log c$):
$$T_{v_N}(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos\left(\frac{2\pi m t}{L}\right).$$
By Parseval's identity on $L^2([0, L])$:
$$\|T_{v_N}\|_{L^2}^2 = \int_0^L |T_{v_N}(t)|^2 \, dt = L \left( v_0^2 + \sum_{m=1}^N v_m^2 \right) = L \|v_N\|_2^2 = L.$$

The Fourier coefficients are recovered by the inversion formula:
$$v_0 = \frac{1}{L} \int_0^L T_{v_N}(t) \, dt, \qquad v_m = \frac{\sqrt{2}}{L} \int_0^L T_{v_N}(t) \cos\left(\frac{2\pi m t}{L}\right) dt \quad (m \ge 1).$$

### 2.2 Boundary Contact and Super-Polynomial Decay
Let $T(t)$ be a smooth function on $[0, L]$. Integrate the Fourier formula by parts:
$$v_m = \frac{\sqrt{2}}{L} \left[ T(t) \frac{\sin(a_m t)}{a_m} \right]_0^L - \frac{\sqrt{2}}{L a_m} \int_0^L T'(t) \sin(a_m t) \, dt, \qquad a_m = \frac{2\pi m}{L}.$$
Since $\sin(a_m 0) = \sin(a_m L) = \sin(2\pi m) = 0$, the boundary term vanishes identically for all $m \ge 1$.

Integrating by parts a second time:
$$v_m = \frac{\sqrt{2}}{L a_m^2} \left[ T'(t) \cos(a_m t) \right]_0^L - \frac{\sqrt{2}}{L a_m^2} \int_0^L T''(t) \cos(a_m t) \, dt.$$
Here the boundary term evaluates to:
$$\left[ T'(t) \cos(a_m t) \right]_0^L = T'(L) \cos(2\pi m) - T'(0) \cos(0) = T'(L) - T'(0).$$

**Theorem 1 (Algebraic Decay under Boundary Contact).**
*Let $k \ge 1$ be an integer. Suppose the continuum solitary wave profile $T_\infty \in C^{2k}([0, L])$ satisfies dual boundary contact:*
$$T_\infty^{(2j+1)}(L) - T_\infty^{(2j+1)}(0) = 0 \qquad (\forall j \in \{0, \dots, k-1\}).$$
*Then its Fourier coefficients satisfy:*
$$\boxed{|v_m| \le \frac{\sqrt{2}}{L} \left(\frac{L}{2\pi m}\right)^{2k} \|T_\infty^{(2k)}\|_{L^1} = \mathcal{O}(m^{-2k}).}$$
*Furthermore, if $T_\infty$ has infinite-order flat boundary contact (Paper NR2 Conjecture 2):*
$$T_\infty^{(j)}(0) = T_\infty^{(j)}(L) = 0 \qquad (\forall j \ge 0),$$
*then its periodic extension $\widetilde{T}_\infty \in C_c^\infty(\mathbb{R})$, and the Fourier coefficients decay faster than every polynomial power:*
$$\boxed{|v_m| = \mathcal{O}(m^{-\infty}) \quad \Longleftrightarrow \quad \forall k \ge 1, \; \exists C_k < \infty : |v_m| \le \frac{C_k}{m^k}.}$$

*Proof.*
Follows directly from $2k$-fold integration by parts, where all boundary terms vanish identically by hypothesis. $\quad \blacksquare$

### 2.3 Complex Analyticity and Exponential Decay (Paley–Wiener)
If the solitary wave $T_\infty(t)$ is real analytic on $[0, L]$ and can be analytically continued into a complex strip of width $\delta > 0$:
$$\mathcal{S}_\delta \equiv \{ z \in \mathbb{C} : 0 \le \operatorname{Re} z \le L, \; |\operatorname{Im} z| < \delta \},$$
with periodic boundary conditions across the strip, then shifting the contour of integration to $t \pm i\delta$ in the complex plane yields the classical Paley–Wiener bound:
$$\boxed{|v_m| \le C_\delta e^{-\frac{2\pi \delta}{L} m} = C_\delta e^{-\sigma m}, \qquad \sigma \equiv \frac{2\pi \delta}{L}.}$$
At $c = 13$, $L = \log 13 \approx 2.5649$. An empirical decay rate of $\sigma \approx 1.0$ corresponds to a strip of analyticity of width:
$$\delta = \frac{\sigma L}{2\pi} \approx \frac{1.0 \times 2.5649}{6.2832} \approx 0.408.$$
This physical strip width is approximately $16\%$ of the box size $L$.

---

## 3. Route B: Discrete Combes–Thomas Resolvent Bounds

### 3.1 Resolvent Inversion of the Complementarity Identity
Theorem 1 of Cell 104 established:
$$w_M \equiv B_M^T v^{(P)} = -(C_M - E_{11} I) v^{(Q)}.$$
Because $E_{11} < \mu_0(C_M)$, the operator $C_M - E_{11} I$ is strictly positive definite, with minimum eigenvalue $\delta_M \ge \delta_\infty \approx 0.892 > 0$.
Inverting this operator yields an **exact formula for the tail vector**:
$$\boxed{v^{(Q)} = -(C_M - E_{11} I)^{-1} w_M.}$$

Let $G_M \equiv (C_M - E_{11} I)^{-1}$ be the high-sector resolvent matrix. For each mode $k \in \{M+1, \dots, N\}$, the component $v_k$ is given by:
$$v_k = -\sum_{j=M+1}^N (G_M)_{k, j} w_{M, j}.$$

### 3.2 Off-Diagonal Decay of the High-Sector Resolvent
In classical condensed matter physics and operator theory, the Combes–Thomas estimate (1973) establishes that the resolvent $(H - E)^{-1}$ of an operator with local or rapidly decaying couplings has exponential off-diagonal decay whenever $E$ lies in a spectral gap.

**Proposition 2 (Combes–Thomas Resolvent Structure).**
*Let $C$ be a bounded symmetric operator on $\ell^2(\mathbb{N})$ with spectral gap $\operatorname{dist}(E, \sigma(C)) = \delta > 0$. If $C$ satisfies the exponential off-diagonal bound:*
$$|C_{j, k}| \le A e^{-\mu |j - k|},$$
*then for any $\eta < \min(\mu, \frac{\delta}{2 \|C\|_{\mathrm{op}}})$, there exists $K(\delta, \eta) < \infty$ such that the resolvent satisfies:*
$$\boxed{|[(C - E I)^{-1}]_{j, k}| \le K e^{-\eta |j - k|}.}$$

### 3.3 Application to Ground-State Boundary Extinction
Because $w_M = B_M^T v^{(P)}$ has entries $w_{M, j} = \sum_{m=0}^M H_{m, j} v_m$, its largest components lie near the interface $j \approx M+1$.
If the resolvent $G_M$ exhibits off-diagonal decay $|(G_M)_{k, j}| \le K e^{-\eta |k - j|}$, then for a mode $k \gg M$ deep in the high sector:
$$|v_k| \le \sum_{j=M+1}^N |(G_M)_{k, j}| |w_{M, j}| \le K \sum_{j=M+1}^N e^{-\eta (k - j)} |w_{M, j}| \le K' e^{-\eta (k - M)}.$$
This provides an autonomous discrete mechanism for exponential tail decay driven by the spectral gap $\delta_M$.

---

## 4. Route C: The Coordinate Inversion Recurrence

### 4.1 Exact Single-Mode Recurrence
Consider the full eigenvalue equation $H v_N = E_{11} v_N$ componentwise:
$$\sum_{m=0}^N H_{km} v_m = E_{11} v_k \qquad (\forall k \in \{0, \dots, N\}).$$
Isolating the diagonal term $m = k$:
$$H_{kk} v_k + \sum_{m \neq k} H_{km} v_m = E_{11} v_k.$$
Since $H_{kk} - E_{11} > 0$ for all $k \ge 1$:
$$\boxed{v_k = -\frac{1}{H_{kk} - E_{11}} \sum_{m \neq k} H_{km} v_m.}$$

### 4.2 Core-Tail Splitting of the Recurrence
For a high mode $k > M$, split the sum into low modes ($m \le M$) and high modes ($m > M$):
$$v_k = -\frac{1}{H_{kk} - E_{11}} \left( \sum_{m=0}^M H_{km} v_m + \sum_{\substack{m=M+1 \\ m \neq k}}^N H_{km} v_m \right).$$
Notice that the first sum is precisely the $k$-th component of the projected coupling vector $w_{M, k} = (B_M^T v^{(P)})_k$:
$$\boxed{v_k = -\frac{1}{H_{kk} - E_{11}} \left( w_{M, k} + \sum_{\substack{m=M+1 \\ m \neq k}}^N H_{km} v_m \right).}$$
This equation shows that the tail components $v_k$ satisfy a linear system driven by the source vector $w_M$, with kernel matrix $\frac{H_{km}}{H_{kk} - E_{11}}$ ($m \neq k$).

---

## 5. Route D: Variational Sobolev Kinetic Energy Bounds

### 5.1 Deterministic Polynomial Enclosure
Even without proving complex analyticity or exponential decay, we can prove **rigorous, unconditional polynomial tail bounds** from the Rayleigh quotient.

Let $s \ge 1$ be a positive integer. Define the discrete Sobolev $H^s$-norm of $v_N$:
$$\|v_N\|_{H^s}^2 \equiv \sum_{m=0}^N (1 + m^2)^s v_m^2.$$

**Theorem 3 (Sobolev Tail-Mass Enclosure).**
*For any cutoff $M \ge 1$ and any Sobolev index $s \ge 1$:*
$$\boxed{\|v^{(Q)}(M)\|^2 \equiv \sum_{m=M+1}^N v_m^2 \le \frac{1}{(1 + (M+1)^2)^s} \|v_N\|_{H^s}^2 \le \frac{\|v_N\|_{H^s}^2}{M^{2s}}.}$$

*Proof.*
Since $m \ge M+1$ for every term in $v^{(Q)}$:
$$(1 + (M+1)^2)^s \sum_{m=M+1}^N v_m^2 \le \sum_{m=M+1}^N (1 + m^2)^s v_m^2 \le \sum_{m=0}^N (1 + m^2)^s v_m^2 = \|v_N\|_{H^s}^2.$$
Dividing by $(1 + (M+1)^2)^s$ yields the result. $\quad \blacksquare$

### 5.2 Kinetic Moments of the Ground State
In physical space, the derivative operator $\partial_t$ corresponds to multiplication by $\frac{2\pi m}{L}$ in Fourier space:
$$\|T_v'\|_{L^2}^2 = \frac{4\pi^2}{L} \sum_{m=1}^N m^2 v_m^2, \qquad \|T_v''\|_{L^2}^2 = \frac{16\pi^4}{L} \sum_{m=1}^N m^4 v_m^2.$$
If the ground-state solitary wave has bounded kinetic energy $\|T_v\|_{H^s} \le K_s < \infty$ independent of $N$:
- At $s = 1$: $\|v^{(Q)}(M)\|^2 \le K_1^2 M^{-2}$ (polynomial rate 2).
- At $s = 2$: $\|v^{(Q)}(M)\|^2 \le K_2^2 M^{-4}$ (polynomial rate 4).
- At $s = 4$: $\|v^{(Q)}(M)\|^2 \le K_4^2 M^{-8}$ (polynomial rate 8).

### 5.3 The Spectral Confinement Mechanism: Archimedean Barrier in Mode Space
Why must the Sobolev kinetic moments $\mathcal{K}_s$ remain finite?
The Connes–van Suijlekom quadratic form contains André Weil's Archimedean distribution:
$$h_+(r) = \operatorname{Re} \psi\left(\frac{1}{4} + \frac{i r}{2}\right) - \log \pi.$$
For large frequencies $r \gg 1$, Stirling's approximation gives:
$$h_+(r) = \frac{1}{2} \log\left(\frac{1}{16} + \frac{r^2}{4}\right) - \log \pi + \mathcal{O}(r^{-2}) \approx \log r - \log(2\pi) > 0.$$
In the discrete Fourier Galerkin representation, the continuous frequency corresponds to $r_m = \frac{2\pi m}{L}$. Consequently, the diagonal entries of the Galerkin matrix grow monotonically with the mode index:
$$H_{mm} \sim \log(m) + \text{const} \quad (m \to \infty).$$
Now consider the ground-state Rayleigh expectation:
$$\langle v_N, H v_N \rangle = \sum_{m=0}^N H_{mm} v_m^2 + \sum_{j \neq k} H_{jk} v_j v_k = E_{11} \approx 0.$$
Because $E_{11}$ is extraordinarily close to zero ($E_{11} \approx -1.06 \times 10^{-51}$), the ground state operates under an extremely rigid energetic budget:
- Every mode $m$ populated by the ground state contributes a strictly positive diagonal cost $H_{mm} v_m^2 > 0$ that grows with $m$.
- To keep the total energy near zero, this diagonal penalty must be compensated by negative off-diagonal pairings $H_{jk} v_j v_k < 0$.
- However, as audited in Step 6, the off-diagonal couplings $|H_{jk}|$ decay away from the diagonal ($|H_{jk}| \sim \frac{1}{|j - k|}$). High modes $m \gg M$ have negligible direct coupling to the dominant low modes $j \in \{0, \dots, 8\}$.
- High modes therefore cannot offset their large positive diagonal cost $H_{mm} v_m^2$ through negative correlations with the low sector.
- In physical analogy with quantum mechanics, the Archimedean growth in Fourier space acts as an **effective confining potential barrier in mode space**, preventing probability leakage into high-frequency modes and enforcing uniform Sobolev boundedness $\sup_N \|v_N\|_{H^s} < \infty$.

---

## 6. Pre-Flight Specification for [`cell105.py`](file:///c:/data/github/connes-cvs-/cell105.py)

The companion computational script performs a comprehensive diagnostic suite at $N=192, c=13, T=600$:
1. **Modal Profile Extraction:** Extract the complete eigenvector components $v_m$ for $m = 0, \dots, 192$.
2. **Asymptotic Decay Fitting:** Fit the sequence $|v_m|$ across $m \in [10, 80]$ against four candidate functional forms:
   - Pure exponential: $\log |v_m| = c_0 - \sigma m$.
   - Modulated exponential: $\log |v_m| = c_0 - \sigma m - \gamma \log m$.
   - Stretched exponential: $\log |v_m| = c_0 - \sigma m^\beta$.
   - Pure power-law: $\log |v_m| = c_0 - p \log m$.
3. **Dense Tail Mass Sweep:** Evaluate the exact tail mass $\|v^{(Q)}(M)\|^2 = \sum_{m=M+1}^N v_m^2$ across a dense grid $M \in \{8, 12, 16, 20, 24, 32, 40, 48, 56, 64, 80, 96\}$.
4. **Sobolev Kinetic Moments:** Compute the moments $\mathcal{K}_s \equiv \sum_{m=1}^N m^{2s} v_m^2$ for $s \in \{1, 2, 3, 4, 6, 8\}$ and test their dimension-independence.
5. **Combes–Thomas Kernel Decay:** Compute the off-diagonal profile $|H_{j, k}|$ of the Galerkin matrix as a function of mode separation $d = |j - k|$ to audit the spatial locality of $H$.
