# Cell 128 Analytical Note: Translation-Defect and Boundary Mismatch Analysis for $Q_{\mathrm{even}}$ on $\Phi^\perp$

**Companion Document:** Analytical Research Note (Gate 1 / Milestone M-G1.5)  
**Status:** Working Analytical Research Note (Theoretical Formulation & Harmonic Analysis)  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.5: Variational Quasimode Codimension Bound](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell127.md`](file:///c:/data/github/connes-cvs-/cell127.md); [`cell126.md`](file:///c:/data/github/connes-cvs-/cell126.md); [`cell125.md`](file:///c:/data/github/connes-cvs-/cell125.md); [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md); [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py); [`Paper-NR1.md`](file:///c:/data/github/connes-cvs-/Paper-NR1.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)  
**Date:** September 2026  

---
## 1. Executive Summary & The Analytical Question

In [`cell127.md`](file:///c:/data/github/connes-cvs-/cell127.md), we established the exact shifted autocorrelation representation of the prime form and its polarization decomposition. Incorporating the exact $1/L$ normalization connecting the matrix form in [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py) to spatial integrals ($\|T_v\|_{L^2}^2 = L \|v\|_2^2$):
$$\langle v, Q_{\mathrm{prime}} v \rangle = -\frac{2}{L} \sum_{q = p^k \le c} w_q \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt = -\frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L T_v(t)^2 \, dt + \frac{1}{L} \sum_{q \le c} w_q \mathcal{D}_q[v] + \frac{1}{L} \sum_{q \le c} w_q \mathcal{B}_q[v],$$
where $w_q \equiv \frac{\Lambda(q)}{\sqrt{q}}$, $\mathcal{D}_q[v] \equiv \int_{\log q}^L |T_v(t) - T_v(t - \log q)|^2 \, dt \ge 0$, and $\mathcal{B}_q[v] \equiv \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt$.

Following the reviewer mandate, this note investigates the direct mathematical interaction between these terms and the Archimedean form $Q_{\mathrm{arch}}[v]$:
> **Can the exact translation-defect representation yield an analytical lower bound on $Q_{\mathrm{even}}[v]$ on $\Phi^\perp$?**
> Specifically, what is the direct relation between the weighted translation defects $\frac{1}{L} \sum w_q \mathcal{D}_q[v]$, the boundary mismatches $\frac{1}{L} \sum w_q \mathcal{B}_q[v]$, and the Archimedean form $Q_{\mathrm{arch}}[v]$, without introducing artificial intermediate scalar Sobolev functionals?

### Summary of Headline Analytical Findings
1. **Identical Vanishing of the Boundary Mismatch (Theorem 2.1):**  
   Every basis function in the even Connes–CvS Galerkin basis, $e_m(t) = \sqrt{2/L} \cos(2\pi m t / L)$, is **strictly symmetric under midpoint reflection** $t \mapsto L - t$. Consequently, for every test state $v \in \mathbb{R}^{N+1}$, $T_v(L - t) \equiv T_v(t)$, which proves algebraically that the boundary mismatch vanishes identically:
   $$\boxed{\mathcal{B}_q[v] \equiv \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt = 0 \quad (\forall q \in \mathcal{P}_c, \; \forall v \in \mathbb{R}^{N+1}).}$$
   The boundary mismatch term is identically zero across the entire Galerkin space, completely eliminating boundary mass asymmetry from the prime quadratic form.
2. **The Step-Potential Decomposition (Theorem 3.1):**  
   The prime form decomposes cleanly into an attractive step potential and a positive non-local kinetic translation defect:
   $$\langle v, Q_{\mathrm{prime}} v \rangle = - \frac{1}{L} \int_0^L W(t) T_v(t)^2 \, dt + \frac{1}{L} \sum_{q \le c} w_q \mathcal{D}_q[v],$$
   where $W(t) \equiv 2 \sum_{\log q \le t} w_q \ge 0$. The spatial average of this potential is an exact algebraic invariant:
   $$\frac{1}{L} \int_0^L W(t) \, dt \equiv -\psi_{\mathrm{prime}}'(0) = 2 \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \left(1 - \frac{\log q}{L}\right) \approx 2.9424 \quad (c = 13).$$
3. **The Effective Coupled Multiplier (Theorem 4.1):**  
   On pure Fourier modes $e_m$, the combined Archimedean and periodic translation-defect operator has multiplier:
   $$\Omega(m) \equiv h_+(a_m) + 4 M(m) = h_+(a_m) + 4 \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \sin^2\left(\frac{\pi m \log q}{L}\right).$$
   For $c = 13$, the translation defect supplies a massive positive injection on low modes:
   $$4 M(1) \approx +9.415, \qquad 4 M(2) \approx +9.940.$$
   This completely overpowers the negative Archimedean dip ($h_+(a_1) \approx -2.20$, $h_+(a_2) \approx -0.91$), yielding strictly positive combined multipliers:
   $$\Omega(1) \approx +7.215 > 0, \qquad \Omega(2) \approx +9.026 > 0.$$
4. **The Incommensurability Limit & The Two-Regime Target:**  
   Because the frequencies $\{\frac{\log p}{\log 13}\}$ are rationally independent, simultaneous Diophantine approximation implies that along suitable high-frequency subsequences $m_j \to \infty$, $M(m_j)$ can become arbitrarily small ($\inf_{m \ge 1} M(m) = 0$). Hence $4M(m) \approx 9.9$ is an empirical observation on low modes, not an asymptotic uniform floor.  
   However, as $m \to \infty$, the Archimedean multiplier grows logarithmically ($h_+(a_m) \sim \log m \to +\infty$). This indicates that coercivity should be proved not by asserting a uniform lower bound on $M(m)$, but via a **two-regime balance**:
   - At low $m$, $M(m)$ overcomes the negative Archimedean dip;
   - At high $m$, $h_+(a_m)$ itself supplies all required positivity.

---

## 2. Identical Vanishing of the Boundary Mismatch via Midpoint Symmetry

In [`cell127.md`](file:///c:/data/github/connes-cvs-/cell127.md) (Proposition 2.2), polarization isolated the boundary mismatch term:
$$\mathcal{B}_q[v] \equiv \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt.$$
We now prove that this term is identically zero for every test state in the Connes–CvS Galerkin truncation.

### Theorem 2.1 (Midpoint Symmetry & Identical Vanishing of $\mathcal{B}_q[v]$)
Let $T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos\left(\frac{2\pi m t}{L}\right)$ for $t \in [0, L]$.  
Then:
1. $T_v(t)$ is identically symmetric under midpoint reflection:
   $$T_v(L - t) \equiv T_v(t) \qquad (\forall t \in [0, L], \; \forall v \in \mathbb{R}^{N+1}).$$
2. Consequently, for every prime power $q \le c$:
   $$\boxed{\mathcal{B}_q[v] \equiv 0 \qquad (\forall v \in \mathbb{R}^{N+1}).}$$

*Proof.*  
1. Evaluate $T_v(L - t)$:
   $$T_v(L - t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos\left(\frac{2\pi m (L - t)}{L}\right) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos\left(2\pi m - \frac{2\pi m t}{L}\right).$$
   For every integer $m \in \mathbb{Z}$:
   $$\cos(2\pi m - \theta) = \cos(2\pi m) \cos\theta + \sin(2\pi m) \sin\theta = 1 \cdot \cos\theta + 0 = \cos\theta.$$
   Therefore, for each basis mode:
   $$\cos\left(\frac{2\pi m (L - t)}{L}\right) = \cos\left(\frac{2\pi m t}{L}\right).$$
   Summing over $m = 0, \dots, N$ yields $T_v(L - t) \equiv T_v(t)$ for all $t \in [0, L]$.
2. In the second integral of $\mathcal{B}_q[v]$, apply the change of variables $s = L - t$, so $ds = -dt$. When $t = L - \log q$, $s = \log q$; when $t = L$, $s = 0$:
   $$\int_{L - \log q}^L T_v(t)^2 \, dt = \int_{\log q}^0 T_v(L - s)^2 (-ds) = \int_0^{\log q} T_v(L - s)^2 \, ds.$$
   Using midpoint symmetry $T_v(L - s) = T_v(s)$:
   $$\int_{L - \log q}^L T_v(t)^2 \, dt = \int_0^{\log q} T_v(s)^2 \, ds.$$
   Substituting this into $\mathcal{B}_q[v]$ gives:
   $$\mathcal{B}_q[v] = \int_0^{\log q} T_v(t)^2 \, dt - \int_0^{\log q} T_v(s)^2 \, ds = 0.$$
   This holds identically for all $q \le c$ and all $v \in \mathbb{R}^{N+1}$. $\blacksquare$

*Significance:* The boundary mismatch term $\mathcal{B}_q[v]$ is not merely small or perturbative—it is an **exact algebraic zero** enforced by the parity structure of the Fourier cosine basis. There is zero boundary asymmetry leakage in the prime quadratic form.

---

## 3. The Step-Potential Decomposition & Spatial Moments

With $\mathcal{B}_q[v] = 0$, the prime quadratic form reduces to:
$$\langle v, Q_{\mathrm{prime}} v \rangle = -\frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L T_v(t)^2 \, dt + \frac{1}{L} \sum_{q \le c} w_q \mathcal{D}_q[v].$$

### Theorem 3.1 (Step-Potential Representation)
Define the step potential $W: [0, L] \to [0, \infty)$ by:
$$W(t) \equiv 2 \sum_{\substack{q \le c \\ \log q \le t}} \frac{\Lambda(q)}{\sqrt{q}} = 2 \sum_{q \le c} w_q \mathbf{1}_{[\log q, L]}(t).$$
Then:
$$\boxed{\langle v, Q_{\mathrm{prime}} v \rangle = - \frac{1}{L} \int_0^L W(t) T_v(t)^2 \, dt + \frac{1}{L} \sum_{q \le c} w_q \mathcal{D}_q[v].}$$

Furthermore:
1. **Total Weight:** $W(L) = 2 \sum_{q \le c} w_q \equiv 2 W$. For $c = 13$, $2 W \approx 9.9446$.
2. **Exact Spatial Mean:**
   $$\frac{1}{L} \int_0^L W(t) \, dt = -\psi_{\mathrm{prime}}'(0) = 2 \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \left(1 - \frac{\log q}{L}\right).$$
   For $c = 13$, $-\psi_{\mathrm{prime}}'(0) \approx 2.9424$.
3. **Exact Fourier Moments:** For every $k \ge 1$:
   $$\frac{1}{L} \int_0^L W(t) \cos\left(\frac{2\pi k t}{L}\right) \, dt = - \frac{1}{k} \psi_{\mathrm{prime}}(k).$$

*Proof.*  
1. Interchanging summation and integration:
   $$-\frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L T_v(t)^2 \, dt = - \frac{1}{L} \int_0^L \left( 2 \sum_{q \le c} w_q \mathbf{1}_{[\log q, L]}(t) \right) T_v(t)^2 \, dt = - \frac{1}{L} \int_0^L W(t) T_v(t)^2 \, dt.$$
2. Compute the mean:
   $$\frac{1}{L} \int_0^L W(t) \, dt = \frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L dt = \frac{2}{L} \sum_{q \le c} w_q (L - \log q) = 2 \sum_{q \le c} w_q \left(1 - \frac{\log q}{L}\right).$$
   From [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py) (eq. 141), $\psi_{\mathrm{prime}}'(0) = - 2 \sum w_q (1 - \log q / L)$, confirming the identity.
3. Compute the Fourier cosine moments for $k \ge 1$:
   $$\int_0^L W(t) \cos(a_k t) \, dt = 2 \sum_{q \le c} w_q \int_{\log q}^L \cos(a_k t) \, dt = 2 \sum_{q \le c} w_q \left[ \frac{\sin(a_k L) - \sin(a_k \log q)}{a_k} \right].$$
   Since $a_k L = 2\pi k$, $\sin(a_k L) = 0$. Using $\sin(a_k \log q) = \sin(2\pi k \log q / L) = -\sin(2\pi k (1 - \log q / L))$:
   $$\int_0^L W(t) \cos(a_k t) \, dt = - \frac{2}{a_k} \sum_{q \le c} w_q \sin(a_k \log q) = \frac{2\pi}{a_k} \left[ -\frac{1}{\pi} \sum_{q \le c} w_q \sin\left(2\pi k \left(1 - \frac{\log q}{L}\right)\right) \right].$$
   Recognizing the bracketed sum as $\psi_{\mathrm{prime}}(k)$ and substituting $a_k = 2\pi k / L$:
   $$\int_0^L W(t) \cos(a_k t) \, dt = \frac{2\pi}{2\pi k / L} \psi_{\mathrm{prime}}(k) = - \frac{L}{k} \psi_{\mathrm{prime}}(k).$$
   Dividing by $L$ gives the moment formula $-\frac{1}{k}\psi_{\mathrm{prime}}(k)$. $\blacksquare$

---

## 4. The Translation-Defect Multiplier $M(m)$ & Archimedean Coupling

We now analyze the positive restoring term:
$$\mathcal{T}_{\mathrm{defect}}[v] \equiv \frac{1}{L} \sum_{q \le c} w_q \mathcal{D}_q[v] = \frac{1}{L} \sum_{q \le c} w_q \int_{\log q}^L |T_v(t) - T_v(t - \log q)|^2 \, dt.$$

### 4.1 The Periodic Translation Defect
Extend $T_v$ periodically to $\mathbb{R}$ with period $L$. The periodic translation defect is:
$$\mathcal{D}_q^{\mathrm{per}}[v] \equiv \int_0^L |T_v(t) - T_v(t - \log q)|^2 \, dt.$$
Substituting the Fourier expansion $T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos(a_m t)$:
$$T_v(t) - T_v(t - \log q) = \sqrt{2} \sum_{m=1}^N v_m \left[ \cos(a_m t) - \cos(a_m (t - \log q)) \right] = -2\sqrt{2} \sum_{m=1}^N v_m \sin\left(\frac{a_m \log q}{2}\right) \sin\left(a_m t - \frac{a_m \log q}{2}\right).$$
Integrating the square over $[0, L]$ using orthogonality:
$$\mathcal{D}_q^{\mathrm{per}}[v] = 8 \sum_{m=1}^N v_m^2 \sin^2\left(\frac{a_m \log q}{2}\right) \int_0^L \sin^2\left(a_m t - \frac{a_m \log q}{2}\right) \, dt = 4 L \sum_{m=1}^N v_m^2 \sin^2\left(\frac{\pi m \log q}{L}\right).$$

### Theorem 4.1 (The Effective Coupled Multiplier)
Dividing by $L$ to match the coefficient quadratic form $\langle v, Q v \rangle$, the periodic translation defect form is diagonal in the Fourier basis:
$$\frac{1}{L} \sum_{q \le c} w_q \mathcal{D}_q^{\mathrm{per}}[v] = 4 \sum_{m=1}^N M(m) v_m^2,$$
where the arithmetic dispersion multiplier $M(m)$ is:
$$\boxed{M(m) \equiv \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \sin^2\left(\frac{\pi m \log q}{L}\right) \ge 0.}$$

Consequently, combining the diagonal Archimedean multiplier $h_+(a_m)$ with the translation defect yields the **effective coupled multiplier**:
$$\boxed{\Omega(m) \equiv h_+(a_m) + 4 M(m) = h_+(a_m) + 4 \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \sin^2\left(\frac{\pi m \log q}{L}\right).}$$

---

### 4.2 Empirical Multiplier Values on Low Modes for $c = 13$
For cutoff $c = 13$ and $L = \log 13 \approx 2.564949$, we evaluate $M(m)$ and $\Omega(m)$ across the critical low-to-moderate modes:

| Mode $m$ | $a_m = 2\pi m / L$ | $h_+(a_m)$ | $M(m)$ | $4 M(m)$ | $\Omega(m) = h_+(a_m) + 4 M(m)$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **$m = 1$** | $2.4497$ | $-2.2001$ | $2.3538$ | **$+9.4152$** | **$+7.2151$** |
| **$m = 2$** | $4.8994$ | $-0.9142$ | $2.4851$ | **$+9.9404$** | **$+9.0262$** |
| **$m = 3$** | $7.3491$ | $+0.3862$ | $2.4782$ | **$+9.9128$** | **$+10.2990$** |
| **$m = 4$** | $9.7988$ | $+0.7208$ | $2.4815$ | **$+9.9260$** | **$+10.6468$** |
| **$m = 6$** | $14.698$ | $+1.1412$ | $2.4862$ | **$+9.9448$** | **$+11.0860$** |
| **$m = 12$** | $29.396$ | $+1.5532$ | $2.4859$ | **$+9.9436$** | **$+11.4968$** |

Across these initial discrete modes, $4M(m) \approx +9.9$ injects substantial positive energy, turning the net multiplier strictly positive ($\ge +7.215$).

---

### 4.3 Diophantine Incommensurability & The Asymptotic Vanishing of $M(m)$
The observation that $4M(m) \approx 9.9$ on the tabulated modes is an empirical property of low frequencies, **not an asymptotic lower bound**.

Because the prime logarithms $\{\log 2, \log 3, \log 5, \dots\}$ are linearly independent over $\mathbb{Q}$, Kronecker's theorem and Dirichlet's simultaneous approximation theorem imply that the sequence of phase vectors:
$$\left( \left\{ m \frac{\log 2}{L} \right\}, \left\{ m \frac{\log 3}{L} \right\}, \dots, \left\{ m \frac{\log 11}{L} \right\} \right) \pmod 1$$
is dense in the multi-dimensional torus $\mathbb{T}^k$. Consequently, for any $\varepsilon > 0$, there exist arbitrarily large integers $m_j \to \infty$ such that all phases $m_j \log q / L$ are simultaneously within $\varepsilon$ of an integer. Along such a subsequence:
$$\sin^2\left(\frac{\pi m_j \log q}{L}\right) \le \pi^2 \varepsilon^2 \quad \Longrightarrow \quad M(m_j) \longrightarrow 0.$$
Therefore:
$$\boxed{\inf_{m \ge 1} M(m) = 0.}$$
Arithmetic incommensurability ensures that the translation defect $M(m)$ **cannot have a uniform positive infimum across all modes**.

However, this does not destroy coercivity:
- For high modes $m \ge M_0$, the Archimedean multiplier grows logarithmically:
  $$h_+(a_m) \sim \log\left(\frac{a_m}{2\pi}\right) = \log\left(\frac{m}{L}\right) \longrightarrow +\infty.$$
  Since $M(m) \ge 0$ unconditionally, $\Omega(m) \ge h_+(a_m) \ge c_0 > 0$ for all $m \ge M_0$.
- For low modes $1 \le m < M_0$, there are only finitely many modes to check, on which $M(m)$ is demonstrably bounded away from zero.
This establishes the **Two-Regime Multiplier Target**, which avoids relying on a false uniform infimum of $M(m)$.

---

## 5. The Truncation Defect & Finite-Interval Correction

We now evaluate the exact finite-interval correction:
$$\Delta\mathcal{D}_q[v] \equiv \mathcal{D}_q[v] - \mathcal{D}_q^{\mathrm{per}}[v] = - \int_0^{\log q} |T_v(t) - T_v(t - \log q)|^2 \, dt \le 0.$$

### Proposition 5.1 (Local Asymptotic Boundary Remainder)
1. The truncation error is strictly non-positive: $\Delta\mathcal{D}_q[v] \le 0$.
2. For any unit vector $v \in \mathbb{R}^{N+1}$:
   $$|\Delta\mathcal{D}_q[v]| \le 2 \int_0^{\log q} T_v(t)^2 \, dt + 2 \int_0^{\log q} T_v(t - \log q)^2 \, dt \le 4 \log q \, \|T_v\|_{L^\infty}^2.$$
3. For a fixed smooth profile with Dirichlet boundary vanishing ($T_v(0) = T_v(L) = 0$), local Taylor expansion yields:
   $$|\Delta\mathcal{D}_q[v]| \le \frac{2}{3} (\log q)^3 |T_v'(0)|^2 + \mathcal{O}((\log q)^4).$$
   *(Cautionary Epistemic Note: This Taylor estimate is a local asymptotic statement for smooth functions. Uniform control across the Galerkin sequence $v_N$ requires uniform boundary derivative bounds, which are governed by the non-circular regularity bridge rather than assumed here.)*

---

## 6. Synthesis: The Coupled Arithmetic-Archimedean Structure

We can now write the complete, exact expression for the full quadratic form $Q_{\mathrm{even}}[v]$ on the Galerkin space:

$$\boxed{\langle v, Q_{\mathrm{even}} v \rangle = \langle v, Q_{\mathrm{arch}} v \rangle + 4 \sum_{m=1}^N M(m) v_m^2 - \frac{1}{L} \int_0^L W(t) T_v(t)^2 \, dt + \frac{1}{L} \sum_{q \le c} w_q \Delta\mathcal{D}_q[v] + \langle v, Q_{\mathrm{pole}} v \rangle.}$$

### The Three Solid Structural Discoveries
1. **Identical Vanishing of Boundary Mismatch ($\mathcal{B}_q \equiv 0$):**  
   Midpoint reflection symmetry $T_v(L - t) \equiv T_v(t)$ eliminates boundary mass asymmetry identically across the entire Galerkin space.
2. **Step-Potential vs Kinetic Defect Representation:**  
   The prime form is not an opaque indefinite matrix, but an exact competition between an attractive local step potential $-W(t)$ and a positive non-local translation defect $\mathcal{D}_q$.
3. **Modewise Coupling with Archimedean Multiplier:**  
   On pure Fourier modes, the translation defect $4M(m)$ combines directly with the Archimedean multiplier $h_+(a_m)$, providing the missing positive energy to defeat the negative low-frequency Archimedean well.

### The Three Remaining Analytical Tasks for Gate 1
Cell 128 exposes the exact operator-theoretic architecture of the continuum margin, clarifying that three distinct analytical tasks remain to be completed:
1. **The Two-Regime Multiplier Theorem:** Rigorously establish $\inf_{m \ge 1}[h_+(a_m) + 4M(m)] \ge c_0 > 0$ by separating finite low-frequency checks ($1 \le m < m_0$) from high-frequency logarithmic Archimedean growth ($m \ge m_0$).
2. **Finite-Interval Defect Control:** Determine the minimal boundary and spectral conditions under which $\frac{1}{L}\sum w_q |\Delta\mathcal{D}_q[v]|$ remains bounded below the bulk positive margin.
3. **Step-Potential Operator Inequality:** Control the non-constant multiplication operator $\frac{1}{L}\int_0^L W(t) T_v(t)^2 dt$ relative to the positive diagonal multiplier, recognizing that $W(t)$ cannot be replaced simply by its spatial mean.

---

---

## References
- [`cell127.md`](file:///c:/data/github/connes-cvs-/cell127.md) — Coupled Regularity Functionals & Scaling Mismatch Obstruction
- [`cell126.md`](file:///c:/data/github/connes-cvs-/cell126.md) — Componentwise Form Domination & The Indefinite Prime Form Obstruction
- [`cell125.md`](file:///c:/data/github/connes-cvs-/cell125.md) — First-Principles Variational Quasimode Codimension Bound
- [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md) — Operator Decomposition Audit & Component Spectra
- [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py) — Galerkin Operator Assembly & Prime Symbol
- [`Paper-NR1.md`](file:///c:/data/github/connes-cvs-/Paper-NR1.md) — Exact Resolvent and Commutator Toolkit
- [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) — The Semiclassical Continuum Programme
- [`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md) — Strategic Roadmap (Gate 1 Pipeline)

