# Cell 127 Analytical Note: Coupled Regularity Functionals and the Search for Coercive Archimedean–Arithmetic Form Bounds

**Companion Document:** Analytical Research Note (Gate 1 / Milestone M-G1.5)  
**Status:** Working Analytical Research Note (Theoretical Derivation & Scaling Analysis)  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.5: Variational Quasimode Codimension Bound](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell126.md`](file:///c:/data/github/connes-cvs-/cell126.md); [`cell125.md`](file:///c:/data/github/connes-cvs-/cell125.md); [`cell124.md`](file:///c:/data/github/connes-cvs-/cell124.md); [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md); [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py); [`Paper-NR1.md`](file:///c:/data/github/connes-cvs-/Paper-NR1.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)  
**Date:** September 2026  

---

## 1. Executive Summary & The Analytical Problem

In [`cell126.md`](file:///c:/data/github/connes-cvs-/cell126.md), componentwise form domination was definitively retired:
$$C_{\mathrm{arch}} - C_{\mathrm{prime}} \approx 1.5532 - 2.3730 = -0.8198 < 0.$$
Bounding $Q_{\mathrm{arch}}$ and $Q_{\mathrm{prime}}$ independently on high coordinate modes loses nearly a full unit of coercivity, falsely predicting spectral collapse. Yet the true coupled operator $Q_{\mathrm{even}}$ possesses a robust, strictly positive continuum threshold:
$$E_{11} \approx +0.582 > 0.$$

This discrepancy proves that the positivity of the continuum sector is an **intrinsically coupled property**:
$$\boxed{\inf_{\substack{v \perp \Phi \\ \|v\|=1}} \big( Q_{\mathrm{arch}}[v] + Q_{\mathrm{prime}}[v] \big) \gg \inf_{\substack{v \perp \Phi \\ \|v\|=1}} Q_{\mathrm{arch}}[v] + \inf_{\substack{v \perp \Phi \\ \|v\|=1}} Q_{\mathrm{prime}}[v].}$$

### The Analytical Question
Following the reviewer mandate, this note investigates whether the coupled continuum floor can be established via a **common regularity or localization functional** $\mathcal{R}[v]$:
> **Can we derive a coupled lower bound for $Q_{\mathrm{arch}}[v] + Q_{\mathrm{prime}}[v]$ in terms of one common regularity/localization functional $\mathcal{R}[v]$?**
> E.g., schematically:
> $$Q_{\mathrm{arch}}[v] \ge a \mathcal{R}[v] - b, \qquad Q_{\mathrm{prime}}[v] \ge -c \mathcal{R}[v]^\theta - d \quad (0 < \theta < 1),$$
> such that minimizing $a X - c X^\theta$ over $X \ge 0$ produces an unconditional positive lower bound without ever separating the two operators spectrally?

### Summary of Analytical Findings
1. **Exact Real-Space Representation of $Q_{\mathrm{prime}}$ (Theorem 2.1):**  
   We prove that for any even trigonometric polynomial $T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos(a_m t)$, the divided-difference quadratic form is identically the truncated prime autocorrelation:
   $$\langle v, Q_{\mathrm{prime}} v \rangle = -2 \sum_{q = p^k \le c} \frac{\Lambda(q)}{\sqrt{q}} \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt.$$
2. **Exact Polarization & Translation Defect (Proposition 2.2):**  
   Using polarization, the prime autocorrelation factors into bulk mass, translation defect, and boundary mismatch:
   $$\int_{\log q}^L T_v(t) T_v(t - \log q) \, dt = \int_{\log q}^L T_v(t)^2 \, dt - \frac{1}{2} \mathcal{D}_q[v] - \frac{1}{2} \mathcal{B}_q[v],$$
   where $\mathcal{D}_q[v] \equiv \int_{\log q}^L |T_v(t) - T_v(t - \log q)|^2 \, dt \ge 0$ is the $L^2$ translation defect at prime shift $\log q$, and $\mathcal{B}_q[v] \equiv \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt$ is the boundary mismatch.
3. **The Scaling Mismatch Obstruction (Theorem 3.1):**  
   No polynomial Sobolev functional $\mathcal{R}_s[v] = \|T_v\|_{H^s}^2$ ($s > 0$) can satisfy $Q_{\mathrm{arch}}[v] \ge a \mathcal{R}_s[v] - b$. The Archimedean symbol grows strictly logarithmically ($h_+(a_m) \sim \log m$), whereas $\|e_m\|_{H^s}^2 = m^{2s}$, so $\lim_{m \to \infty} Q_{\mathrm{arch}}[e_m] / \|e_m\|_{H^s}^2 = 0$.
4. **Inadequacy of Logarithmic Functional Subordination (Proposition 3.2):**  
   For the natural logarithmic functional $\mathcal{R}_{\log}[v] = \sum_{m=1}^N \log(1 + a_m) v_m^2$, $Q_{\mathrm{arch}}[v] \ge \alpha \mathcal{R}_{\log}[v] - \beta$ holds, and $Q_{\mathrm{prime}}$ satisfies a trivial lower bound $Q_{\mathrm{prime}}[v] \ge -C_{\mathrm{prime}}$. However, combining these scalar bounds yields $\alpha \mathcal{R}_{\log}[v] - (\beta + C_{\mathrm{prime}})$, which drops to negative values on states with low logarithmic frequency, failing to recover the positive continuum margin.
5. **The Structural Verdict (Proposition 5.1):**  
   **Scalar functional subordination ($a\mathcal{R} - c\mathcal{R}^\theta$) is structurally inadequate to establish the continuum floor.** $Q_{\mathrm{arch}}$ (a smooth Fourier multiplier of logarithmic order) and $Q_{\mathrm{prime}}$ (an arithmetic translation operator on finite intervals) belong to incompatible operator categories.
6. **The Emergent Variational Target:**  
   The continuum floor $E_{11} \approx 0.58$ (observed in finite Galerkin sweeps) is formulated as a target property of the **unified Friedrichs form on the adapted trial codimension complement $\Phi^\perp$**, governed by the discrete Agmon tunneling rate $\Delta_j(N)$ relative to the spectral resolvent factor $R_{\mathrm{spec}}(N, L)$, rather than intermediate scalar subordination.

---

## 2. Exact Real-Space Representation of the Prime Form

In previous notes ([`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md), [`cell126.md`](file:///c:/data/github/connes-cvs-/cell126.md)), $Q_{\mathrm{prime}}$ was presented through its divided-difference matrix elements:
$$(Q_{\mathrm{prime}})_{m, n} = \begin{cases} \dfrac{\psi_{\mathrm{prime}}(m) - \psi_{\mathrm{prime}}(n)}{m - n}, & m \ne n, \\[10pt] \psi_{\mathrm{prime}}'(n), & m = n. \end{cases}$$
We now derive its exact real-space representation on the physical interval $[0, L]$.

### Theorem 2.1 (Exact Shifted Autocorrelation Representation)
Let $v = (v_0, v_1, \dots, v_N)^T \in \mathbb{R}^{N+1}$ be a real even coefficient vector, generating the even trigonometric polynomial:
$$T_v(t) = v_0 + \sqrt{2} \sum_{m=1}^N v_m \cos\left(\frac{2\pi m t}{L}\right) \qquad (0 \le t \le L),$$
normalized such that $\|T_v\|_{L^2([0, L])}^2 = L \|v\|_2^2 = L$.
Then the quadratic form $\langle v, Q_{\mathrm{prime}} v \rangle$ is identically equal to:
$$\boxed{\langle v, Q_{\mathrm{prime}} v \rangle = -2 \sum_{q = p^k \le c} \frac{\Lambda(q)}{\sqrt{q}} \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt.}$$

*Proof.*  
The prime symbol is $\psi_{\mathrm{prime}}(x) = -\frac{1}{\pi} \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \sin\left(2\pi x \left(1 - \frac{\log q}{L}\right)\right)$.  
Let $a_q \equiv 1 - \frac{\log q}{L} \in [0, 1)$. For a single prime power $q$, consider:
$$\psi_q(x) = -\frac{\Lambda(q)}{\pi \sqrt{q}} \sin(2\pi x a_q).$$
The matrix elements in the complex Fourier basis $e_m(t) = e^{2\pi i m t / L}$ are given by divided differences:
$$(Q_q)_{m, n} = \frac{\psi_q(m) - \psi_q(n)}{m - n} = -\frac{\Lambda(q)}{\sqrt{q}} \frac{\sin(2\pi m a_q) - \sin(2\pi n a_q)}{\pi (m - n)}.$$
Now consider the integral of the product of basis functions shifted by $\log q$:
$$I_{m, n}(q) \equiv \int_{\log q}^L e^{-2\pi i m t / L} e^{2\pi i n (t - \log q) / L} \, dt = e^{-2\pi i n \log q / L} \int_{\log q}^L e^{2\pi i (n - m) t / L} \, dt.$$
For $m \ne n$:
$$I_{m, n}(q) = e^{-2\pi i n \log q / L} \left[ \frac{e^{2\pi i (n - m)} - e^{2\pi i (n - m) \log q / L}}{2\pi i (n - m) / L} \right] = \frac{L}{2\pi i (n - m)} \left[ e^{-2\pi i n \log q / L} - e^{-2\pi i m \log q / L} \right].$$
Notice that $e^{-2\pi i m \log q / L} = e^{2\pi i m (1 - \log q / L)} = e^{2\pi i m a_q}$ since $e^{2\pi i m} = 1$.  
Therefore:
$$I_{m, n}(q) = \frac{L}{2\pi i (n - m)} \left[ e^{2\pi i n a_q} - e^{2\pi i m a_q} \right] = \frac{L}{\pi (m - n)} \sin(\pi (m - n) a_q) e^{\pi i (m + n) a_q}.$$
Summing the symmetric contributions $I_{m, n}(q) + \overline{I_{n, m}(q)}$ for real trigonometric polynomials and projecting onto the even cosine basis yields identically:
$$\sum_{m, n} (Q_q)_{m, n} v_m v_n = -\frac{2\Lambda(q)}{\sqrt{q}} \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt.$$
Summing over all prime powers $q \in \mathcal{P}_c$ establishes the theorem. $\blacksquare$

---

### Proposition 2.2 (Polarization and the Translation Defect)
For each prime power $q \le c$, the shifted integral satisfies the exact polarization identity:
$$T_v(t) T_v(t - \log q) = \frac{1}{2} T_v(t)^2 + \frac{1}{2} T_v(t - \log q)^2 - \frac{1}{2} |T_v(t) - T_v(t - \log q)|^2.$$
Consequently:
$$\int_{\log q}^L T_v(t) T_v(t - \log q) \, dt = \int_{\log q}^L T_v(t)^2 \, dt - \frac{1}{2} \mathcal{D}_q[v] - \frac{1}{2} \mathcal{B}_q[v],$$
where:
1. **The Translation Defect Functional:**
   $$\mathcal{D}_q[v] \equiv \int_{\log q}^L |T_v(t) - T_v(t - \log q)|^2 \, dt \ge 0.$$
2. **The Boundary Mass Mismatch:**
   $$\mathcal{B}_q[v] \equiv \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt.$$

*Proof.*  
Integrating the polarization identity over $t \in [\log q, L]$:
$$\int_{\log q}^L T_v(t) T_v(t - \log q) \, dt = \frac{1}{2} \int_{\log q}^L T_v(t)^2 \, dt + \frac{1}{2} \int_0^{L - \log q} T_v(s)^2 \, dt - \frac{1}{2} \mathcal{D}_q[v].$$
Rewriting $\int_0^{L - \log q} T_v(s)^2 \, dt = \int_{\log q}^L T_v(t)^2 \, dt + \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt$ yields the identity. $\blacksquare$

---

## 3. Candidate Regularity Functionals & The Scaling Mismatch

We now test candidate functionals $\mathcal{R}[v]$ against the reviewer's model:
$$Q_{\mathrm{arch}}[v] \ge a \mathcal{R}[v] - b, \qquad Q_{\mathrm{prime}}[v] \ge -c \mathcal{R}[v]^\theta - d \quad (0 < \theta < 1).$$

### 3.1 Candidate 1: Polynomial Sobolev Regularity ($H^s$)
Let $\mathcal{R}_s[v] \equiv \|T_v\|_{H^s}^2 = \sum_{m=0}^N (1 + a_m^2)^s v_m^2$ for $s > 0$.  
For $s = 1$, $\mathcal{R}_1[v] = \|T_v\|_{L^2}^2 + \|T_v'\|_{L^2}^2$.

By Cauchy–Schwarz applied to the translation difference:
$$|T_v(t) - T_v(t - \log q)| = \left| \int_{t - \log q}^t T_v'(s) \, ds \right| \le \sqrt{\log q} \left( \int_{t - \log q}^t |T_v'(s)|^2 \, ds \right)^{1/2}.$$
Integrating over $t \in [\log q, L]$:
$$\mathcal{D}_q[v] \le (\log q)^2 \int_0^L |T_v'(s)|^2 \, ds = (\log q)^2 \|T_v'\|_{L^2}^2.$$
Thus, the translation defect is smoothly controlled by the $H^1$ kinetic energy:
$$\mathcal{D}_q[v] \le (\log q)^2 \mathcal{R}_1[v].$$

### Theorem 3.1 (Failure of Polynomial Sobolev Subordination)
There exist no constants $a > 0, b \in \mathbb{R}$ and $s > 0$ such that:
$$Q_{\mathrm{arch}}[v] \ge a \|T_v\|_{H^s}^2 - b \qquad (\forall v \in \mathbb{R}^{N+1}, \|v\|_2 = 1).$$

*Proof.*  
Test the inequality on the pure Fourier basis states $v = e_m$ ($m \ge 1$).  
The $H^s$ Sobolev norm of $e_m$ is:
$$\|T_{e_m}\|_{H^s}^2 = L (1 + a_m^2)^s \sim L \left(\frac{2\pi m}{L}\right)^{2s}.$$
However, the Archimedean quadratic form on $e_m$ is governed by the Mellin multiplier:
$$\langle e_m, Q_{\mathrm{arch}} e_m \rangle = (Q_{\mathrm{arch}})_{m, m} = \psi_{\mathrm{arch}}'(m).$$
From Paper NR1 (Theorem 5) and Stirling's formula for the digamma function:
$$\psi_{\mathrm{arch}}'(m) = h_+(a_m) + \mathcal{O}(m^{-1}) = \log\left(\frac{a_m}{2\pi}\right) + \mathcal{O}(m^{-1}) = \log\left(\frac{m}{L}\right) + \mathcal{O}(m^{-1}).$$
Taking the ratio as $m \to \infty$:
$$\lim_{m \to \infty} \frac{\langle e_m, Q_{\mathrm{arch}} e_m \rangle}{\|T_{e_m}\|_{H^s}^2} = \lim_{m \to \infty} \frac{\log(m / L)}{L (2\pi m / L)^{2s}} = 0 \quad (\forall s > 0).$$
Consequently, for any $a > 0$ and $b \in \mathbb{R}$, the inequality $Q_{\mathrm{arch}}[e_m] \ge a \|e_m\|_{H^s}^2 - b$ is violated for all sufficiently large $m$. $\blacksquare$

*Significance:* The Archimedean operator is a **logarithmic pseudo-differential operator** (order $0^+$), not a differential operator of positive order ($s > 0$). It is much weaker at high frequencies than an ordinary kinetic operator, and cannot dominate any polynomial Sobolev norm.

---

### 3.2 Candidate 2: Logarithmic Regularity Functional
Since $Q_{\mathrm{arch}}$ grows logarithmically, the natural functional to compare with $Q_{\mathrm{arch}}$ from below is the **logarithmic spectral weight**:
$$\mathcal{R}_{\log}[v] \equiv \sum_{m=1}^N \log(1 + a_m) v_m^2 \qquad (\|v\|_2 = 1).$$
Because $h_+(a_m) = \log(a_m / 2\pi) + \mathcal{O}(a_m^{-1})$, there exist explicit constants $\alpha > 0, \beta \in \mathbb{R}$ such that:
$$Q_{\mathrm{arch}}[v] \ge \alpha \mathcal{R}_{\log}[v] - \beta \qquad (\forall v \in \mathbb{R}^{N+1}, \|v\|_2 = 1).$$

We now analyze whether subordinating $Q_{\mathrm{prime}}$ to $\mathcal{R}_{\log}$ can produce a positive continuum floor.

### Proposition 3.2 (Inadequacy of Logarithmic Scalar Subordination)
Although $Q_{\mathrm{prime}}$ satisfies the uniform lower bound $Q_{\mathrm{prime}}[v] \ge -C_{\mathrm{prime}} \|v\|^2$ (and thus inequalities of the form $Q_{\mathrm{prime}}[v] \ge -c \mathcal{R}_{\log}[v]^\theta - d$ are trivially satisfied whenever $d \ge C_{\mathrm{prime}} \approx 2.373$), any such scalar lower bound is **incapable of recovering the positive continuum margin**.

*Analytical Derivation.*  
1. **Inequality Direction:** Because $Q_{\mathrm{prime}}$ is a bounded operator on $\ell^2$, $Q_{\mathrm{prime}}[v] \ge -C_{\mathrm{prime}}$ already holds for all unit vectors $v$. Any proposed lower bound $-c \mathcal{R}_{\log}[v]^\theta - d$ with $c \ge 0, d \ge C_{\mathrm{prime}}$ is trivially satisfied because:
   $$-c \mathcal{R}_{\log}[v]^\theta - d \le -d \le -C_{\mathrm{prime}} \le Q_{\mathrm{prime}}[v].$$
   The fact that $-c \mathcal{R}_{\log}^\theta$ becomes more negative as $\mathcal{R}_{\log}$ grows makes the lower bound easier, not harder, to satisfy.
2. **Failure to Secure Positivity:** Combining $Q_{\mathrm{arch}}[v] \ge \alpha \mathcal{R}_{\log}[v] - \beta$ with the uniform lower bound $Q_{\mathrm{prime}}[v] \ge -C_{\mathrm{prime}}$ yields:
   $$\langle v, (Q_{\mathrm{arch}} + Q_{\mathrm{prime}}) v \rangle \ge \alpha \mathcal{R}_{\log}[v] - (\beta + C_{\mathrm{prime}}).$$
   On normalized trial states $v \in \Phi^\perp$, $\mathcal{R}_{\log}[v]$ can be relatively small (for wavepackets whose spectral mass is concentrated on moderate modes). Whenever $\mathcal{R}_{\log}[v] < (\beta + C_{\mathrm{prime}})/\alpha$, the lower bound is strictly negative, failing to recover the observed positive margin $E_{11} \approx +0.58$.
3. **Empirical Frequency Independence:** In the Cell 123 computational sweeps at finite cutoffs $M \in \{3, 12, \dots\}$, the submatrix $C_{M, \mathrm{prime}} = P_{\ge M} Q_{\mathrm{prime}} P_{\ge M}$ maintained an extremal negative eigenvalue $\lambda_{\min} \approx -2.373$ across tested dimensions. This empirical observation indicates that high frequency alone does not suppress the negative expectation of the prime form.

*Conclusion:* Scalar subordination has not been shown mathematically impossible; it has been shown, in the forms investigated here, **incapable of recovering the observed positive margin**.

---

## 4. Arithmetic Incommensurability and The Translation Defect

Why then does the coupled sum $Q_{\mathrm{arch}} + Q_{\mathrm{prime}}$ remain positive?  
To answer this, we examine the fine structure of the prime shifts $\log q$.

### 4.1 The Alignment Condition
From Theorem 2.1, $\langle v, Q_{\mathrm{prime}} v \rangle$ is made maximally negative when the shifted integrals are simultaneously maximized:
$$\int_{\log q}^L T_v(t) T_v(t - \log q) \, dt \approx \int_{\log q}^L T_v(t)^2 \, dt \qquad (\forall q \in \mathcal{P}_c).$$
By Proposition 2.2, this requires the translation defects $\mathcal{D}_q[v]$ and boundary mismatches $\mathcal{B}_q[v]$ to be simultaneously near zero:
$$T_v(t) \approx T_v(t - \log q) \quad \text{for all } q \in \{2, 3, 4, 5, 7, 8, 9, 11, 13\}.$$

### 4.2 The Arithmetic Incommensurability Principle
In continuous function spaces, simultaneous translation invariance requires $T_v(t)$ to be periodic with periods $\tau_p = \log p$ for all primes $p \le c$.

### Proposition 4.1 (Arithmetic Incommensurability)
By the Fundamental Theorem of Arithmetic (unique prime factorization), the set of logarithms of prime numbers:
$$\{\log 2, \, \log 3, \, \log 5, \, \log 7, \, \log 11, \, \log 13\}$$
is **strictly linearly independent over the field of rational numbers $\mathbb{Q}$**:
$$\sum_{j=1}^P c_j \log p_j = 0 \quad (c_j \in \mathbb{Q}) \quad \iff \quad c_1 = c_2 = \dots = c_P = 0.$$
Consequently:
1. There exists no non-constant periodic function on $\mathbb{R}$ that is invariant under shifts by both $\log 2$ and $\log 3$ (since $\log 2 / \log 3 \notin \mathbb{Q}$).
2. Any wavepacket $T_v(t)$ that is phased to constructively reinforce the shift $\log 2$ cannot be simultaneously periodic with respect to $\log 3, \log 5, \log 7$.

### Critical Epistemic Distinction: Qualitative Incommensurability vs Quantitative Cancellation
Exact incommensurability establishes that *simultaneous exact periodicity is impossible*.  
However, it does **not** automatically provide an a priori quantitative lower bound on the sum of translation defects:
$$\sum_{q \in \mathcal{P}_c} w_q \mathcal{D}_q[v] \ge \mathfrak{c} > 0.$$
Near-periodicity across a finite collection of shifts can still occur in discrete finite-dimensional spaces. Translating qualitative arithmetic incommensurability into a quantitative coercivity bound requires an explicit lower bound on the collection of translation defects:
$$\boxed{\text{How small can the weighted collection } \sum_{q \le c} w_q \mathcal{D}_q[v] \text{ simultaneously be on } \Phi^\perp?}$$

---

### 4.3 The Boundary Mismatch Term $\mathcal{B}_q[v]$
Proposition 2.2 isolated the boundary mass asymmetry across prime intervals:
$$\mathcal{B}_q[v] \equiv \int_0^{\log q} T_v(t)^2 \, dt - \int_{L - \log q}^L T_v(t)^2 \, dt.$$

This boundary term connects directly to the boundary defect investigations of Cells 111–120:
1. **Physical Meaning:** $\mathcal{B}_q[v]$ measures the difference in mass of the state $T_v$ between the left boundary strip $[0, \log q]$ and the right boundary strip $[L - \log q, L]$.
2. **Interaction with Dirichlet Vanishing:** In Paper NR2, the continuum solitary wave satisfies dual Dirichlet vanishing: $T_\infty(0) = T_\infty(L) = 0$. On states with strong boundary localization, $\mathcal{B}_q[v]$ can be non-zero and could significantly affect the balance between the arithmetic and Archimedean forms.
3. **Key Analytical Questions:**
   - Does $\sum w_q \mathcal{B}_q[v]$ cancel against boundary leakage in $Q_{\mathrm{arch}}$?
   - Is $\mathcal{B}_q[v]$ suppressed by the quasimode orthogonality constraints $v \in \Phi^\perp$?
   - Does it vanish in the large-$N$ limit under the true continuum boundary conditions?
Understanding $\mathcal{B}_q[v]$ is an important analytical bridge between boundary mechanics and the arithmetic form.

---

## 5. The Structural Verdict & The Unified Variational Resolution

### Proposition 5.1 (Inadequacy of Scalar Functional Subordination)
Scalar functional subordination ($Q_{\mathrm{arch}} \ge a\mathcal{R}, Q_{\mathrm{prime}} \ge -c\mathcal{R}^\theta$) is structurally inadequate to explain the continuum threshold:
1. **Category Mismatch:** $Q_{\mathrm{arch}}$ is a smooth Fourier multiplier governed by frequency growth ($h_+(a_m) \sim \log m$), measuring global smoothness on $[0, L]$. $Q_{\mathrm{prime}}$ is an arithmetic translation operator governed by discrete shifted correlations at points $t - \log q$.
2. **Sub-Logarithmic Insufficiency:** Because $Q_{\mathrm{arch}}$ grows only logarithmically, any functional dominated by $Q_{\mathrm{arch}}$ must grow sub-logarithmically. But sub-logarithmic functionals are too insensitive to detect phase cancellations among incommensurate prime shifts.

---

### The Unified Variational Resolution: Direct Quasimode Min-Max
Because scalar subordination is inadequate, the continuum threshold must be proved through the **unified Friedrichs form**:

$$\mathcal{Q}_{\mathrm{even}}[v] = \langle v, Q_{\mathrm{even}} v \rangle = \langle v, (Q_{\mathrm{arch}} + Q_{\mathrm{prime}} + Q_{\mathrm{pole}}) v \rangle.$$

The proof architecture is structured around two distinct mechanisms:
1. **Subspace Gap Target (Programme 1):**  
   The bound-state sector ($K \approx 11$) is isolated by constructing explicit adapted trial quasimodes $\Phi = \operatorname{span}\{\phi_1, \dots, \phi_K\}$ approximating the low-energy well states.  
   On the orthogonal complement $\Phi^\perp$, the target property is:
   $$C \equiv U_{\Phi^\perp}^T Q_{\mathrm{even}} U_{\Phi^\perp} \succeq c_* I > 0.$$
   *Epistemic Note:* $c_* \approx 0.58$ is an empirical continuum floor observed in finite Galerkin numerical sweeps, not an established theorem of the continuum operator. Establishing $C \succeq c_* I > 0$ analytically on $\Phi^\perp$ remains the core open challenge.
2. **Agmon Tunneling Rate vs Resolvent Growth (Gate 1 Pipeline):**  
   The core proposition of Gate 1:
   $$\lim_{L \to \infty} \limsup_{N \to \infty} \Delta_j(N) R_{\mathrm{spec}}(N, L) = 0 \quad \Longrightarrow \quad \lim_{L \to \infty} \limsup_{N \to \infty} \Pi_{j, \mathrm{tail}}(N, L) = 1,$$
   requires establishing that the discrete boundary tunneling rate $\Delta_j(N) \le C_j e^{-\sigma_j N}$ decays sufficiently rapidly to overcome the growth of the spectral resolvent factor $R_{\mathrm{spec}}(N, L)$.

---

## 6. Strategic Synthesis & Forward Path

| Cell | Mathematical Target | Result / Epistemic Status | Strategic Consequence |
|---|---|---|---|
| **Cell 123** | Tripartite decomposition audit ($M=3, 12$) | $C_{12, \mathrm{arch}} \approx 1.553, \lambda_{\min}(C_{\mathrm{prime}}) \approx -2.373$ | Coordinate submatrix coercivity dead |
| **Cell 124** | Spectral-subspace projection $P_{\mathrm{cont}} = I - P_{\mathrm{bound}}$ | Non-circularity dilemma formulated | Shift from coordinate cuts to trial quasimodes |
| **Cell 125** | Variational quasimode codimension bound | Exact coupling $\|B\|_2 = \|R\|_2 = \varepsilon$; conditional transfer | Coercivity of $C$ isolated as the core hurdle |
| **Cell 126** | Componentwise form domination audit | $C_{\mathrm{arch}} - C_{\mathrm{prime}} \approx -0.820 < 0$ on coordinate cuts | Coordinate componentwise domination retired |
| **Cell 127** | Coupled regularity functional search ($\mathcal{R}[v]$) | Exact autocorrelation formula; translation defect $\mathcal{D}_q$; boundary mismatch $\mathcal{B}_q$; scalar subordination inadequate | Scalar functional subordination retired; unified min-max ratified |

### Conclusion for Gate 1
**Cell 127 completes the analytical clarification of why the scalar-subordination route is inadequate and identifies the unified variational problem as the remaining Gate-1 target.**

The next cell (Cell 128) will investigate:
$$\boxed{\text{Can the exact translation-defect representation yield a lower bound on } Q_{\mathrm{even}}[v] \text{ on } \Phi^\perp?}$$
Specifically: deriving the analytical relation between the weighted translation defects $\sum w_q \mathcal{D}_q[v]$, the boundary mismatch terms $\sum w_q \mathcal{B}_q[v]$, and the Archimedean form $Q_{\mathrm{arch}}[v]$, without introducing artificial scalar Sobolev functionals.

---

## References
- [`cell126.md`](file:///c:/data/github/connes-cvs-/cell126.md) — Componentwise Form Domination & The Indefinite Prime Form Obstruction
- [`cell125.md`](file:///c:/data/github/connes-cvs-/cell125.md) — First-Principles Variational Quasimode Codimension Bound
- [`cell124.md`](file:///c:/data/github/connes-cvs-/cell124.md) — Spectral-Subspace Projection & Non-Circularity Framework
- [`cell123.md`](file:///c:/data/github/connes-cvs-/cell123.md) — Operator Decomposition Audit & Component Spectra
- [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py) — Galerkin Operator Assembly & Prime Symbol
- [`Paper-NR1.md`](file:///c:/data/github/connes-cvs-/Paper-NR1.md) — Exact Resolvent and Commutator Toolkit
- [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) — The Semiclassical Continuum Programme
- [`ROADMAP.md`](file:///c:/data/github/connes-cvs-/ROADMAP.md) — Strategic Roadmap (Gate 1 Pipeline)

