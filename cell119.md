# Cell 119 Analytical Note: Exact Algebraic Modal Decomposition ($A_{N, k}$ vs $B_{N, k}$ Cancellation Balance)

**Companion Computational Script:** [`cell119.py`](file:///c:/data/github/connes-cvs-/cell119.py) | **Verification Log:** `cell119.out` (pending compute node execution)  
**Status:** Pre-Flight Formulation & Analytical Architecture (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell118.md`](file:///c:/data/github/connes-cvs-/cell118.md); [`cell118.out`](file:///c:/data/github/connes-cvs-/cell118.out); [`cell117.md`](file:///c:/data/github/connes-cvs-/cell117.md); [`cell116.md`](file:///c:/data/github/connes-cvs-/cell116.md); [`cell112.md`](file:///c:/data/github/connes-cvs-/cell112.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Theoretical Context

### 1. The Post-Cell 118 Landscape
Cell 118 established three critical empirical facts regarding the finite-$N$ boundary flux $(H u_N)_N = \sum_{k=1}^N H_{Nk} k^2 v_{N, k}$:
1. **Low-Mode Peak at $k = 3$:** The cumulative modal flux partial sum $S_N(m) = \sum_{k=1}^m F_{N, k}$ universally peaks at $k = 3$ across all tested dimensions $N \in \{64, 96, 128, 160, 192\}$, reaching peak magnitudes $M_{\mathrm{peak}} \sim 10^{-3} - 10^{-2}$. However, the individual sign patterns $F_{N, k}$ vary with $N$.
2. **Early Extinction into a Quasi-Stable Residual:** The normalized master-curve collapse $S_N(x)/M_{\mathrm{peak}} \to G(x)$ is falsified. However, the *unnormalized* cumulative sums $S_N(x)$ converge remarkably rapidly: for $N = 192$, $S_N(0.1) \approx 4.99 \times 10^{-10}$, $S_N(0.2) \approx -3.38 \times 10^{-19}$, and by $x = 0.3$, $S_N(0.3) = -1.28 \times 10^{-22}$. Throughout $x \in [0.4, 0.95]$, $S_N(x)$ enters a quasi-stable plateau $\sim -1.3 \times 10^{-22}$, receiving only minor endpoint corrections at $x \approx 1$.
3. **Decades of Cancellation:** The dynamic range between the $k = 3$ excursion ($1.14 \times 10^{-2}$) and the final boundary flux ($-8.35 \times 10^{-23}$) spans **approximately 20 decades of cancellation** (a factor of $1.37 \times 10^{20}$).
4. **Failure of Simple Arithmetic Remainder Correlations:** Pearson correlations between boundary residuals and the simplest prime-power sums ($\Sigma_{\mathrm{cos}}(N), \psi_{\mathrm{prime}}'(N)$) are weak ($|r| \le 0.24$).

### 2. The Analytical Pivot: Why Does Cancellation Occur?
Having retired both the local Airy boundary-layer model (Cell 117) and the universal continuum curve collapse (Cell 118), the investigation turns from descriptive numerical reconnaissance to **exact algebraic analysis**.

The boundary flux is the endpoint value of a highly nonlocal cancellation across the modal spectrum. Cell 119 isolates the exact algebraic structure of this cancellation by splitting the divided-difference kernel into its two constitutive components:
$$F_{N, k} = A_{N, k} - B_{N, k},$$
where $A_{N, k}$ carries the boundary symbol $a_N = 2 N \psi(N)$ and $B_{N, k}$ carries the internal modal symbol $a_k = 2 k \psi(k)$.

---

## 1. Mathematical Derivations & Algebraic Structure

### 1.1 Exact Divided-Difference Decomposition
Recall that in the canonical even basis ($v$-basis of size $N+1$), the off-diagonal matrix elements of the Hamiltonian $H$ for $k \in \{1, \dots, N-1\}$ are:
$$H_{N, k} = \frac{2 N \psi(N) - 2 k \psi(k)}{N^2 - k^2} = \frac{a_N - a_k}{N^2 - k^2},$$
where $a_m \equiv 2 m \psi(m) = \sqrt{2} m^2 H_{0, m}$.

For each off-diagonal mode $k \in \{1, \dots, N-1\}$, the modal flux summand $F_{N, k} \equiv H_{N, k} k^2 v_{N, k}$ splits identically into:
$$\boxed{F_{N, k} = A_{N, k} - B_{N, k},}$$
where:
$$\boxed{A_{N, k} \equiv a_N \frac{k^2 v_{N, k}}{N^2 - k^2} = \frac{2 N \psi(N) k^2 v_{N, k}}{N^2 - k^2},}$$
$$\boxed{B_{N, k} \equiv a_k \frac{k^2 v_{N, k}}{N^2 - k^2} = \frac{2 k^3 \psi(k) v_{N, k}}{N^2 - k^2}.}$$

At the upper boundary endpoint $k = N$, the divided difference possesses a removable singularity, and the diagonal entry is given by:
$$F_{N, N} \equiv H_{N, N} N^2 v_{N, N} = \left(\psi'(N) + \frac{\psi(N)}{N}\right) N^2 v_{N, N}.$$

The total boundary flux is therefore:
$$(H u_N)_N = \sum_{k=1}^{N-1} F_{N, k} + F_{N, N} = \Sigma_A^{\mathrm{off}}(N) - \Sigma_B^{\mathrm{off}}(N) + F_{N, N},$$
where the total off-diagonal sums are:
$$\Sigma_A^{\mathrm{off}}(N) \equiv \sum_{k=1}^{N-1} A_{N, k}, \qquad \Sigma_B^{\mathrm{off}}(N) \equiv \sum_{k=1}^{N-1} B_{N, k}.$$

---

### 1.2 Low-Mode Asymptotic Expansion ($k \ll N$)
For low modes $k \ll N$, we expand the geometric denominator:
$$\frac{k^2}{N^2 - k^2} = \frac{k^2}{N^2} \frac{1}{1 - k^2/N^2} = \frac{k^2}{N^2} + \frac{k^4}{N^4} + \frac{k^6}{N^6} + \dots$$

Substituting this expansion into $A_{N, k}$ and $B_{N, k}$:
1. **The $A_{N, k}$ Piece:**
   $$A_{N, k} = 2 N \psi(N) \left[ \frac{k^2 v_{N, k}}{N^2} + \frac{k^4 v_{N, k}}{N^4} + \dots \right] = \frac{2\psi(N)}{N} k^2 v_{N, k} + \frac{2\psi(N)}{N^3} k^4 v_{N, k} + \dots$$
2. **The $B_{N, k}$ Piece:**
   $$B_{N, k} = 2 k \psi(k) \left[ \frac{k^2 v_{N, k}}{N^2} + \frac{k^4 v_{N, k}}{N^4} + \dots \right] = \frac{2 k^3 \psi(k) v_{N, k}}{N^2} + \frac{2 k^5 \psi(k) v_{N, k}}{N^4} + \dots$$

Notice the crucial relative scaling between $A_{N, k}$ and $B_{N, k}$ at low modes:
$$\frac{B_{N, k}}{A_{N, k}} = \frac{a_k}{a_N} = \frac{k \psi(k)}{N \psi(N)} \sim \mathcal{O}\left(\frac{k}{N}\right) \ll 1 \qquad (\text{for } k \ll N).$$

**Key Insight:**  
At low modes ($k \in \{1, 2, 3, \dots, 10\}$), $B_{N, k}$ is suppressed by an explicit factor of $k/N$ relative to $A_{N, k}$. Consequently:
$$F_{N, k} \approx A_{N, k} \approx \frac{2\psi(N)}{N} k^2 v_{N, k} \qquad (k \ll N).$$
This explains why the initial low-mode excursion is driven almost entirely by the boundary symbol $a_N$ acting on the kinetic weights $k^2 v_k$.

However, as $k$ approaches macroscopic fractions $x = k/N \sim 0.2 - 0.3$, $a_k = 2 k \psi(k)$ grows to become comparable in magnitude to $a_N = 2 N \psi(N)$.

---

### 1.3 Exact Connection to Theorem 1 and Boundary Contact
Recall Theorem 1 (from [`cell112.md`](file:///c:/data/github/connes-cvs-/cell112.md) and [`cell116.md`](file:///c:/data/github/connes-cvs-/cell116.md)):
$$(H u_N)_N = \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + E_1 N^2 v_{N, N}.$$

We now establish the exact algebraic link between the sums $\Sigma_A^{\mathrm{off}}, \Sigma_B^{\mathrm{off}}$ and the boundary terms in Theorem 1.

Using the exact rational identity $\frac{k^2}{N^2 - k^2} = -1 + \frac{N^2}{N^2 - k^2}$:
1. **For $\Sigma_A^{\mathrm{off}}$:**
   $$\Sigma_A^{\mathrm{off}}(N) = a_N \sum_{k=1}^{N-1} \frac{k^2 v_{N, k}}{N^2 - k^2} = -a_N \sum_{k=1}^{N-1} v_{N, k} + a_N N^2 \sum_{k=1}^{N-1} \frac{v_{N, k}}{N^2 - k^2}.$$
   Since $T_v(0) = v_{N, 0} + \sqrt{2} \sum_{k=1}^N v_{N, k}$, we have:
   $$\sum_{k=1}^{N-1} v_{N, k} = \frac{T_v(0) - v_{N, 0}}{\sqrt{2}} - v_{N, N}.$$
   Therefore:
   $$\boxed{\Sigma_A^{\mathrm{off}}(N) = -\frac{a_N}{\sqrt{2}} T_v(0) + \frac{a_N}{\sqrt{2}} v_{N, 0} + a_N v_{N, N} + a_N N^2 \sum_{k=1}^{N-1} \frac{v_{N, k}}{N^2 - k^2}.}$$

2. **For $\Sigma_B^{\mathrm{off}}$:**
   $$\Sigma_B^{\mathrm{off}}(N) = \sum_{k=1}^{N-1} a_k \frac{k^2 v_{N, k}}{N^2 - k^2} = -\sum_{k=1}^{N-1} a_k v_{N, k} + N^2 \sum_{k=1}^{N-1} \frac{a_k v_{N, k}}{N^2 - k^2}.$$
   Since $\alpha_N \equiv \sum_{k=1}^N a_k v_{N, k} = \sum_{k=1}^{N-1} a_k v_{N, k} + a_N v_{N, N}$, we have:
   $$\sum_{k=1}^{N-1} a_k v_{N, k} = \alpha_N - a_N v_{N, N}.$$
   Therefore:
   $$\boxed{\Sigma_B^{\mathrm{off}}(N) = -\alpha_N + a_N v_{N, N} + N^2 \sum_{k=1}^{N-1} \frac{a_k v_{N, k}}{N^2 - k^2}.}$$

3. **Subtracting $\Sigma_B^{\mathrm{off}}$ from $\Sigma_A^{\mathrm{off}}$:**
   $$\Sigma_A^{\mathrm{off}}(N) - \Sigma_B^{\mathrm{off}}(N) = \left( -\frac{a_N}{\sqrt{2}} T_v(0) + \frac{a_N}{\sqrt{2}} v_{N, 0} + a_N v_{N, N} + a_N N^2 \sum_{k=1}^{N-1} \frac{v_{N, k}}{N^2 - k^2} \right)$$
   $$- \left( -\alpha_N + a_N v_{N, N} + N^2 \sum_{k=1}^{N-1} \frac{a_k v_{N, k}}{N^2 - k^2} \right)$$
   $$= \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + \frac{a_N}{\sqrt{2}} v_{N, 0} + N^2 \sum_{k=1}^{N-1} \frac{a_N - a_k}{N^2 - k^2} v_{N, k}.$$

From the $N$-th row of the eigenvalue equation $(H v_N)_N = E_1 v_{N, N}$:
$$N^2 \sum_{k=1}^{N-1} \frac{a_N - a_k}{N^2 - k^2} v_{N, k} = N^2 \sum_{k=1}^{N-1} H_{Nk} v_{N, k} = E_1 N^2 v_{N, N} - \frac{a_N}{\sqrt{2}} v_{N, 0} - H_{NN} N^2 v_{N, N}.$$

The macroscopic cross-term $\frac{a_N}{\sqrt{2}} v_{N, 0} \sim \mathcal{O}(10)$ cancels **identically**:
$$\Sigma_A^{\mathrm{off}}(N) - \Sigma_B^{\mathrm{off}}(N) = \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + E_1 N^2 v_{N, N} - H_{NN} N^2 v_{N, N}.$$
Adding the diagonal term $F_{N, N} = H_{NN} N^2 v_{N, N}$ recovers Theorem 1 exactly:
$$(H u_N)_N = \Sigma_A^{\mathrm{off}}(N) - \Sigma_B^{\mathrm{off}}(N) + F_{N, N} = \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + E_1 N^2 v_{N, N}.$$

---

## 2. The Core Scientific Hypotheses for Cell 119

Cell 119 is designed to resolve which of the following mathematical mechanisms governs the 20 decades of cancellation:

### Hypothesis H1 (Massive Inter-Sum Cancellation)
*Both total off-diagonal sums $\Sigma_A^{\mathrm{off}}(N)$ and $\Sigma_B^{\mathrm{off}}(N)$ are individually macroscopic ($\sim 10^{-2}$ or $\mathcal{O}(1)$), and they cancel against each other to produce the $10^{-22}$ residual:*
$$|\Sigma_A^{\mathrm{off}}(N)| \gg 10^{-5}, \quad |\Sigma_B^{\mathrm{off}}(N)| \gg 10^{-5}, \qquad |\Sigma_A^{\mathrm{off}} - \Sigma_B^{\mathrm{off}}| \sim 10^{-22}.$$
*Falsification Criterion:* If either $|\Sigma_A^{\mathrm{off}}|$ or $|\Sigma_B^{\mathrm{off}}|$ is already of order $\le 10^{-15}$, H1 is falsified.

### Hypothesis H2 (Pre-Cancellation in the Discrete Moments)
*The individual sums $\Sigma_A^{\mathrm{off}}$ and $\Sigma_B^{\mathrm{off}}$ are already suppressed to tiny values because the solitary wave eigenvector $v_N$ satisfies discrete moment cancellation:*
$$M_2(N) \equiv \sum_{k=1}^{N-1} k^2 v_{N, k} \approx 0, \qquad M_{\psi, 3}(N) \equiv \sum_{k=1}^{N-1} 2 k^3 \psi(k) v_{N, k} \approx 0.$$
*In this case, the cancellation occurs internally within $\Sigma_A$ and within $\Sigma_B$ across modes, rather than between $\Sigma_A$ and $\Sigma_B$.*

### Hypothesis H3 (Trajectory Phase-Lock & Crossover)
*In the cumulative trajectory $S_A(x) = \sum_{k \le xN} A_{N, k}$ and $S_B(x) = \sum_{k \le xN} B_{N, k}$, $S_A(x)$ surges early at $x \le 0.05$ (generating the $k = 3$ peak), while $S_B(x)$ rises gradually as $x$ increases, exactly balancing $S_A(x)$ by the critical threshold $x^* \sim 0.2 - 0.3$.*

---

## 3. Module Layout in `cell119.py`

1. **Module 1: Low-Mode Term-by-Term Anatomy ($k \le 10$):**  
   Evaluates $A_{N, k}$, $B_{N, k}$, $F_{N, k} = A_{N, k} - B_{N, k}$, and the ratio $B_{N, k}/A_{N, k} = a_k/a_N$ across $N \in \{64, 96, 128, 160, 192\}$. Tests whether the $k = 3$ peak is driven purely by $A_{N, 3}$.
2. **Module 2: Total Sums Separation Across All 8 Dimensions:**  
   Evaluates $\Sigma_A^{\mathrm{off}}(N)$, $\Sigma_B^{\mathrm{off}}(N)$, their difference $\Delta \Sigma(N)$, the endpoint term $F_{N, N}$, and the cancellation factor $\mathcal{C}_{AB}(N) \equiv \max(|\Sigma_A|, |\Sigma_B|) / |\Delta \Sigma|$ across $N \in \{64, 80, 96, 112, 128, 144, 160, 192\}$.
3. **Module 3: Macroscopic Cumulative Trajectory Tracking ($x = k/N \in (0, 1)$):**  
   Samples $S_A(x)$, $S_B(x)$, and $S_F(x) = S_A(x) - S_B(x)$ on the continuum grid $x \in \{0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95\}$ for $N = 192$.
4. **Module 4: Moments Balance and Asymptotic Expansions:**  
   Evaluates the discrete solitary moments $M_2(N) = \sum k^2 v_k$, $M_4(N) = \sum k^4 v_k$, $M_{\psi, 3}(N) = 2 \sum k^3 \psi(k) v_k$, and $M_{\psi, 5}(N) = 2 \sum k^5 \psi(k) v_k$. Compares the leading asymptotic estimates $\Sigma_{A, \mathrm{lead}} = \frac{2\psi(N)}{N} M_2$ and $\Sigma_{B, \mathrm{lead}} = \frac{1}{N^2} M_{\psi, 3}$ against the full sums.
5. **Module 5: Exact Theorem 1 Synthesis & Residual Audit:**  
   Audits the exact boundary identity $(H u_N)_N = \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + E_1 N^2 v_{N, N}$, evaluating the individual contact and resolvent pieces to verify identity residuals to $> 65$ decimal digits.

---

## 4. Expected Impact on Gate 1

If **Hypothesis H1** holds, the $10^{20}$-scale cancellation is an inter-symbol cancellation between the global boundary weight $a_N$ and the distributed modal symbol $a_k$. This would reduce Gate 1 boundary extinction to an exact operator-integral duality in the continuum limit.

If **Hypothesis H2** holds, the cancellation is intrinsic to the solitary wave profile itself (curvature and higher moment vanishing at the boundary), establishing that boundary extinction is driven by the internal boundary smoothness of the solitary wave.
