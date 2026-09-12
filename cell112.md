# Cell 112 Analytical Note: High-Precision Extinction Audit & Multi-Route Boundary Flux Cancellation

**Companion Computational Script:** [`cell112.py`](file:///c:/data/github/connes-cvs-/cell112.py) | **Verification Log:** `cell112.out`  
**Status:** Working Research Note (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell111.md`](file:///c:/data/github/connes-cvs-/cell111.md); [`cell111.out`](file:///c:/data/github/connes-cvs-/cell111.out); [`cell110.md`](file:///c:/data/github/connes-cvs-/cell110.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary

Cell 111 established the exact row-wise resolvent identity (Theorem 1, certified to $1.09 \times 10^{-66}$ residual):
$$\boxed{\alpha_N = (H u_N)_m + \frac{T_{v_N}(0)}{\sqrt{2}} a_m - E_{11} m^2 v_{N, m} \qquad (\forall m \in \{1, \dots, N\}),}$$
and its upper-boundary mode specialization ($m = N$):
$$\boxed{\alpha_N = (H u_N)_N + \frac{T_{v_N}(0)}{\sqrt{2}} a_N - E_{11} N^2 v_{N, N}.}$$

In numerical calculations at 70 dps (`cell111.out`), this decomposition revealed that at $N = 192$:
$$\alpha_N = 1.923 \times 10^{-22}, \qquad (H u_N)_N = 9.540 \times 10^{-23}, \qquad \frac{T(0)}{\sqrt{2}} a_N = 9.687 \times 10^{-23}.$$
The boundary kinetic flux $(H u_N)_N$ is of the **exact same magnitude** as the boundary contact term $\frac{T(0)}{\sqrt{2}} a_N$, and the two terms sum precisely to $\alpha_N$.

However, the reviewer's audit of Cell 111 exposed a critical limitation:
Across $N \in \{64, 96, 128, 192\}$, the scalar $\alpha_N \approx 1.92 \times 10^{-22}$, the boundary defect $T_{v_N}(0) \approx 6.66 \times 10^{-25}$, and the second moment $S_2(N) \approx 2.48 \times 10^{-19}$ all effectively stabilized. Consequently, the extinction product $P_\alpha(N) = |\alpha_N| \sqrt{N}$ turned slightly upward (from $2.02 \times 10^{-21}$ at $N=96$ to $2.66 \times 10^{-21}$ at $N=192$), leaving the asymptotic extinction hypothesis $\alpha_N = o(N^{-1/2})$ unproven at 70 dps.

Cell 112 directly addresses this frontier through two decisive advancements:
1. **The Algebraic Origin of the $(H u_N)_N$ Coupling:** We prove analytically why $(H u_N)_N$ is algebraically locked to $\frac{T(0)}{\sqrt{2}} a_N$, originating from the exact rational identity $\frac{k^2}{N^2 - k^2} = -1 + \frac{N^2}{N^2 - k^2}$.
2. **High-Precision Extinction Audit at 110 dps:** We evaluate the ground state, boundary flux, and $\alpha_N$ across three independent mathematical routes at 110 decimal digits of precision to determine whether the plateau at $1.92 \times 10^{-22}$ is an eigensolver noise floor or a genuine asymptotic feature.

---

## 1. Algebraic Origin of Boundary Flux Cancellation

### 1.1 The Exact Decomposition of $(H u_N)_N$
Recall the definition of the boundary kinetic flux:
$$(H u_N)_N = \sum_{k=1}^N H_{Nk} k^2 v_{N, k} = H_{NN} N^2 v_{N, N} + \sum_{k=1}^{N-1} \frac{a_N - a_k}{N^2 - k^2} k^2 v_{N, k},$$
where $a_k = 2 k \psi(k) = \sqrt{2} k^2 H_{0k}$.

**Theorem 1 (Algebraic Redirection of Boundary Flux — Rigorous).**  
*The boundary kinetic flux $(H u_N)_N$ identically satisfies:*
$$(H u_N)_N = \alpha_N - \frac{T_{v_N}(0)}{\sqrt{2}} a_N + E_{11} N^2 v_{N, N},$$
*arising from the exact rational partial fraction decomposition:*
$$\frac{k^2}{N^2 - k^2} = -1 + \frac{N^2}{N^2 - k^2}.$$

*Proof.*  
Substitute $\frac{k^2}{N^2 - k^2} = -1 + \frac{N^2}{N^2 - k^2}$ into the off-diagonal sum:
$$\sum_{k=1}^{N-1} \frac{a_N - a_k}{N^2 - k^2} k^2 v_{N, k} = -\sum_{k=1}^{N-1} (a_N - a_k) v_{N, k} + N^2 \sum_{k=1}^{N-1} \frac{a_N - a_k}{N^2 - k^2} v_{N, k}.$$

We evaluate the two terms separately:
1. **The Constant Term ($-1$):**
   $$-\sum_{k=1}^{N-1} (a_N - a_k) v_{N, k} = -a_N \sum_{k=1}^{N-1} v_{N, k} + \sum_{k=1}^{N-1} a_k v_{N, k}.$$
   Using $\sum_{k=1}^{N-1} v_{N, k} = \frac{T_{v_N}(0) - v_{N, 0}}{\sqrt{2}} - v_{N, N}$ and $\sum_{k=1}^{N-1} a_k v_{N, k} = \alpha_N - a_N v_{N, N}$:
   $$-\sum_{k=1}^{N-1} (a_N - a_k) v_{N, k} = -a_N \left( \frac{T_{v_N}(0) - v_{N, 0}}{\sqrt{2}} - v_{N, N} \right) + \alpha_N - a_N v_{N, N}$$
   $$= -\frac{a_N}{\sqrt{2}} T_{v_N}(0) + \frac{a_N v_{N, 0}}{\sqrt{2}} + \alpha_N.$$

2. **The $N^2$ Term:**
   $$N^2 \sum_{k=1}^{N-1} \frac{a_N - a_k}{N^2 - k^2} v_{N, k} = N^2 \sum_{k=1}^{N-1} H_{Nk} v_{N, k}.$$
   From the $N$-th row of the ground-state eigenvalue equation $(H v_N)_N = E_{11} v_{N, N}$:
   $$H_{N0} v_{N, 0} + \sum_{k=1}^{N-1} H_{Nk} v_{N, k} + H_{NN} v_{N, N} = E_{11} v_{N, N} \implies \sum_{k=1}^{N-1} H_{Nk} v_{N, k} = E_{11} v_{N, N} - H_{N0} v_{N, 0} - H_{NN} v_{N, N}.$$
   Multiplying by $N^2$:
   $$N^2 \sum_{k=1}^{N-1} H_{Nk} v_{N, k} = E_{11} N^2 v_{N, N} - N^2 H_{N0} v_{N, 0} - H_{NN} N^2 v_{N, N}.$$
   Since $H_{N0} = \frac{\sqrt{2}}{N} \psi(N)$, we have $N^2 H_{N0} = \sqrt{2} N \psi(N) = \frac{a_N}{\sqrt{2}}$. Therefore:
   $$-N^2 H_{N0} v_{N, 0} = -\frac{a_N v_{N, 0}}{\sqrt{2}}.$$

Combining (1) and (2), the large cross-term $\frac{a_N v_{N, 0}}{\sqrt{2}}$ cancels identically:
$$\sum_{k=1}^{N-1} H_{Nk} k^2 v_{N, k} = \left( -\frac{a_N}{\sqrt{2}} T_{v_N}(0) + \frac{a_N v_{N, 0}}{\sqrt{2}} + \alpha_N \right) + \left( E_{11} N^2 v_{N, N} - \frac{a_N v_{N, 0}}{\sqrt{2}} - H_{NN} N^2 v_{N, N} \right)$$
$$= \alpha_N - \frac{a_N}{\sqrt{2}} T_{v_N}(0) + E_{11} N^2 v_{N, N} - H_{NN} N^2 v_{N, N}.$$

Adding the diagonal term $H_{NN} N^2 v_{N, N}$ gives $(H u_N)_N = \alpha_N - \frac{T_{v_N}(0)}{\sqrt{2}} a_N + E_{11} N^2 v_{N, N}$. $\quad \blacksquare$

### 1.2 The Meaning of the Numerical Equipartition
Theorem 1 proves that $(H u_N)_N$ is **not an independent dynamical variable**. Rather, the boundary flux is precisely the resolvent residual:
$$(H u_N)_N = \alpha_N - \frac{T_{v_N}(0)}{\sqrt{2}} a_N + \mathcal{O}(E_{11}).$$
When $\alpha_N \approx 2 \times \frac{T(0)}{\sqrt{2}} a_N$ (as observed in Cell 111 where $\kappa_\alpha \approx 288$ while $\frac{a_N}{\sqrt{2}} \approx \sqrt{2} N \psi(N) \approx 145$), the difference is automatically:
$$(H u_N)_N \approx \frac{T(0)}{\sqrt{2}} a_N \approx \frac{1}{2} \alpha_N.$$
The equipartition observed in Cell 111 is therefore a direct algebraic consequence of the fact that the proportionality ratio $\kappa_\alpha(N) \equiv |\alpha_N| / |T(0)|$ is approximately twice $\frac{a_N}{\sqrt{2}}$!

---

## 2. Multi-Route Evaluation of $\alpha_N$

To eliminate any ambiguity and detect potential eigensolver noise, Cell 112 evaluates $\alpha_N$ across three distinct mathematical pathways:

1. **Route 1 (Direct Modal Definition):**
   $$\alpha_N^{(1)} = \sum_{k=1}^N a_k v_{N, k} = \sum_{k=1}^N 2 k \psi(k) v_{N, k}.$$
2. **Route 2 (Boundary Row Specialization, $m = N$):**
   $$\alpha_N^{(2)} = (H u_N)_N + \frac{T_{v_N}(0)}{\sqrt{2}} a_N - E_{11} N^2 v_{N, N}.$$
3. **Route 3 (High-Sector Resolvent Average, $m > M$):**
   Averaging Theorem 1 across the upper sector $Q = \{M+1, \dots, N\}$:
   $$\alpha_N^{(3)} = \frac{1}{N - M} \sum_{m=M+1}^N \left( (H u_N)_m + \frac{T_{v_N}(0)}{\sqrt{2}} a_m - E_{11} m^2 v_{N, m} \right).$$

*Multi-Route Consistency Test:*  
If Routes 1, 2, and 3 agree to 90+ digits, the computed value of $\alpha_N$ is algebraically certified against the computed vector $v_N$.

---

## 3. Independent Eigensolver Residual Certification

To verify that the computed ground state $v_N$ is a genuine eigenvector and not distorted by ill-conditioning, Cell 112 computes:
1. **Rayleigh Quotient Error:** $|E_{11} - \langle v_N, H v_N \rangle / \|v_N\|_2^2|$.
2. **Operator Residual Norm:** $\|(H - E_{11} I) v_N\|_2$.
3. **Commutator Orthogonality Identity:**
   $$\langle v_N, \xi_N \rangle = v_0 \frac{\alpha_N}{\sqrt{2}} + \sum_{m=1}^N v_m \left( \alpha_N - \frac{T(0)}{\sqrt{2}} a_m \right) \equiv 0.$$

---

## 4. Testable Hypotheses for Cell 112

**Hypothesis 1 (Precision Floor Resolution at 110 dps).**  
The plateau $|\alpha_N| \approx 1.92 \times 10^{-22}$ observed at 70 dps is an artifact of the finite eigensolver precision floor. When evaluated at 110 dps, $|\alpha_N|$ continues to decrease for $N \ge 64$, restoring the decay of $P_\alpha(N) = |\alpha_N| \sqrt{N} \to 0$.

**Hypothesis 2 (Spectral Doublet / Ground State Branch Tracking).**  
At $N \ge 64$, the leading component $v_0$ dropped from $0.54$ to $0.064$, and $E_{11}$ transitioned from $+2.1 \times 10^{-50}$ to $-1.06 \times 10^{-51}$. We test whether this transition is an artifact of negative Archimedean finite-$T$ tail leakage or represents an excited tunneling doublet state, and track the lowest three eigenpairs $(E_0, E_1, E_2)$.

---

## 5. Pre-Flight Specification for `cell112.py`

1. **Parameters:** $c = 13$, $T = 600$, $N \in \{32, 48, 64, 80, 96, 128, 192\}$, precision `mpmath dps = 110`, `GROUND_DPS = 110`.
2. **Audit 1 (Eigensolver Residual & Spectrum Tracking):** Certify $\|(H - E_{11} I) v_N\|_2 < 10^{-100}$ and track the lowest three eigenvalues $E_0, E_1, E_2$.
3. **Audit 2 (Multi-Route $\alpha_N$ Agreement):** Verify $|\alpha_N^{(1)} - \alpha_N^{(2)}| < 10^{-95}$ and $|\alpha_N^{(1)} - \alpha_N^{(3)}| < 10^{-95}$.
4. **Audit 3 (Extinction Product Progression at 110 dps):** Measure $P_T(N) = |T(0)| N^{3/2}$, $P_\alpha(N) = |\alpha_N| \sqrt{N}$, and $\kappa_\alpha(N) = |\alpha_N| / |T(0)|$.
5. **Audit 4 (High-Sector Defect Norm $\|\xi_N^{(Q)}\|_2$):** Track exact source norm and triangle bound across $M \in \{24, 32, 48, 64\}$.
6. **Audit 5 (Scorecard & Sentinel):** Clean tabular output and explicit completion sentinel.
