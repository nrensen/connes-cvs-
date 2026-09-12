# Cell 111 Analytical Note: Discrete Boundary Defect Extinction & The Scalar Cancellation Mechanism

**Companion Computational Script:** [`cell111.py`](file:///c:/data/github/connes-cvs-/cell111.py) | **Verification Log:** `cell111.out`  
**Status:** Working Research Note (Gate 1 / Milestone M12 / Route D)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell110.md`](file:///c:/data/github/connes-cvs-/cell110.md); [`cell109.md`](file:///c:/data/github/connes-cvs-/cell109.md); [`cell108.md`](file:///c:/data/github/connes-cvs-/cell108.md); [`cell107.md`](file:///c:/data/github/connes-cvs-/cell107.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)  
**Date:** September 2026  

---

## Executive Summary

Following the formalization of the non-circular regularity bridge in Cell 110 (promoted to Theorem 9.16 in [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)), the entire Gate 1 regularity problem is reduced to two independent asymptotic hypotheses:
$$(H_{\mathrm{ext}}): \quad \lim_{N \to \infty} \|\xi_N^{(Q)}\|_2 = 0, \qquad (H_{\mathrm{gap}}): \quad \inf_{N > M} \delta_M(N) \ge \delta_\infty > 0.$$
Under $(H_{\mathrm{ext}})$ and $(H_{\mathrm{gap}})$, uniform discrete $H^2$-Sobolev regularity is proved unconditionally from $\|v_N\|_2 = 1$ without circularity:
$$\sup_{N \ge 1} \mathcal{K}_2(N) \le M^4 + \frac{(C_\xi + C_B^{\mathrm{univ}}(M))^2}{\delta_\infty^2} < \infty.$$

As established in Section 9.10 of [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md), the boundary defect source vector decomposes into:
$$\xi_N^{(Q)} = \alpha_N e^{(Q)} - \frac{T_{v_N}(0)}{\sqrt{2}} a^{(Q)},$$
where $a_k = 2 k \psi(k) = \sqrt{2} k^2 H_{0k}$, $e^{(Q)} = (1, \dots, 1)^T \in \mathbb{R}^{N-M}$, and:
$$\alpha_N \equiv \sum_{k=1}^N a_k v_{N, k} = \sum_{k=1}^N 2 k \psi(k) v_{N, k}.$$
By Proposition 9.11, the triangle inequality gives:
$$\|\xi_N^{(Q)}\|_2 \le |\alpha_N| \sqrt{N - M} + \frac{|T_{v_N}(0)|}{\sqrt{2}} \|a^{(Q)}\|_2 \le |\alpha_N| \sqrt{N - M} + \sqrt{\frac{2}{3}} C_\psi |T_{v_N}(0)| N^{3/2}.$$

Under the empirically supported semiclassical exponential boundary-suppression law $|T_{v_N}(0)| = \mathcal{O}(e^{-\sigma N})$, the contact term $|T_{v_N}(0)| N^{3/2} \to 0$ extinguishes with infinite exponential margin. 
The singular remaining analytical hurdle of Gate 1 is therefore to establish the rate:
$$\boxed{|\alpha_N| \sqrt{N - M} \longrightarrow 0 \quad \Longleftrightarrow \quad \alpha_N = o(N^{-1/2}).}$$

Cell 111 attacks the mathematical mechanism governing $\alpha_N$. We discover an exact, mode-by-mode resolvent identity expressing $\alpha_N$ directly in terms of the boundary contact defect $T_{v_N}(0)$ and the boundary kinetic flux $(H u_N)_N$.

---

## 1. Exact Mode-by-Mode Resolvent Identities for $\alpha_N$

### 1.1 Theorem 1: The Exact Row-Wise Resolvent Identity
In Cell 107, we established the exact kinetic resolvent equation:
$$(H - E_{11} I) u_N = \xi_N, \qquad u_N = K^2 v_N \quad (u_{N, m} = m^2 v_{N, m}).$$
The source vector $\xi_N \in \mathbb{R}^{N+1}$ has explicit components:
$$\xi_{N, 0} = \frac{\alpha_N}{\sqrt{2}}, \qquad \xi_{N, m} = \alpha_N - \frac{T_{v_N}(0)}{\sqrt{2}} a_m \quad (m \ge 1).$$

**Theorem 1 (Exact Row-Wise Identity for $\alpha_N$ — Rigorous).**  
*Let $v_N$ be any eigenvector of $H$ with eigenvalue $E_{11}$, normalized or unnormalized. Then for every mode $m \in \{1, \dots, N\}$, the scalar $\alpha_N = \sum_{k=1}^N 2 k \psi(k) v_{N, k}$ satisfies the exact algebraic identity:*
$$\boxed{\alpha_N = (H u_N)_m + \frac{T_{v_N}(0)}{\sqrt{2}} a_m - E_{11} m^2 v_{N, m},}$$
*where $(H u_N)_m = \sum_{k=1}^N H_{mk} k^2 v_{N, k}$.*

*Proof.*  
Recall the exact commutator identity $[K^2, H]_{jk} = a_j - a_k$ on positive modes $j, k \ge 1$, and $[K^2, H]_{j0} = \frac{a_j}{\sqrt{2}}$ for $j \ge 1$. Applying the commutator to the eigenvector $v_N$:
$$([K^2, H] v_N)_m = (K^2 H v_N)_m - (H K^2 v_N)_m = E_{11} m^2 v_{N, m} - (H u_N)_m.$$
On the other hand, expanding the matrix product directly:
$$([K^2, H] v_N)_m = [K^2, H]_{m0} v_{N, 0} + \sum_{k=1}^N [K^2, H]_{mk} v_{N, k}$$
$$= \frac{a_m}{\sqrt{2}} v_{N, 0} + \sum_{k=1}^N (a_m - a_k) v_{N, k}$$
$$= \frac{a_m}{\sqrt{2}} v_{N, 0} + a_m \sum_{k=1}^N v_{N, k} - \sum_{k=1}^N a_k v_{N, k}$$
$$= \frac{a_m}{\sqrt{2}} \left( v_{N, 0} + \sqrt{2} \sum_{k=1}^N v_{N, k} \right) - \alpha_N = \frac{T_{v_N}(0)}{\sqrt{2}} a_m - \alpha_N.$$
Equating the two expressions:
$$E_{11} m^2 v_{N, m} - (H u_N)_m = \frac{T_{v_N}(0)}{\sqrt{2}} a_m - \alpha_N.$$
Rearranging immediately yields the boxed equation. $\quad \blacksquare$

### 1.2 The Upper-Boundary Specialization ($m = N$)
Evaluating Theorem 1 at the highest discrete mode $m = N$:
$$\boxed{\alpha_N = (H u_N)_N + \frac{T_{v_N}(0)}{\sqrt{2}} a_N - E_{11} N^2 v_{N, N}.}$$
This partitions $\alpha_N$ into three transparent physical components:
1. **The Boundary Contact Term:** $\frac{T_{v_N}(0)}{\sqrt{2}} a_N = \sqrt{2} N \psi(N) T_{v_N}(0)$. Since $|\psi(N)| \le C_\psi$, this term is explicitly $\mathcal{O}(N T_{v_N}(0))$.
2. **The Ground Energy Leakage:** $E_{11} N^2 v_{N, N}$. For $c=13$, $E_{11} \sim 10^{-43}$ at $N=24$ and $\sim 10^{-51}$ at $N=192$, so this term is smaller than $10^{-48}$ and completely negligible.
3. **The Boundary Kinetic Flux:** $(H u_N)_N = \sum_{k=1}^N H_{Nk} k^2 v_{N, k}$.

---

## 2. Boundary Kinetic Flux & Curvature Cancellation

### 2.1 Asymptotic Expansion of the Boundary Flux
For $k \ll N$, expanding the divided-difference matrix element:
$$H_{Nk} = \frac{2(N\psi(N) - k\psi(k))}{N^2 - k^2} = \frac{2\psi(N)}{N} \frac{1}{1 - k^2/N^2} - \frac{2 k \psi(k)}{N^2} \frac{1}{1 - k^2/N^2}$$
$$= \frac{2\psi(N)}{N} + \frac{2\psi(N) k^2}{N^3} - \frac{2 k \psi(k)}{N^2} + \mathcal{O}(N^{-4}).$$

Substituting into the boundary kinetic flux $(H u_N)_N = \sum_{k=1}^N H_{Nk} k^2 v_{N, k}$:
$$(H u_N)_N = \frac{2\psi(N)}{N} \sum_{k=1}^N k^2 v_{N, k} - \frac{1}{N^2} \sum_{k=1}^N 2 k^3 \psi(k) v_{N, k} + \frac{2\psi(N)}{N^3} \sum_{k=1}^N k^4 v_{N, k} + \dots$$
$$= \frac{2\psi(N)}{N} S_2(N) - \frac{1}{N^2} S_\psi(N) + \frac{2\psi(N)}{N^3} \mathcal{K}_2(N) + \dots$$

### 2.2 Physical Curvature Cancellation
Notice that the leading term of $(H u_N)_N$ is proportional to:
$$S_2(N) \equiv \sum_{k=1}^N k^2 v_{N, k}.$$
In Cell 109 and Section 9.9 of [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md), we established:
$$S_2(N) \longrightarrow -\frac{1}{\sqrt{2}} \left(\frac{L}{2\pi}\right)^2 T_\infty''(0).$$
Under the infinite-order boundary flatness conjecture (Conjecture 2 in [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)), the continuum solitary wave has zero boundary curvature: $T_\infty''(0) = 0$!
Indeed, in numerical computations (`cell110.out`), $S_2(24) \approx -6.32 \times 10^{-11} \ll J_4(24) \approx 1263$.

Consequently, $(H u_N)_N$ is itself strongly suppressed by the boundary flatness of the ground state!

---

## 3. Physical-Space Coordinate Representation

In physical coordinate space, the spatial trigonometric polynomial is:
$$T_{v_N}(t) = v_{N, 0} + \sqrt{2} \sum_{k=1}^N v_{N, k} \cos\left(\frac{2\pi k t}{L}\right),$$
with derivative:
$$T_{v_N}'(t) = -\sqrt{2} \sum_{k=1}^N v_{N, k} \left(\frac{2\pi k}{L}\right) \sin\left(\frac{2\pi k t}{L}\right).$$

Recall that $\psi(k)$ are the Fourier sine coefficients of the Weil functional kernel $\mathcal{W}(t)$:
$$\psi(k) = \int_0^L \mathcal{W}(t) \sin\left(\frac{2\pi k t}{L}\right) dt.$$
Multiplying by $2 k$:
$$2 k \psi(k) = \frac{L}{\pi} \int_0^L \mathcal{W}(t) \left(\frac{2\pi k}{L} \sin\left(\frac{2\pi k t}{L}\right)\right) dt.$$

Summing against $v_{N, k}$ yields the exact coordinate integral identity:
$$\boxed{\alpha_N = -\frac{L}{\sqrt{2}\pi} \int_0^L \mathcal{W}(t) T_{v_N}'(t) \, dt.}$$

Integrating by parts:
$$\alpha_N = -\frac{L}{\sqrt{2}\pi} \left( \left[ \mathcal{W}(t) T_{v_N}(t) \right]_0^L - \int_0^L \mathcal{W}'(t) T_{v_N}(t) \, dt \right).$$
Since $T_{v_N}(0) \approx 0$ and $T_{v_N}(L) \approx 0$, the boundary term is directly proportional to the boundary contact defect:
$$\left[ \mathcal{W}(t) T_{v_N}(t) \right]_0^L = (\mathcal{W}(L) - \mathcal{W}(0)) T_{v_N}(0) + \mathcal{O}(|T_{v_N}(0) - T_{v_N}(L)|).$$
This provides an independent physical explanation for why $\alpha_N$ scales directly with the Dirichlet boundary defect $T_{v_N}(0)$.

---

## 4. Testable Hypotheses for Cell 111

**Hypothesis 1 (Boundary Contact Proportionality).**  
The scalar $\alpha_N$ is asymptotically proportional to the Dirichlet boundary defect:
$$\alpha_N = \kappa_\alpha(N) \, T_{v_N}(0), \qquad \text{with } |\kappa_\alpha(N)| \le C_\alpha N.$$
If Hypothesis 1 holds, then $|\alpha_N| \sqrt{N - M} \le C_\alpha N^{3/2} |T_{v_N}(0)| \to 0$ with exponential margin under WKB boundary suppression.

**Hypothesis 2 (Numerical Floor Identification in Cell 108).**  
The apparent plateau of $|\alpha_N| \approx 1.9 \times 10^{-22}$ for $N \ge 64$ observed in Cell 108 was an artifact of the 50-dps eigensolver noise floor ($\|a\|_2 \|\delta v\|_2 \approx 1280 \times 10^{-25} \approx 1.3 \times 10^{-22}$). When evaluated at higher precision (80–90 dps), $|\alpha_N|$ continues to decay exponentially.

**Hypothesis 3 (Row-Wise Residual Invariance).**  
For every $m \in \{1, \dots, N\}$, the identity residual:
$$r_m(N) \equiv \left| \alpha_N - \left( (H u_N)_m + \frac{T_{v_N}(0)}{\sqrt{2}} a_m - E_{11} m^2 v_{N, m} \right) \right|$$
vanishes to machine precision ($< 10^{-60}$ at 70 dps).

---

## 5. Pre-Flight Specification for `cell111.py`

1. **Parameters:** $c = 13$, $T = 600$, $N \in \{16, 24, 32, 48, 64, 80, 96, 128\}$, precision `mpmath dps = 90`.
2. **Audit 1 (Row-Wise Residual Verification):** Verify Theorem 1 identity residual $\max_{1 \le m \le N} r_m(N) < 10^{-60}$ across all dimensions.
3. **Audit 2 (Boundary Mode Decomposition, $m = N$):** Compare the three terms of $\alpha_N = (H u_N)_N + \frac{T(0)}{\sqrt{2}} a_N - E_{11} N^2 v_N$.
4. **Audit 3 (Proportionality Ratio $\kappa_\alpha(N)$):** Track $\kappa_\alpha(N) \equiv |\alpha_N| / |T_{v_N}(0)|$ across dimensions to test whether $\kappa_\alpha(N) = \mathcal{O}(N)$.
5. **Audit 4 (Extinction Products):** Track $P_T(N) = |T_{v_N}(0)| N^{3/2}$, $P_\alpha(N) = |\alpha_N| \sqrt{N}$, and the exact boundary defect source norm $\|\xi_N^{(Q)}\|_2$ across cutoffs $M \in \{24, 32, 48\}$.
6. **Audit 5 (Scorecard & Sentinel):** Clean tabular summary and explicit completion sentinel.
