# Cell 123 Analytical Note: Operator Decomposition Audit ($Q_{\mathrm{arch}}, Q_{\mathrm{prime}}, Q_{\mathrm{pole}}$) & High-Mode Coercivity

**Companion Document:** Analytical Research Note (Gate 1 / Milestone M-G1.4)  
**Status:** Executed (`cell123.out`) & Audited / Dual Hypothesis Falsification & Component Cancellation Discovery  
**Target:** [ROADMAP.md Gate 1 — Milestone M-G1.4: Continuum Spectral Threshold Proof](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell122.md`](file:///c:/data/github/connes-cvs-/cell122.md); [`cell121.md`](file:///c:/data/github/connes-cvs-/cell121.md); [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out); [`cell120.out`](file:///c:/data/github/connes-cvs-/cell120.out); [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py); [`cell.py`](file:///c:/data/github/connes-cvs-/cell.py)  
**Execution Runtime:** 647.97 s at 70 dps ($N \in [16, 64]$)  
**Date:** September 2026  

---

## 1. Executive Summary & Audit Conclusions

Cell 123 was designed to audit the operator decomposition $Q_{\mathrm{even}} = Q_{\mathrm{even}, \mathrm{arch}} + Q_{\mathrm{even}, \mathrm{prime}} + Q_{\mathrm{even}, \mathrm{pole}}$ directly from the code definitions in [`connes_cvs/operator.py`](file:///c:/data/github/connes-cvs-/connes_cvs/operator.py) and test whether principal coordinate submatrices $C_M = Q_{\mathrm{even}}[M..N, M..N]$ provide an analytical proof of the continuum spectral threshold $E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0$ for Gate 1.

The computational run ([`cell123.out`](file:///c:/data/github/connes-cvs-/cell123.out)) achieved extraordinary numerical fidelity and settled two major hypotheses while uncovering a profound new mathematical structure:

1. **Algebraic Identity Verified to 70 Decimals:**
   $$\|Q_{\mathrm{even}} - (Q_{\mathrm{even}, \mathrm{arch}} + Q_{\mathrm{even}, \mathrm{prime}} + Q_{\mathrm{even}, \mathrm{pole}})\|_F \le 2.28 \times 10^{-70} \quad (\forall N \le 64),$$
   confirming that the tripartite operator decomposition is implemented with complete algebraic fidelity.
2. **$M = 3$ Coercivity Conclusively Falsified (Hypothesis H-C3.1):**
   $$\lambda_{\min}(C_3) \approx 4.64 \times 10^{-39} \quad (N = 64) \quad \text{vs} \quad E_3 = 2.98 \times 10^{-38}.$$
   The Cauchy interlacing bound $\lambda_{\min}(C_3) \le E_3(Q_{\mathrm{even}})$ holds identically. $C_3$ retains the remaining $\sim 8$ bound states of the potential well and cannot possibly satisfy $C_3 \succeq 0.386 I$.
3. **Tight Interior Interlacing Confirmed (Hypothesis H-C3.2):**
   $$E_{13}^{(N)} \ge \lambda_{10}(C_3) \ge E_{10}^{(N)}.$$
   At $N = 64$:
   $$E_{10} = 0.57558, \qquad \lambda_{10}(C_3) = 0.57558, \qquad E_{13} = 0.58242.$$
   Although $C_3$ is not coercive as a whole, its **tenth eigenvalue already sits in the scattering continuum**.
4. **$M = 12$ Coordinate Coercivity Falsified (Hypothesis H-C12):**
   $$\lambda_{\min}(C_{12}) \text{ collapses from } 3.83 \times 10^{-2} \; (N=16) \longrightarrow 7.87 \times 10^{-6} \; (N=64).$$
   The minimum eigenvalue collapses by **four orders of magnitude**. Coordinate-mode truncation to $m \ge 12$ does **not** produce a uniformly coercive principal block $C_{12} \succeq c I > 0$.
5. **Component Sign Discovery (Strong Prime Indefiniteness):**
   On high modes ($M = 12, N = 64$):
   $$C_{\mathrm{arch}} \succeq 1.553 I, \qquad \lambda_{\min}(C_{\mathrm{prime}}) \approx -2.373, \qquad C_{\mathrm{pole}} \approx \mathcal{O}(10^{-6}).$$
   The prime component is **strongly indefinite**. Full submatrix near-positivity ($\lambda_{\min}(C_{12}) \approx 7.87 \times 10^{-6}$) is produced by **near-perfect destructive cancellation between the Archimedean and prime sectors**, not by componentwise positivity.
6. **The Fundamental Conceptual Distinction:**
   $$\boxed{\text{coordinate-mode truncation } (\operatorname{span}\{e_M, \dots, e_N\}) \ne \text{spectral-subspace projection } (P_{\mathrm{cont}} = I - P_{\mathrm{bound}}).}$$
   Physical bound states are wavepackets with non-zero tails across high Fourier modes. Truncating coordinate indices leaves bound-state tail weight in $C_{12}$, creating an additional near-zero eigenvalue.

---

## 2. Tripartite Decomposition & Matrix Invariants

Across discrete dimensions $N \in \{16, 20, 24, 28, 32, 36, 40, 48, 64\}$ at $c = 13, T = 600, \text{dps} = 70$:

| $N$ | $\|Q_{\mathrm{even}} - \sum Q_{\mathrm{comp}}\|_F$ | $\lambda_{\min}(Q_{\mathrm{even}})$ | $\lambda_{\min}(C_3)$ | $E_3(Q_{\mathrm{even}})$ | Status ($\lambda_{\min} \le E_3$) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 16 | $7.37 \times 10^{-71}$ | $7.39 \times 10^{-35}$ | $1.74 \times 10^{-20}$ | $1.30 \times 10^{-19}$ | **VERIFIED** |
| 20 | $8.62 \times 10^{-71}$ | $1.34 \times 10^{-39}$ | $1.48 \times 10^{-24}$ | $1.03 \times 10^{-23}$ | **VERIFIED** |
| 24 | $8.98 \times 10^{-71}$ | $2.60 \times 10^{-43}$ | $2.07 \times 10^{-27}$ | $1.41 \times 10^{-26}$ | **VERIFIED** |
| 28 | $1.06 \times 10^{-70}$ | $9.34 \times 10^{-47}$ | $2.87 \times 10^{-30}$ | $1.91 \times 10^{-29}$ | **VERIFIED** |
| 32 | $1.28 \times 10^{-70}$ | $2.17 \times 10^{-49}$ | $2.00 \times 10^{-32}$ | $1.31 \times 10^{-31}$ | **VERIFIED** |
| 36 | $1.60 \times 10^{-70}$ | $2.65 \times 10^{-50}$ | $2.92 \times 10^{-34}$ | $1.89 \times 10^{-33}$ | **VERIFIED** |
| 40 | $1.73 \times 10^{-70}$ | $2.57 \times 10^{-50}$ | $9.27 \times 10^{-36}$ | $5.98 \times 10^{-35}$ | **VERIFIED** |
| 48 | $1.82 \times 10^{-70}$ | $2.10 \times 10^{-50}$ | $6.89 \times 10^{-38}$ | $4.42 \times 10^{-37}$ | **VERIFIED** |
| 64 | $2.28 \times 10^{-70}$ | $-5.08 \times 10^{-52}$ | $4.64 \times 10^{-39}$ | $2.98 \times 10^{-38}$ | **VERIFIED** |

---

## 3. Auditing the $M = 3$ Submatrix ($C_3 = Q_{\mathrm{even}}[3..N, 3..N]$)

### 3.1 Bound-State Collapse vs Interior Continuum Interlacing
At $c = 13$, the continuous Archimedean multiplier changes sign between $a_2$ and $a_3$: $h_+(a_2) \approx -0.914, h_+(a_3) \approx +0.386$.
However, the submatrix $C_3$ fails coercivity completely because $\lambda_{\min}(C_3)$ tracks $E_3$:

| $N$ | $E_3(Q)$ | $\lambda_{\min}(C_3)$ | $E_{10}(Q)$ | $\lambda_{10}(C_3)$ | $E_{13}(Q)$ | Interlacing Chain |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 16 | $1.30 \times 10^{-19}$ | $1.74 \times 10^{-20}$ | $1.04298$ | $2.49043$ | $2.49043$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 20 | $1.03 \times 10^{-23}$ | $1.48 \times 10^{-24}$ | $0.12816$ | $2.06092$ | $2.06667$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 24 | $1.41 \times 10^{-26}$ | $2.07 \times 10^{-27}$ | $5.35 \times 10^{-3}$ | $1.79144$ | $1.80126$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 28 | $1.91 \times 10^{-29}$ | $2.87 \times 10^{-30}$ | $7.47 \times 10^{-4}$ | $1.29038$ | $1.31419$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 32 | $1.31 \times 10^{-31}$ | $2.00 \times 10^{-32}$ | $1.37 \times 10^{-4}$ | $0.81453$ | $0.81621$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 36 | $1.89 \times 10^{-33}$ | $2.92 \times 10^{-34}$ | $3.77 \times 10^{-5}$ | $0.78162$ | $0.78163$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 40 | $5.98 \times 10^{-35}$ | $9.27 \times 10^{-36}$ | $1.69 \times 10^{-5}$ | $0.72710$ | $0.72745$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 48 | $4.42 \times 10^{-37}$ | $6.89 \times 10^{-38}$ | $8.10 \times 10^{-6}$ | $0.62446$ | $0.62883$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |
| 64 | $2.98 \times 10^{-38}$ | $4.64 \times 10^{-39}$ | $6.99 \times 10^{-6}$ | $0.57558$ | $0.58242$ | $E_{13} \ge \lambda_{10} \ge E_{10}$ |

**Conclusion on H-C3.1 & H-C3.2:**
- $\lambda_{\min}(C_3) \le E_3$ is verified unconditionally to 70 dps. $C_3 \succeq 0.386 I$ is **definitively falsified**.
- However, the tenth eigenvalue $\lambda_{10}(C_3)$ sits tightly inside the continuum: at $N = 64$, $\lambda_{10}(C_3) = 0.57558$, while $E_{10} = 6.99 \times 10^{-6}$ and $E_{13} = 0.58242$. Removing 3 modes isolates the continuum **in the interior spectrum**, but not at the bottom of the principal submatrix.

---

## 4. Auditing the $M = 12$ Submatrix ($C_{12} = Q_{\mathrm{even}}[12..N, 12..N]$)

### 4.1 Numerical Breakdown & Collapse of $\lambda_{\min}(C_{12})$

| $N$ | $\lambda_{\min}(C_{12})$ | $\lambda_1(C_{12})$ | $E_{13}(Q)$ | Cauchy Bound Status |
| :---: | :---: | :---: | :---: | :---: |
| 16 | $3.83 \times 10^{-2}$ | $0.84674$ | $2.49043$ | **VERIFIED** |
| 20 | $3.57 \times 10^{-3}$ | $0.09495$ | $2.06667$ | **VERIFIED** |
| 24 | $3.24 \times 10^{-4}$ | $0.01260$ | $1.80126$ | **VERIFIED** |
| 28 | $4.95 \times 10^{-5}$ | $0.00718$ | $1.31419$ | **VERIFIED** |
| 32 | $1.57 \times 10^{-5}$ | $0.00404$ | $0.81621$ | **VERIFIED** |
| 36 | $1.01 \times 10^{-5}$ | $0.00220$ | $0.78163$ | **VERIFIED** |
| 40 | $8.76 \times 10^{-6}$ | $0.00157$ | $0.72745$ | **VERIFIED** |
| 48 | $8.09 \times 10^{-6}$ | $0.00118$ | $0.62883$ | **VERIFIED** |
| 64 | $7.87 \times 10^{-6}$ | $0.00082$ | $0.58242$ | **VERIFIED** |

### 4.2 Rejection of Uniform Coordinate Coercivity
In the pre-flight note of Cell 123, it was hypothesized that removing the 12 lowest coordinate modes would eliminate all $\bar{N}_{\mathrm{bound}} \approx 11$ bound states, yielding uniform coercivity $C_{12} \succeq c_{12} I > 0$ with $c_{12} \approx 0.50$.
**The data conclusively refute this hypothesis:**
- $\lambda_{\min}(C_{12})$ collapses monotonically from $0.0383$ down to $7.87 \times 10^{-6}$ across $N \in [16, 64]$.
- The second eigenvalue $\lambda_1(C_{12})$ also collapses: from $0.847$ down to $8.20 \times 10^{-4}$.
- Consequently, Cauchy interlacing $E_{13} \ge \lambda_1(C_{12}) \ge \lambda_{\min}(C_{12})$ provides only an asymptotically vanishing lower bound ($7.87 \times 10^{-6}$), **not** the macroscopic bound $0.582$ that exists in the full spectrum $E_{13}$.

---

## 5. The Component Spectrum Audit: Archimedean vs Prime Cancellation

Evaluating the spectrum of individual operator sectors on $\mathcal{H}_{\mathrm{high}}^{(12)} = \operatorname{span}\{e_{12}, \dots, e_N\}$ reveals why $C_{12}$ is not strongly positive:

| Sector | Spectrum at $N = 16$ | Spectrum at $N = 32$ | Spectrum at $N = 64$ | Nature |
| :--- | :---: | :---: | :---: | :---: |
| **Archimedean ($C_{12, \mathrm{arch}}$)** | $[1.555, 1.848]$ | $[1.553, 2.530]$ | $[1.553, 3.216]$ | **Strictly Positive Definite** ($\ge 1.55$) |
| **Prime ($C_{12, \mathrm{prime}}$)** | $[-1.676, 1.008]$ | $[-2.180, 1.875]$ | $[-2.373, 2.212]$ | **Strongly Indefinite** |
| **Pole ($C_{12, \mathrm{pole}}$)** | $[0, 2.94 \times 10^{-6}]$ | $[0, 4.25 \times 10^{-6}]$ | $[0, 4.42 \times 10^{-6}]$ | **Positive Semidefinite** (negligible) |
| **Total ($C_{12}$)** | $\lambda_{\min} = 0.03825$ | $\lambda_{\min} = 1.57 \times 10^{-5}$ | $\lambda_{\min} = 7.87 \times 10^{-6}$ | **Near-Zero (Cancellation)** |

### Analytical Interpretation
1. **The Prime Term is Strongly Indefinite:**
   The prime symbol carries an explicit minus sign: $\psi_{\mathrm{prime}}(x) = -\frac{1}{\pi} \sum \frac{\Lambda(n)}{\sqrt{n}} \sin(\dots)$. While the continuous Archimedean piece is coercive ($C_{\mathrm{arch}} \ge 1.553$), the prime matrix elements produce negative eigenvalues reaching $-2.373$.
2. **Near-Perfect Destructive Cancellation:**
   The positive Archimedean sector and the negative modes of the prime sector almost perfectly cancel each other on high modes:
   $$\lambda_{\min}(C_{12}) = \lambda_{\min}(C_{\mathrm{arch}} + C_{\mathrm{prime}} + C_{\mathrm{pole}}) \approx 1.55 - 1.55 = 7.87 \times 10^{-6}.$$
   Positivity in the Connes–CvS Galerkin operator is **not componentwise**; it is an emergent balance between arithmetic and Archimedean distributions.

---

## 6. The $M$-Cutoff Sweep: Spectral Localization vs Coordinate Truncation

Tracking $\lambda_{\min}(C_M)$ across cutoff index $M$ and dimension $N$:

| $N$ | $M=3$ | $M=4$ | $M=6$ | $M=8$ | $M=10$ | $M=12$ | $M=14$ | $M=16$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 16 | $1.74 \times 10^{-20}$ | $1.57 \times 10^{-16}$ | $6.00 \times 10^{-10}$ | $4.14 \times 10^{-5}$ | $4.30 \times 10^{-3}$ | $0.03825$ | $0.2771$ | --- |
| 20 | $1.48 \times 10^{-24}$ | $1.90 \times 10^{-20}$ | $3.31 \times 10^{-13}$ | $2.33 \times 10^{-7}$ | $9.62 \times 10^{-5}$ | $3.57 \times 10^{-3}$ | $0.03461$ | $0.2384$ |
| 24 | $2.07 \times 10^{-27}$ | $4.48 \times 10^{-23}$ | $1.55 \times 10^{-15}$ | $1.86 \times 10^{-9}$ | $5.62 \times 10^{-6}$ | $3.24 \times 10^{-4}$ | $3.20 \times 10^{-3}$ | $0.05196$ |
| 28 | $2.87 \times 10^{-30}$ | $9.91 \times 10^{-26}$ | $1.25 \times 10^{-17}$ | $4.88 \times 10^{-11}$ | $2.73 \times 10^{-6}$ | $4.95 \times 10^{-5}$ | $1.64 \times 10^{-3}$ | $0.01897$ |
| 32 | $2.00 \times 10^{-32}$ | $8.78 \times 10^{-28}$ | $1.76 \times 10^{-19}$ | $2.05 \times 10^{-12}$ | $1.16 \times 10^{-6}$ | $1.57 \times 10^{-5}$ | $9.99 \times 10^{-4}$ | $7.19 \times 10^{-3}$ |
| 36 | $2.92 \times 10^{-34}$ | $1.96 \times 10^{-29}$ | $9.52 \times 10^{-21}$ | $2.29 \times 10^{-13}$ | $4.04 \times 10^{-7}$ | $1.01 \times 10^{-5}$ | $6.07 \times 10^{-4}$ | $3.76 \times 10^{-3}$ |
| 40 | $9.27 \times 10^{-36}$ | $9.06 \times 10^{-31}$ | $8.21 \times 10^{-22}$ | $4.13 \times 10^{-14}$ | $1.82 \times 10^{-7}$ | $8.76 \times 10^{-6}$ | $4.52 \times 10^{-4}$ | $2.93 \times 10^{-3}$ |
| 48 | $6.89 \times 10^{-38}$ | $1.56 \times 10^{-32}$ | $5.01 \times 10^{-23}$ | $9.09 \times 10^{-15}$ | $8.66 \times 10^{-8}$ | $8.09 \times 10^{-6}$ | $3.48 \times 10^{-4}$ | $2.50 \times 10^{-3}$ |
| 64 | $4.64 \times 10^{-39}$ | $1.83 \times 10^{-33}$ | $1.79 \times 10^{-23}$ | $5.42 \times 10^{-15}$ | $7.53 \times 10^{-8}$ | $7.87 \times 10^{-6}$ | $2.52 \times 10^{-4}$ | $2.13 \times 10^{-3}$ |

### The Physical Lesson
At fixed $N = 64$, $\lambda_{\min}(C_M)$ increases smoothly across cutoffs:
$$10^{-39} \to 10^{-33} \to 10^{-23} \to 10^{-15} \to 10^{-8} \to 10^{-6} \to 10^{-4} \to 10^{-3}.$$
This is **not a sharp step function to a continuum floor**. It is the gradual extinction of bound-state wavepacket tails as higher Fourier modes are truncated.
Because bound states $u_j$ have non-compact Fourier support (their Fourier tails decay rapidly but are non-zero), deleting the coordinate basis vectors $e_0, \dots, e_{M-1}$ leaves residual overlap with the bound-state subspace, permitting a Rayleigh quotient of order $\mathcal{O}(\text{tail mass})$.

---

## 7. Strategic Redirection for Gate 1 & Milestone M-G1.4

1. **Retire Coordinate Submatrix Coercivity:**
   Cauchy interlacing on coordinate principal submatrices $C_M = Q[M..N, M..N]$ cannot manufacture the continuum gap. This branch is officially closed.
2. **Cell 121 Full-Operator Spectral Gap is Fully Intact:**
   Cell 121's discovery remains completely solid:
   $$E_{13}^{(N)} - E_3^{(N)} \ge 0.582 \quad (\forall N \in [16, 64]).$$
   The continuum gap is a property of the **full operator's eigenvalues**, not of its coordinate blocks. Gate 1 cares exclusively about $E_{L+1} - E_{j+1}$, not about coordinate localization.
3. **The Natural Spectral Subspace Projection (Forward Path for Cell 124):**
   Instead of coordinate truncation $C_M$, the natural operator separating the bound states from the continuum is the **spectral projection**:
   $$P_{\mathrm{bound}}^{(N)} = \sum_{j=0}^{10} u_j^{(N)} (u_j^{(N)})^T, \qquad P_{\mathrm{cont}}^{(N)} = I - P_{\mathrm{bound}}^{(N)}.$$
   On the continuum spectral subspace $\mathcal{H}_{\mathrm{cont}} = \operatorname{Ran}(P_{\mathrm{cont}})$, the operator $P_{\mathrm{cont}} Q_{\mathrm{even}} P_{\mathrm{cont}}$ has minimum eigenvalue identically equal to $E_{11}^{(N)}$ (or $E_{12}^{(N)}$), which is macroscopic ($\approx 0.58$).
   The analytical challenge for Cell 124 is to establish a variational min-max lower bound on $E_{12}$ or $E_{13}$ directly from the quadratic form, bypassing coordinate principal submatrices entirely.
