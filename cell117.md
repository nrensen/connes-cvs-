# Cell 117 Analytical Note: Boundary-Row Modal Profiling & Boundary Layer Scaling Collapse Diagnostic

**Companion Computational Script:** [`cell117.py`](file:///c:/data/github/connes-cvs-/cell117.py) | **Verification Log:** [`cell117.out`](file:///c:/data/github/connes-cvs-/cell117.out) (runtime: 1412.63 s at 70 dps)  
**Status:** Executed & Audited / Decisive Negative Result (Airy Hypothesis Refuted & Retired) (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell116.md`](file:///c:/data/github/connes-cvs-/cell116.md); [`cell116.out`](file:///c:/data/github/connes-cvs-/cell116.out); [`cell115.md`](file:///c:/data/github/connes-cvs-/cell115.md); [`cell114.md`](file:///c:/data/github/connes-cvs-/cell114.md); [`cell112a.md`](file:///c:/data/github/connes-cvs-/cell112a.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary & Calibrated Verdict

Cell 117 was authored to directly test whether the empirical $N^{-1/3}$ scaling observed in Cells 115–116 arises from a true physical boundary-layer mechanism (the Airy Dirichlet boundary layer hypothesis) or was merely an empirical coincidence.

### Headline Verdict

> **Classification: Decisive Negative Diagnostic. The Airy boundary-layer hypothesis is refuted and formally retired.**
>
> 1. **Complete Absence of Scaling Collapse:**
>    - **Discrete Fixed-Depth Test ($\theta = 0$):** The rescaled boundary modes $S_{N, j} = N^{1/3} F_{N, N-j}$ fail to collapse, exhibiting massive relative spreads between $465\%$ and $1073\%$ ($j = 0$: $692.5\%$, $j = 1$: $1073.0\%$, $j = 2$: $853.2\%$).
>    - **Airy Coordinate Test ($\theta = 1/3$):** Across matched dimensionless coordinates $\eta = j / N^{1/3}$, the ratio $S_{192}(\eta) / S_{64}(\eta)$ swings wildly from $-1.48$ to $+0.008$, including multiple sign inversions. There is zero evidence of a universal boundary profile $f_{\mathrm{Airy}}(\eta)$.
>
> 2. **Analytical & Numerical Audit of $H_{0, N}$ and $a_N$ (Module 1):**
>    - Confirmed the exact operator identities $H_{0, N} = \sqrt{2}\psi(N)/N$ and $a_N = 2 N \psi(N)$ to machine precision ($10^{-70}$).
>    - Demonstrated that the prime contribution $\psi_{\mathrm{prime}}(N)$ is an almost-periodic sum of high-frequency sines with incommensurate frequencies that does not decay, oscillating indefinitely with an $\mathcal{O}(1)$ envelope.
>    - Consequently, $a_N / \sqrt{2} = \sqrt{2} N \psi(N)$ has an $\mathcal{O}(N)$ oscillatory envelope swinging from $+24.17$ to $-145.38$. This definitively disproves the Cell 116 hypothesis that the boundary behavior is governed by smooth asymptotic scaling of $a_N$.
>
> 3. **Discovery of Massive 20-Order Destructive Cancellation (Module 5):**
>    - The net boundary flux $(H u_N)_N \sim 10^{-22}$ is **NOT** governed by a localized boundary layer.
>    - The cumulative partial sums $S_{\mathrm{core}}(M) = \sum_{k=1}^M F_{N, k}$ reach a peak magnitude $M_{\mathrm{peak}} \sim 10^{-3} - 10^{-2}$ that **always occurs at $k_{\mathrm{peak}} = 3$** for every tested dimension $N \in \{64, 96, 128, 160, 192\}$.
>    - The cancellation factor $\mathcal{C}_{\mathrm{cancel}}(N) \equiv M_{\mathrm{peak}}(N) / |(H u_N)_N|$ reaches an astonishing:
>      $$\mathcal{C}_{\mathrm{cancel}}(192) = 1.3685 \times 10^{20}.$$
>    - The boundary flux is the residual of approximately **20 orders of magnitude of destructive cancellation** across the modal spectrum.
>
> 4. **Epistemic Re-Interpretation of the Cell 115 Exponent ($\beta \approx -0.29$):**
>    The sub-critical decay observed in Cell 115 is **not** the scaling exponent of a coherent local boundary layer. It is the residual of an oscillatory arithmetic cancellation. The Airy hypothesis is retired, and the investigation pivots to the global arithmetic and spectral cancellation structure.

---

## 1. Audited Numerical Results (`cell117.out`)

### 1.1 Module 1: Exact Analytical & Numerical Audit of $H_{0, N}$ and $a_N$

The symbol decomposition $\psi(N) = \psi_{\mathrm{prime}}(N) + \psi_{\mathrm{pole}}(N) + \psi_{\mathrm{arch}}(N)$ was evaluated across $N \in [32, 192]$:

| $N$ | $H_{0, N}$ | $\psi_{\mathrm{prime}}(N)$ | $\psi_{\mathrm{pole}}(N)$ | $\psi_{\mathrm{arch}}(N)$ | $\psi_{\mathrm{total}}(N)$ | $a_N / \sqrt{2}$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **32** | $+0.0100745$ | $+0.0443757$ | $+0.0076456$ | $+0.1759392$ | $+0.2279604$ | $+10.3163$ |
| **48** | $+0.0044931$ | $+0.0091107$ | $+0.0050972$ | $+0.1382934$ | $+0.1525012$ | $+10.3521$ |
| **64** | $+0.0012861$ | $-0.0453486$ | $+0.0038229$ | $+0.0997298$ | $+0.0582040$ | $+5.2680$ |
| **80** | $-0.0007180$ | $-0.1035612$ | $+0.0030583$ | $+0.0598843$ | $-0.0406185$ | $-4.5955$ |
| **96** | $-0.0018937$ | $-0.1494010$ | $+0.0025486$ | $+0.0183051$ | $-0.1285473$ | $-17.4522$ |
| **112** | $-0.0010399$ | $-0.0589561$ | $+0.0021845$ | $-0.0255844$ | $-0.0823560$ | $-13.0445$ |
| **128** | $-0.0027270$ | $-0.1761793$ | $+0.0019115$ | $-0.0725530$ | $-0.2468208$ | $-44.6793$ |
| **144** | $-0.0020769$ | $-0.0894985$ | $+0.0016991$ | $-0.1236788$ | $-0.2114782$ | $-43.0669$ |
| **160** | $+0.0009443$ | $+0.2858707$ | $+0.0015292$ | $-0.1805665$ | $+0.1068333$ | $+24.1736$ |
| **176** | $-0.0022166$ | $-0.0314502$ | $+0.0013902$ | $-0.2457929$ | $-0.2758529$ | $-68.6602$ |
| **192** | $-0.0039437$ | $-0.2127402$ | $+0.0012743$ | $-0.3239441$ | $-0.5354100$ | $-145.3794$ |

- **Analytical Finding:** $\psi_{\mathrm{prime}}(N)$ fluctuates between $-0.21$ and $+0.29$ due to incommensurate prime-power frequencies. Multiplying by $\sqrt{2} N$ produces the linear envelope growth $|a_N / \sqrt{2}| \sim \mathcal{O}(N)$ with irregular sign changes.

---

### 1.2 Module 2: Localized Branch Eigenvector Extraction & Identity Verification

The localized branch (state $k = 1$) was solved across $N \in \{64, 96, 128, 160, 192\}$:

| $N$ | Ground Weight $v_0$ | Normalization Residual | Boundary Flux $(H u_N)_N$ | Boundary Coupling $\alpha_N$ | Contact Term $T_v(0)$ | Theorem 1 Residual |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | $0.665695$ | $9.06 \times 10^{-72}$ | $-2.4678 \times 10^{-22}$ | $-2.4219 \times 10^{-22}$ | $8.7212 \times 10^{-25}$ | $9.55 \times 10^{-71}$ |
| **96** | $0.666227$ | $3.62 \times 10^{-71}$ | $-1.7041 \times 10^{-22}$ | $-1.8163 \times 10^{-22}$ | $6.4269 \times 10^{-25}$ | $7.75 \times 10^{-70}$ |
| **128** | $0.666319$ | $1.81 \times 10^{-71}$ | $-1.4577 \times 10^{-22}$ | $-1.7287 \times 10^{-22}$ | $6.0647 \times 10^{-25}$ | $8.77 \times 10^{-70}$ |
| **160** | $0.666380$ | $1.81 \times 10^{-71}$ | $-1.8553 \times 10^{-22}$ | $-1.7116 \times 10^{-22}$ | $5.9439 \times 10^{-25}$ | $9.04 \times 10^{-71}$ |
| **192** | $0.666394$ | $9.06 \times 10^{-72}$ | $-8.3491 \times 10^{-23}$ | $-1.6798 \times 10^{-22}$ | $5.8117 \times 10^{-25}$ | $2.33 \times 10^{-69}$ |

- **Verification:** Theorem 1 identity $(H u_N)_N = \alpha_N - \frac{a_N}{\sqrt{2}} T_v(0) + E_{11} N^2 v_{N, N}$ is verified unconditionally to residuals below $2.4 \times 10^{-69}$ across all dimensions.

---

### 1.3 Module 3 & 4: Refutation of the Scaling Collapse

#### Fixed-Depth Test ($\theta = 0$): Rescaled Terms $S_{N, j} = N^{1/3} F_{N, N-j}$

| $j = N - k$ | $N = 64$ | $N = 96$ | $N = 128$ | $N = 160$ | $N = 192$ | Relative Spread |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **0** | $-1.63 \times 10^{-22}$ | $-2.85 \times 10^{-23}$ | $-5.68 \times 10^{-23}$ | $-1.95 \times 10^{-22}$ | $+1.76 \times 10^{-22}$ | **$692.5 \%$** |
| **1** | $-3.48 \times 10^{-23}$ | $-2.27 \times 10^{-24}$ | $-1.73 \times 10^{-25}$ | $-6.73 \times 10^{-23}$ | $+4.99 \times 10^{-23}$ | **$1073.0 \%$** |
| **2** | $+1.42 \times 10^{-24}$ | $+8.02 \times 10^{-24}$ | $+3.90 \times 10^{-26}$ | $-2.19 \times 10^{-23}$ | $-5.12 \times 10^{-24}$ | **$853.2 \%$** |
| **3** | $-4.50 \times 10^{-23}$ | $+2.02 \times 10^{-24}$ | $-2.86 \times 10^{-25}$ | $-2.38 \times 10^{-24}$ | $-2.10 \times 10^{-24}$ | **$492.5 \%$** |
| **4** | $+3.96 \times 10^{-25}$ | $-6.03 \times 10^{-24}$ | $-1.29 \times 10^{-24}$ | $-9.80 \times 10^{-24}$ | $+5.94 \times 10^{-24}$ | **$729.6 \%$** |
| **5** | $-1.49 \times 10^{-24}$ | $+7.73 \times 10^{-24}$ | $-4.31 \times 10^{-24}$ | $-8.56 \times 10^{-24}$ | $-1.67 \times 10^{-24}$ | **$980.5 \%$** |
| **8** | $+8.36 \times 10^{-23}$ | $+7.91 \times 10^{-25}$ | $-5.33 \times 10^{-24}$ | $+1.62 \times 10^{-23}$ | $+3.14 \times 10^{-25}$ | **$465.3 \%$** |
| **20** | $+2.92 \times 10^{-20}$ | $+2.47 \times 10^{-25}$ | $-1.08 \times 10^{-24}$ | $-4.60 \times 10^{-24}$ | $+1.69 \times 10^{-27}$ | **$500.2 \%$** |

#### Airy Coordinate Test ($\theta = 1/3$): Matched $\eta = j / N^{1/3}$

| Target $\eta$ | $(N=64, j)$ | $S(N=64)$ | $(N=192, j)$ | $S(N=192)$ | Ratio $S_{192} / S_{64}$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **$0.0$** | $j = 0$ | $-1.63 \times 10^{-22}$ | $j = 0$ | $+1.76 \times 10^{-22}$ | **$-1.077$** (Sign inversion) |
| **$0.5$** | $j = 2$ | $+1.42 \times 10^{-24}$ | $j = 3$ | $-2.10 \times 10^{-24}$ | **$-1.485$** (Sign inversion) |
| **$1.0$** | $j = 4$ | $+3.96 \times 10^{-25}$ | $j = 6$ | $-1.21 \times 10^{-25}$ | **$-0.306$** (Sign inversion) |
| **$1.5$** | $j = 6$ | $-3.95 \times 10^{-24}$ | $j = 9$ | $-1.25 \times 10^{-25}$ | **$+0.032$** |
| **$2.0$** | $j = 8$ | $+8.36 \times 10^{-23}$ | $j = 12$ | $+6.74 \times 10^{-25}$ | **$+0.008$** |
| **$2.5$** | $j = 10$ | $+1.31 \times 10^{-23}$ | $j = 14$ | $+2.54 \times 10^{-24}$ | **$+0.194$** |
| **$3.0$** | $j = 12$ | $-1.83 \times 10^{-23}$ | $j = 17$ | $+4.81 \times 10^{-25}$ | **$-0.026$** (Sign inversion) |

**Definitive Conclusion:**  
Under the Airy rescaling, ratios between $N = 64$ and $N = 192$ vary from $-1.48$ to $+0.008$, accompanied by frequent sign flips. There is **no scaling collapse whatsoever**. The individual boundary terms do not conform to an Airy profile.

---

### 1.4 Module 5: Binned Decomposition & The 20-Order Cancellation

#### Sector Contributions Across Dimensions

| Sector / Bin | $N = 64$ | $N = 96$ | $N = 128$ | $N = 160$ | $N = 192$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Bin 0: Endpoint ($j=0$)** | $-4.08 \times 10^{-23}$ | $-6.22 \times 10^{-24}$ | $-1.13 \times 10^{-23}$ | $-3.60 \times 10^{-23}$ | $+3.05 \times 10^{-23}$ |
| **Bin 1: Bdy Layer ($1 \le j \le 4$)** | $-1.95 \times 10^{-23}$ | $+3.82 \times 10^{-25}$ | $-3.40 \times 10^{-25}$ | $-1.87 \times 10^{-23}$ | $+8.43 \times 10^{-24}$ |
| **Bin 2: Near-Bdy ($5 \le j \le 10$)** | $+2.30 \times 10^{-23}$ | $+6.30 \times 10^{-24}$ | $-3.02 \times 10^{-24}$ | $-2.20 \times 10^{-24}$ | $+3.81 \times 10^{-25}$ |
| **Bin 3: Intermed ($11 \le j \le 20$)** | $+1.91 \times 10^{-20}$ | $-1.21 \times 10^{-24}$ | $-3.52 \times 10^{-24}$ | $-2.74 \times 10^{-24}$ | $+9.33 \times 10^{-25}$ |
| **Bin 4: Transition ($21 \le j \le 40$)** | **$+4.98734 \times 10^{-12}$** | $+8.77 \times 10^{-25}$ | $-6.74 \times 10^{-25}$ | $-1.36 \times 10^{-24}$ | $+4.54 \times 10^{-24}$ |
| **Bin 5: Bulk Core ($j > 40$)** | **$-4.98734 \times 10^{-12}$** | $-1.7054 \times 10^{-22}$ | $-1.2693 \times 10^{-22}$ | $-1.2457 \times 10^{-22}$ | **$-1.2827 \times 10^{-22}$** |
| **Total Flux $(H u_N)_N$** | **$-2.4678 \times 10^{-22}$** | **$-1.7041 \times 10^{-22}$** | **$-1.4577 \times 10^{-22}$** | **$-1.8553 \times 10^{-22}$** | **$-8.3491 \times 10^{-23}$** |

#### Cancellation Factor Tracking Across Dimensions

| Dimension $N$ | Net Flux $(H u_N)_N$ | Peak Partial Sum $M_{\mathrm{peak}}(N)$ | Mode of Peak $k_{\mathrm{peak}}$ | Cancellation Ratio $\mathcal{C}_{\mathrm{cancel}}(N)$ |
| :---: | :---: | :---: | :---: | :---: |
| **64** | $-2.4678 \times 10^{-22}$ | $3.2637 \times 10^{-3}$ | **$k = 3$** | **$1.3225 \times 10^{19}$** |
| **96** | $-1.7041 \times 10^{-22}$ | $5.7156 \times 10^{-3}$ | **$k = 3$** | **$3.3540 \times 10^{19}$** |
| **128** | $-1.4577 \times 10^{-22}$ | $8.0075 \times 10^{-3}$ | **$k = 3$** | **$5.4933 \times 10^{19}$** |
| **160** | $-1.8553 \times 10^{-22}$ | $2.6419 \times 10^{-3}$ | **$k = 3$** | **$1.4240 \times 10^{19}$** |
| **192** | $-8.3491 \times 10^{-23}$ | $1.1426 \times 10^{-2}$ | **$k = 3$** | **$1.3685 \times 10^{20}$** |

---

## 2. Theoretical Deconstruction: Why the Airy Analogy Failed

1. **Non-Local Kernel vs Local Schrödinger Barrier:**  
   The Airy boundary layer model is derived from a local differential equation $-\frac{1}{N^2}\phi'' + F t \phi = E\phi$ with a hard Dirichlet wall at $t = 0$. The Connes–CvS Galerkin operator $H_{jk} = \frac{2(j\psi(j) - k\psi(k))}{j^2 - k^2}$ is an intrinsically non-local divided-difference operator whose entries couple every mode to the entire spectrum.
2. **Arithmetic Incommensurability vs Continuous Phase:**  
   The symbol $\psi(N)$ contains discrete high-frequency oscillations from prime powers $\sin(2\pi N \frac{\log p^k}{L})$. These trigonometric terms fluctuate rapidly without phase alignment, disrupting any smooth semiclassical boundary envelope.
3. **The Locus of Cancellation ($k_{\mathrm{peak}} = 3$):**  
   In a true boundary-layer phenomenon, the dominant contribution lives in the high-frequency modes $k \approx N$, with the bulk contributing negligible background. In `cell117.out`, the exact opposite occurs: the partial sum peaks at the low-frequency core mode $k = 3$ with magnitude $\sim 10^{-2}$, and the remaining 189 modes participate in a vast destructive cancellation that quenches the total flux down to $10^{-22}$.

---

## 3. The New Pivot: Arithmetic & Continuum Phase Cancellation

The failure of the Airy mechanism leads to a conceptual pivot:

$$\boxed{\text{Old Question: What boundary layer produces the } N^{-0.3} \text{ defect?}}$$
$$\boxed{\text{New Question: What arithmetic/spectral cancellation produces the } 10^{-22} \text{ residual?}}$$

The boundary defect $\alpha_N = \sum_k a_k v_k$ is the residual of an oscillatory sum whose individual summands are 20 orders of magnitude larger. Fitting a power law to such a cancellation residual yields an effective pre-asymptotic exponent that reflects the rate of phase cancellation, not a local boundary-layer geometry.

---

## 4. Forward Mandate for Cell 118

To decode the global cancellation mechanism:

1. **Dual-Coordinate Profiling ($x = k/N$ vs $j = N - k$):**  
   Analyze the modal flux $F_{N, k} = H_{Nk} k^2 v_k$ in the continuous macroscopic coordinate $x = k/N \in (0, 1)$.
2. **Limiting Cancellation Profile:**  
   Compute the cumulative profile $S_N(x) \equiv \sum_{k \le xN} F_{N, k}$ across $x \in (0, 1)$ for $N \in \{64, 96, 128, 160, 192\}$ to test whether $S_N(x)$ collapses onto a universal continuum cancellation curve $G(x)$.
3. **Arithmetic Correlation Test:**  
   Test whether the residual $\alpha_N$ or $(H u_N)_N$ correlates with the arithmetic symbol fluctuations $\psi'(N)$ or prime-power sums $\sum_{p^k \le c} \frac{\log p}{p^{k/2}} \cos(2\pi N \frac{\log p^k}{L})$.
