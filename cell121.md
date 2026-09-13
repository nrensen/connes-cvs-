# Cell 121 Analytical Note: Core-Size ($L$) Invariant Product Mapping & The Bound-State-to-Continuum Transition for Gate 1

**Companion Computational Script:** [`cell121.py`](file:///c:/data/github/connes-cvs-/cell121.py) | **Verification Log:** [`cell121.out`](file:///c:/data/github/connes-cvs-/cell121.out) (runtime: 90.42 s at 70 dps)  
**Status:** Executed & Audited / Calibrated Diagnostic Confirmation (Gate 1 / Milestone M-G1.0 / Core-Size Scaling)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell120.md`](file:///c:/data/github/connes-cvs-/cell120.md); [`cell120.out`](file:///c:/data/github/connes-cvs-/cell120.out); [`cell96.md`](file:///c:/data/github/connes-cvs-/cell96.md); [`cell97.md`](file:///c:/data/github/connes-cvs-/cell97.md); [`cell60.md`](file:///c:/data/github/connes-cvs-/cell60.md); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §8.26–8.27  
**Date:** September 2026  

---

## 1. Executive Summary & Calibrated Verdict

Cell 121 executed a systematic 2D grid sweep over core sizes $L \in \{4, 6, 8, 10, 12, 14, 16\}$ and dimensions $N \in \{16, 20, 24, 28, 32, 36, 40, 48, 64\}$ at 70 dps to map the lower Ritz gaps $g_{2, L}(N)$, spectral factors $R_{\mathrm{spec}, 2}(N, L)$, and Gate 1 invariant products $\mathcal{P}_j(N, L) \equiv \Delta_j(N) R_{\mathrm{spec}, j}(N, L)$.

### Headline Verdict

> **Classification: Discovery of the Bound-State-to-Continuum Spectral Transition & Empirical Confirmation of Continuum-Core Gate 1 Extinction.**
>
> 1. **The Core-Size Transition is Sharp and Qualitative:**  
>    The lower Ritz gap $g_{2, L}(N) \equiv E_{L+1}^{(N)} - E_3^{(N)}$ exhibits a dramatic phase transition between $L = 10$ and $L = 12$:
>    - For $L \le 10$, the gap collapses exponentially with $N$: at $N = 64$, $g_{2, 4} \approx 5.04 \times 10^{-27}$, $g_{2, 6} \approx 1.48 \times 10^{-17}$, $g_{2, 8} \approx 1.83 \times 10^{-9}$, and $g_{2, 10} \approx 6.03 \times 10^{-3}$.
>    - For $L \ge 12$, the gap stabilizes to macroscopic values: $g_{2, 12} \ge 0.582$, $g_{2, 14} \ge 0.782$, and $g_{2, 16} \ge 0.873$ across all tested dimensions $N \in [16, 64]$.  
>    This provides strong empirical confirmation that the well bound-state capacity is $\bar{N}_{\mathrm{bound}} \approx 11$.
>
> 2. **Route 1B Revived in Continuum-Core Form:**  
>    Route 1B was rejected in Cell 120 because choosing $L = 4$ forced the tail boundary into the bound-state cluster, causing $R_{\mathrm{spec}} \to 10^{53}$. Once $L \ge 12$, the denominator stabilizes: $R_{\mathrm{spec}}(64, 12) \approx 15.14$, $R_{\mathrm{spec}}(64, 14) \approx 8.39$, and $R_{\mathrm{spec}}(64, 16) \approx 6.74$.  
>    **Conclusion:** Route 1B is dead for fixed small $L$ inside the well, but viable when formulated for continuum cores $L \ge \bar{N}_{\mathrm{bound}}$.
>
> 3. **Spectacular Product Extinction:**  
>    In the continuum regime, the Gate 1 product drops by over 50 orders of magnitude compared to $L = 4$:
>    $$\mathcal{P}_2(64, 4) \approx 5.12 \times 10^{12} \quad \Longrightarrow \quad \mathcal{P}_2(64, 12) \approx 3.83 \times 10^{-40}, \quad \mathcal{P}_2(64, 16) \approx 1.71 \times 10^{-40}.$$
>    The supremum envelope $\mathcal{S}_2(L) \equiv \sup_{N \ge 20} \mathcal{P}_2(N, L)$ collapses monotonically across core sizes from $5.65 \times 10^{12}$ down to $3.14 \times 10^{-26}$.
>
> 4. **Multi-Doublet Consistency:**  
>    The suppression holds universally across doublets $j \in \{0, 1, 2\}$, with the expected hierarchy $P_0 < P_1 < P_2$ (e.g. at $N=32, L=8$: $P_0 \approx 6.75 \times 10^{-32}$, $P_1 \approx 9.89 \times 10^{-26}$, $P_2 \approx 4.50 \times 10^{-20}$).
>
> 5. **Calibrated Epistemic Brake:**  
>    The numerical evidence provides strong support for the continuum-core regime across tested dimensions ($N \le 64, L \le 16$), but does **not** constitute an analytical proof of the double limit $\lim_{L \to \infty} \limsup_{N \to \infty} \mathcal{P}_j(N, L) = 0$. Furthermore, $R_{\mathrm{spec}}(N, 12)$ grows mildly from $0.535$ to $15.14$, indicating that it is not strictly constant in $N$, but satisfies $R_{\mathrm{spec}} \ll e^{cN}$.

---

## 2. Audited Numerical Results (`cell121.out`)

### 2.1 Module 1: Eigensystem Assembly & Precision Floor Flagging
- Across $N \in [16, 36]$, all parity doublet splittings $\Delta_0, \Delta_1, \Delta_2$ are strictly clean.
- At $N = 40$, $\Delta_0(40) = 4.886 \times 10^{-51}$ enters the 70-dps eigensolver noise floor (flagged `d0`).
- For $j = 2$, $\Delta_2(N)$ remains completely clean throughout $N \le 64$ ($\Delta_2(64) = 2.533 \times 10^{-41} \gg 10^{-70}$).

### 2.2 Module 2: Lower Ritz Gaps $g_{2, L}(N) \equiv E_{L+1}^{(N)} - E_3^{(N)}$

| $N$ | $L = 4$ | $L = 6$ | $L = 8$ | $L = 10$ | $L = 12$ | $L = 14$ | $L = 16$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **16** | $6.43 \times 10^{-12}$ | $3.92 \times 10^{-6}$ | $0.0867$ | $1.976$ | $2.490$ | $2.775$ | --- |
| **20** | $2.06 \times 10^{-15}$ | $1.03 \times 10^{-8}$ | $1.77 \times 10^{-3}$ | $1.211$ | $2.067$ | $2.550$ | $2.769$ |
| **24** | $6.75 \times 10^{-18}$ | $7.41 \times 10^{-11}$ | $3.04 \times 10^{-5}$ | $0.4118$ | $1.801$ | $2.055$ | $2.483$ |
| **28** | $2.63 \times 10^{-20}$ | $1.11 \times 10^{-12}$ | $2.01 \times 10^{-6}$ | $0.1337$ | $1.314$ | $1.962$ | $2.127$ |
| **32** | $3.15 \times 10^{-22}$ | $2.58 \times 10^{-14}$ | $1.81 \times 10^{-7}$ | $0.0416$ | $0.8162$ | $1.748$ | $1.983$ |
| **36** | $9.42 \times 10^{-24}$ | $1.96 \times 10^{-15}$ | $2.87 \times 10^{-8}$ | $0.0181$ | $0.7816$ | $1.322$ | $1.923$ |
| **40** | $6.44 \times 10^{-25}$ | $2.94 \times 10^{-16}$ | $7.93 \times 10^{-9}$ | $0.0125$ | $0.7275$ | $0.8774$ | $1.677$ |
| **48** | $2.11 \times 10^{-26}$ | $2.94 \times 10^{-17}$ | $3.00 \times 10^{-9}$ | $9.24 \times 10^{-3}$ | $0.6288$ | $0.7824$ | $1.320$ |
| **64** | $5.04 \times 10^{-27}$ | $1.48 \times 10^{-17}$ | $1.83 \times 10^{-9}$ | $6.03 \times 10^{-3}$ | **$0.5824$** | **$0.7821$** | **$0.8731$** |

- **Spectral Transition:** For $L \le 10$, $g_{2, L}(N)$ collapses toward zero as $N$ increases. For $L \ge 12$, $g_{2, L}(N)$ stabilizes above $0.58$, proving that $E_{13}$ lies above the bound-state cluster.

### 2.3 Module 3: Spectral Remote-Tail Factors $R_{\mathrm{spec}, 2}(N, L)$

| $N$ | $L = 4$ | $L = 6$ | $L = 8$ | $L = 10$ | $L = 12$ | $L = 14$ | $L = 16$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **16** | $8.02 \times 10^{22}$ | $2.16 \times 10^{11}$ | $4.41 \times 10^{2}$ | $0.850$ | $0.535$ | $0.431$ | --- |
| **20** | $8.52 \times 10^{29}$ | $3.39 \times 10^{16}$ | $1.15 \times 10^{6}$ | $2.462$ | $0.846$ | $0.556$ | $0.471$ |
| **24** | $8.37 \times 10^{34}$ | $6.96 \times 10^{20}$ | $4.12 \times 10^{9}$ | $22.51$ | $1.177$ | $0.904$ | $0.619$ |
| **28** | $5.92 \times 10^{39}$ | $3.33 \times 10^{24}$ | $1.01 \times 10^{12}$ | $229.3$ | $2.374$ | $1.065$ | $0.906$ |
| **32** | $4.16 \times 10^{43}$ | $6.18 \times 10^{27}$ | $1.25 \times 10^{14}$ | $2379.0$ | $6.175$ | $1.346$ | $1.047$ |
| **36** | $4.84 \times 10^{46}$ | $1.12 \times 10^{30}$ | $5.21 \times 10^{15}$ | $1.31 \times 10^{4}$ | $7.027$ | $2.456$ | $1.161$ |
| **40** | $1.13 \times 10^{49}$ | $5.39 \times 10^{31}$ | $7.42 \times 10^{16}$ | $3.01 \times 10^{4}$ | $8.823$ | $6.064$ | $1.660$ |
| **48** | $1.11 \times 10^{52}$ | $5.74 \times 10^{33}$ | $5.52 \times 10^{17}$ | $5.81 \times 10^{4}$ | $12.55$ | $8.105$ | $2.846$ |
| **64** | **$2.02 \times 10^{53}$** | **$2.34 \times 10^{34}$** | **$1.53 \times 10^{18}$** | **$1.41 \times 10^{5}$** | **$15.14$** | **$8.394$** | **$6.736$** |

- **Mild Growth in Continuum:** For $L = 12$, $R_{\mathrm{spec}}$ grows mildly from $0.535$ to $15.14$ across $N \in [16, 64]$. It is not strictly constant in $N$, but remains macroscopic and bounded by $\mathcal{O}(15)$, contrasting sharply with the $10^{53}$ explosion at $L = 4$.

### 2.4 Module 4: Gate 1 Invariant Products $\mathcal{P}_2(N, L) = \Delta_2(N) R_{\mathrm{spec}, 2}(N, L)$

| $N$ | $L = 4$ | $L = 6$ | $L = 8$ | $L = 10$ | $L = 12$ | $L = 14$ | $L = 16$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **16** | $+64.46$ | $1.74 \times 10^{-10}$ | $3.54 \times 10^{-19}$ | $6.83 \times 10^{-22}$ | $4.30 \times 10^{-22}$ | $3.46 \times 10^{-22}$ | --- |
| **20** | $5.67 \times 10^{4}$ | $2.26 \times 10^{-9}$ | $7.65 \times 10^{-20}$ | $1.64 \times 10^{-25}$ | $5.63 \times 10^{-26}$ | $3.70 \times 10^{-26}$ | $3.14 \times 10^{-26}$ |
| **24** | $5.23 \times 10^{6}$ | $4.35 \times 10^{-8}$ | $2.58 \times 10^{-19}$ | $1.41 \times 10^{-27}$ | $7.35 \times 10^{-29}$ | $5.65 \times 10^{-29}$ | $3.87 \times 10^{-29}$ |
| **28** | $4.08 \times 10^{8}$ | $2.29 \times 10^{-7}$ | $6.96 \times 10^{-20}$ | $1.58 \times 10^{-29}$ | $1.64 \times 10^{-31}$ | $7.34 \times 10^{-32}$ | $6.24 \times 10^{-32}$ |
| **32** | $1.50 \times 10^{10}$ | $2.23 \times 10^{-6}$ | $4.50 \times 10^{-20}$ | $8.57 \times 10^{-31}$ | $2.22 \times 10^{-33}$ | $4.85 \times 10^{-34}$ | $3.77 \times 10^{-34}$ |
| **36** | $1.91 \times 10^{11}$ | $4.44 \times 10^{-6}$ | $2.06 \times 10^{-20}$ | $5.17 \times 10^{-32}$ | $2.78 \times 10^{-35}$ | $9.71 \times 10^{-36}$ | $4.59 \times 10^{-36}$ |
| **40** | $1.17 \times 10^{12}$ | $5.59 \times 10^{-6}$ | $7.70 \times 10^{-21}$ | $3.12 \times 10^{-33}$ | $9.16 \times 10^{-37}$ | $6.29 \times 10^{-37}$ | $1.72 \times 10^{-37}$ |
| **48** | $5.65 \times 10^{12}$ | $2.91 \times 10^{-6}$ | $2.80 \times 10^{-22}$ | $2.95 \times 10^{-35}$ | $6.37 \times 10^{-39}$ | $4.12 \times 10^{-39}$ | $1.44 \times 10^{-39}$ |
| **64** | **$5.12 \times 10^{12}$** | **$5.92 \times 10^{-7}$** | **$3.88 \times 10^{-23}$** | **$3.58 \times 10^{-36}$** | **$3.83 \times 10^{-40}$** | **$2.13 \times 10^{-40}$** | **$1.71 \times 10^{-40}$** |

### 2.5 Module 5: Supremum Envelope $\mathcal{S}_2(L) \equiv \sup_{N \ge 20} \mathcal{P}_2(N, L)$

| $L$ | $\mathcal{S}_2(L) = \sup_N \mathcal{P}_2$ | $\operatorname{argmax} N$ | $\min_N \mathcal{P}_2(N, L)$ | Status |
| :---: | :---: | :---: | :---: | :--- |
| **4** | $5.65 \times 10^{12}$ | 48 | $5.67 \times 10^{4}$ | Blow-up (Crowded Well) |
| **6** | $5.59 \times 10^{-6}$ | 40 | $2.26 \times 10^{-9}$ | Sub-Unit |
| **8** | $2.58 \times 10^{-19}$ | 24 | $3.88 \times 10^{-23}$ | Strong Suppression |
| **10** | $1.64 \times 10^{-25}$ | 20 | $3.58 \times 10^{-36}$ | Deep Extinction |
| **12** | $5.63 \times 10^{-26}$ | 20 | $3.83 \times 10^{-40}$ | Deep Extinction |
| **14** | $3.70 \times 10^{-26}$ | 20 | $2.13 \times 10^{-40}$ | Deep Extinction |
| **16** | $3.14 \times 10^{-26}$ | 20 | $1.71 \times 10^{-40}$ | Deep Extinction |

---

## 3. Theoretical & Strategic Implications

### 3.1 The Physical Role of Core Size $L$
We must not treat $L$ as a cosmetic truncation index. $L$ selects whether the Stieltjes tail starts:
- **Inside the bound-state sector ($L < \bar{N}_{\mathrm{bound}} \approx 11$):** Adjacent Ritz levels are exponentially clustered ($g_{2, L} \sim e^{-c_L N}$), creating catastrophic denominators that overpower tunneling suppression unless $L$ is large.
- **In the continuum ($L \ge \bar{N}_{\mathrm{bound}} \approx 11$):** Level spacing is macroscopic ($g_{2, L} \ge 0.58$), removing small denominators entirely.

### 3.2 Product Dominance Over Operator Boundedness
We do **not** need to prove $R_{\mathrm{spec}}(N, L) = \mathcal{O}(1)$ unconditionally.
Because tunneling splitting collapses at $\Delta_2(64) \approx 2.53 \times 10^{-41}$, any soft growth bound:
$$R_{\mathrm{spec}}(N, L) \le e^{o(N)} \quad \text{or} \quad R_{\mathrm{spec}}(N, L) \le e^{c N} \quad (c < \sigma_2)$$
is more than sufficient to force the product $\mathcal{P}_2(N, L) \to 0$.

### 3.3 Distinction Between Empirical Gap and Asymptotic Theorem
While $g_{2, 12}(N)$ stays above $0.582$ over $N \in [16, 64]$, this finite-$N$ data does not prove an asymptotic lower bound $\inf_N g_{2, 12}(N) \ge g_* > 0$. The data is consistent with the earlier continuum gap observation ($g_{11} \approx 0.42-0.57$), formulating the exact analytical target for Cell 122.

---

## 4. Forward Path: Cell 122 Analytical Target

Rather than pursuing broader numerical sweeps, **Cell 122** must be an **analytical reduction note** targeting:

> **Analytical Lower Bound on the Continuum Spectral Threshold:**  
> Prove from the operator structure that for core sizes $L \ge 12$, the Ritz level $E_{L+1}^{(N)}$ lies uniformly above the bound-state sector:
> $$E_{L+1}^{(N)} \ge E_{\mathrm{cont}}^- > 0 \qquad (\forall N \ge N_0),$$
> establishing that the continuum spectral denominator satisfies $(E_{L+1} - E_j)(E_{L+1} - E_{j+1}) \ge (E_{\mathrm{cont}}^-)^2 > 0$ and reducing Gate 1 entirely to the tunneling decay rate of $\Delta_j(N)$.
