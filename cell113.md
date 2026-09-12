# Cell 113 Analytical Note: Eigenvector Overlap Continuation, Finite-T Separation & Solitary Branch Extinction

**Companion Computational Script:** [`cell113.py`](file:///c:/data/github/connes-cvs-/cell113.py) | **Verification Log:** `cell113.out`  
**Status:** Working Research Note (Gate 1 / Milestone M12 / Route D / Phase IX)  
**Target:** [ROADMAP.md Gate 1 — Finite-N Spectral Mechanism & Joint-Limit Tail Extinction](file:///c:/data/github/connes-cvs-/ROADMAP.md)  
**Dependencies:** [`cell112a.md`](file:///c:/data/github/connes-cvs-/cell112a.md); [`cell112a.out`](file:///c:/data/github/connes-cvs-/cell112a.out); [`cell111.md`](file:///c:/data/github/connes-cvs-/cell111.md); [`cell111.out`](file:///c:/data/github/connes-cvs-/cell111.out); [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md) §9.10  
**Date:** September 2026  

---

## Executive Summary

Cell 112a established that the scalar plateau $\alpha_N \approx 1.92 \times 10^{-22}$ and defect plateau $T_{v_N}(0) \approx 6.66 \times 10^{-25}$ observed for $N \ge 64$ in Cell 111 are not eigensolver precision noise. Rather, the spectrum undergoes a structural reordering:
- For $N \le 32$, the literal ground state ($k=0$) is the localized solitary wave ($E_0 \approx 2.17 \times 10^{-49}, v_0 \approx 0.543$).
- For $N \ge 64$, the literal lowest eigenvalue drops below zero ($E_0 = -5.08 \times 10^{-52} \to -1.063 \times 10^{-51}$) and corresponds to a delocalized edge/background mode with $v_0 \approx 0.064$.
- Concurrently, the localized solitary wave persists as the first excited eigenstate ($k=1$), with stable positive energy $E_1 \approx 4.71 \times 10^{-50}$ and dominant central amplitude $v_0(E_1) \approx 0.666$.

### The Epistemic Realignment
1. **The Eigensolver Did Not "Switch Incorrectly":**
   The eigensolver was queried for the smallest eigenvalue and correctly returned it. What changed was our interpretation: the lowest-eigenvalue branch ceases to coincide with the localized solitary-wave branch once the finite-$T$ edge mode drops below it.
2. **Consequence for Gate 1 and Theorem 9.16:**
   Theorem 9.16 (the Non-Circular Regularity Bridge in [`Paper-NR2.md`](file:///c:/data/github/connes-cvs-/Paper-NR2.md)) requires uniform $H^2$ control of the state entering the continuum limit. If the finite-$(N, T)$ ground state is a cutoff-induced leakage mode while the continuum solitary wave is an excited state, Gate 1 must be formulated for the **branch-selected solitary eigenvector** $v_N^{\mathrm{sol}}$, rather than casually identifying it with the literal finite-$N$ ground state.

Cell 113 executes the threefold investigation required to resolve this frontier:
1. **Part A (Overlap Branch Continuation):** Tracks $v_N^{\mathrm{sol}}$ across a fine grid $N \in [24, 96]$ ($\Delta N = 4$) via inter-dimensional eigenvector overlap to determine whether the reordering is an avoided crossing or a level crossing.
2. **Part B (Finite-$T$ Leakage Intervention):** Tests whether $E_{\mathrm{edge}}(T)$ scales with the continuous Archimedean tail defect $\delta_T^{\mathrm{tail}}$ across $T \in \{400, 500, 600\}$ at fixed $N = 48$, while $E_{\mathrm{sol}}(T)$ remains stable.
3. **Part C (Extinction on $v_N^{\mathrm{sol}}$):** Evaluates $T_{\mathrm{sol}}(0)$, $\alpha_N^{\mathrm{sol}}$, and the extinction product $P_\alpha^{\mathrm{sol}}(N) = |\alpha_N^{\mathrm{sol}}| \sqrt{N}$ specifically on the localized solitary branch.

---

## 1. Overlap Branch Continuation (Part A)

### 1.1 Mathematical Formulation
Let $v_{N_{\mathrm{prev}}}^{\mathrm{sol}} \in \mathbb{R}^{N_{\mathrm{prev}}+1}$ be the established solitary wave state at dimension $N_{\mathrm{prev}}$. When stepping to $N > N_{\mathrm{prev}}$, we embed $v_{N_{\mathrm{prev}}}^{\mathrm{sol}}$ into $\mathbb{R}^{N+1}$ by zero-padding:
$$v_{\mathrm{embed}} = \big( v_{N_{\mathrm{prev}}}^{\mathrm{sol}}[0], \dots, v_{N_{\mathrm{prev}}}^{\mathrm{sol}}[N_{\mathrm{prev}}], \underbrace{0, \dots, 0}_{N - N_{\mathrm{prev}}} \big)^T.$$

For each of the lowest $K = 5$ orthonormal eigenvectors $\{v_N^{(0)}, \dots, v_N^{(4)}\}$ of $H_N$, we compute the absolute overlap:
$$\mathcal{O}_k(N) \equiv \big| \langle v_{\mathrm{embed}}, v_N^{(k)} \rangle \big| = \left| \sum_{m=0}^{N_{\mathrm{prev}}} v_{N_{\mathrm{prev}}}^{\mathrm{sol}}[m] v_N^{(k)}[m] \right|.$$

The solitary wave branch index at dimension $N$ is defined uniquely by:
$$\boxed{k_{\mathrm{sol}}(N) \equiv \operatorname{argmax}_{0 \le k < K} \mathcal{O}_k(N).}$$

The tracked solitary eigenvector is $v_N^{\mathrm{sol}} \equiv v_N^{(k_{\mathrm{sol}}(N))}$.

### 1.2 Diagnostic Targets
- **Continuity Threshold:** With a fine step $\Delta N = 4$, if the branch is continuous, the overlap satisfies $\mathcal{O}_{k_{\mathrm{sol}}} > 0.95$.
- **Crossing Anatomy:**
  - If $\min_N |E_1(N) - E_0(N)| > 0$ and the overlap smoothly transfers between $k=0$ and $k=1$, the transition is an **avoided crossing** driven by weak boundary coupling.
  - If $E_0$ and $E_1$ intersect with vanishing gap and orthogonal states, it is an **exact level crossing** permitted by differing spatial or parity symmetries.

---

## 2. Finite-$T$ Separation Hypothesis (Part B)

### 2.1 The Physical Mechanism
In [`Paper-NR1.md`](file:///c:/data/github/connes-cvs-/Paper-NR1.md) (Theorem 5.5), the truncated Galerkin matrix $Q^{(T)}$ satisfies:
$$v^T Q^{(T)} v - \mathcal{Q}_{\mathrm{total}}^{(\infty)}(v) \equiv -\frac{1}{\pi} \int_T^\infty h_+(r) K_{\mathrm{Fourier}}(v, r, L) \, dr \equiv -\delta_T^{\mathrm{tail}}(v).$$

For large cutoff $T > a_N = 2\pi N/L$, the tail defect expands as:
$$\delta_T^{\mathrm{tail}}(v) = \sum_{k=0}^\infty A_k(v, N) \mathcal{J}_k(T, L),$$
where $\mathcal{J}_0(T, L) = \frac{1}{\pi} \int_T^\infty \frac{h_+(r)(1 - \cos rL)}{r^2} dr$.

Because $h_+(r) \sim \log(r / 2\pi)$ for large $r$, $\mathcal{J}_0(T, L)$ is positive and decays with $T$. Consequently, the finite-$T$ truncation subtracts a small positive quantity $-\delta_T^{\mathrm{tail}}$ from the continuum quadratic form.

### 2.2 The Intervention Experiment
At fixed $N = 48$ (the threshold where $E_0 \approx 2.1 \times 10^{-50}$ and $E_1 \approx 1.49 \times 10^{-49}$), we vary $T \in \{400, 500, 600\}$:
- **Hypothesis B1 (Edge Mode Tracking):** The edge eigenvalue $E_{\mathrm{edge}}(T)$ shifts in direct correlation with the Archimedean cutoff tail $\mathcal{J}_0(T, L)$.
- **Hypothesis B2 (Solitary Invariance):** The solitary eigenvalue $E_{\mathrm{sol}}(T)$ is dominated by the localized potential well and remains stable across varying $T$.

---

## 3. Extinction Metrics on the Solitary Branch (Part C)

### 3.1 Boundary Observables
On the certified solitary wave branch $v_N^{\mathrm{sol}}$, we compute:
1. **Boundary Contact Defect:**
   $$T_{\mathrm{sol}}(0) = v_{N, 0}^{\mathrm{sol}} + \sqrt{2} \sum_{m=1}^N v_{N, m}^{\mathrm{sol}}.$$
2. **Boundary Multiplier Sum:**
   $$\alpha_N^{\mathrm{sol}} = \sum_{k=1}^N a_k v_{N, k}^{\mathrm{sol}} = \sum_{k=1}^N 2 k \psi(k) v_{N, k}^{\mathrm{sol}}.$$
3. **Second Curvature Moment:**
   $$S_2^{\mathrm{sol}}(N) = \sum_{k=1}^N k^2 v_{N, k}^{\mathrm{sol}}.$$
4. **Extinction Product:**
   $$P_\alpha^{\mathrm{sol}}(N) \equiv |\alpha_N^{\mathrm{sol}}| \sqrt{N}.$$
5. **Boundary Ratio Calibration:**
   Per the reviewer's calibration, the proportionality is formulated as:
   $$\kappa_\alpha^{\mathrm{sol}}(N) \equiv \frac{|\alpha_N^{\mathrm{sol}}|}{|T_{\mathrm{sol}}(0)|} \approx \sqrt{2} a_N = 2 N \psi(N).$$

### 3.2 Hypothesis C1 (Restoration of Extinction Decay)
The plateau $P_\alpha(N) \approx 2.66 \times 10^{-21}$ observed in Cell 111 is an artifact of measuring the delocalized edge mode ($k=0$ for $N \ge 64$). On the true localized solitary branch $v_N^{\mathrm{sol}}$, $P_\alpha^{\mathrm{sol}}(N) \to 0$ as $N \to \infty$, satisfying hypothesis $(H_{\mathrm{ext}})$.

---

## 4. Pre-Flight Specification for `cell113.py`

| Parameter / Module | Specification | Purpose |
| :--- | :--- | :--- |
| **Precision** | `mp.mp.dps = 70`, `GROUND_DPS = 70` | Leverages cached $N=192$ matrix ($T=600$) in 1.7 s |
| **Fine Sweep Grid** | $N \in \{24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 80, 96\}$ | Resolves transition region with $\Delta N = 4$ |
| **Overlap Algorithm** | Zero-padded dot product $\mathcal{O}_k = \|\langle v_{\mathrm{prev}}^{\mathrm{sol}}, v_N^{(k)} \rangle\|$ | Robust, mathematical branch continuation |
| **Part B Cutoff Grid** | $T \in \{400, 500, 600\}$ at fixed $N = 48$ | Tests edge eigenvalue response to Archimedean tail |
| **Part C Extinction Grid** | $T_{\mathrm{sol}}(0)$, $\alpha_N^{\mathrm{sol}}$, $P_\alpha^{\mathrm{sol}}(N)$, $\|\xi_N^{(Q)}\|_2$ | Measures genuine solitary wave extinction |
| **Output Style** | Dispassionate, dry numeric reports | Eliminates "certified" and interpretive narratives from stdout |
