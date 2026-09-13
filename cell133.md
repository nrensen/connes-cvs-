# CELL 133 — COORDINATE-SPACE UNIFICATION OF COUPLED CANCELLATION & THE EXCESS POSITIVITY MARGIN ON $\Phi^\perp$

**Target Gate:** Gate 1 (Finite-$N$ Spectral Mechanism & Asymptotic Tail Extinction, Milestone M-G1.6)  
**Target Proposition:** Structural and Coordinate-Space Origin of the Macroscopic Positivity Margin:
$$R_{\mathrm{net}}(w_{\mathrm{bad}}) \equiv R_{\mathrm{comp}}(w_{\mathrm{bad}}) + R_{\mathrm{arch}}(w_{\mathrm{bad}}) + R_{\mathrm{pole}}(w_{\mathrm{bad}}) \approx -0.4870 + 0.4398 + 0.2121 = +0.1649 > 0$$
on the continuum subspace $\Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}})$, demonstrating that:
1. The discrete competition deficit $R_{\mathrm{comp}} \approx -0.487$ is the variational minimum of an indefinite functional whose negative potential $\frac{1}{L}\int_0^L W(t) |T_v(t)|^2 dt$ is intrinsically constrained by the orthogonality of $v$ to the 11 bound states.
2. The restoring terms $R_{\mathrm{arch}} \approx +0.440$ and $R_{\mathrm{pole}} \approx +0.212$ are direct physical consequences of the oscillatory Fourier localization ($m^* = 26$) and spatial envelope of the physical wavepacket $T_{v_{\mathrm{bad}}}(t)$.
3. Together, they strictly overpower the competition deficit by a macroscopic margin $R_{\mathrm{net}} \approx +0.165 > 0$, explaining the emergent positivity of the Friedrichs continuum operator.

**Verification / Falsification Criteria:**
1. **Continuous Coordinate Integration Residual Audit:** Certify dynamically to 40 decimal digits that:
   $$\big| v_{\mathrm{bad}}^T \widetilde{W} v_{\mathrm{bad}} - \frac{1}{L}\int_0^L W(t) |T_{v_{\mathrm{bad}}}(t)|^2 dt \big| < 10^{-40}$$
   via piecewise high-order Gauss–Legendre quadrature on prime intervals $[\log q_k, \log q_{k+1}]$, using the exact Cell 129 step potential $W(t) = 2 \sum_{q \le c} w_q \mathbf{1}_{[\log q, L]}(t)$.
2. **Archimedean Spectral Partition Audit:** Verify that the continuous Fourier transform $\Phi_{v_{\mathrm{bad}}}(r) = \sum_m (v_{\mathrm{bad}})_m \phi_m(r)$ has a suppressed negative Archimedean component:
   $$\frac{|I_{\mathrm{arch}}^{(-)}|}{I_{\mathrm{arch}}^{\mathrm{tot}}} \approx 5.8\%, \qquad I_{\mathrm{arch}}^{(-)} = \frac{1}{\pi} \int_0^{r_*} |\Phi_{v_{\mathrm{bad}}}(r)|^2 h_+(r) dr \approx -0.0935,$$
   and that the continuous total integral matches $v_{\mathrm{bad}}^T Q_{\mathrm{even}}^{\mathrm{arch}} v_{\mathrm{bad}}$ to high precision.
3. **Asymptotic Margin Stabilization:** Track the net surplus margin:
   $$\Delta_{\mathrm{surplus}}(N) \equiv (R_{\mathrm{arch}} + R_{\mathrm{pole}}) - |R_{\mathrm{comp}}| = R_{\mathrm{net}}(w_{\mathrm{bad}})$$
   across $N \in [24, 64]$ and demonstrate that it stabilizes to a strictly positive macroscopic constant $\approx 0.165 > 0$.

**Companion Computational Script:** [`cell133.py`](file:///c:/data/github/connes-cvs-/cell133.py)  

---

## 1. Executive Context & The Mystery of the Excess Margin

In Cell 131, the hypothesis of independent operator domination was decisively falsified:
$$\lambda_{\min}(\widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}) \to -0.486979 < 0.$$
Yet the full compressed operator remains strictly positive:
$$\lambda_{\min}(\widehat{Q}_{\mathrm{even}}) = E_{11}(N) = 0.00602464 > 0 \quad (N = 64).$$

In Cell 132, the exact operator decomposition was audited and certified to an operator norm residual of $\|\mathcal{R}\|_{\mathrm{op}} \le 1.0319 \times 10^{-48} < 10^{-45}$. Cell 132 revealed that:
1. **Strict Solitary Instability:** The competition operator $\widehat{\mathcal{Q}}_{\mathrm{comp}} \equiv \widehat{\Omega} - \widehat{\mathcal{K}}_{\mathrm{neg}}$ possesses **precisely one** negative eigenvalue across all tested dimensions $N \ge 24$ ($k_{\mathrm{neg}} \equiv 1$). All remaining 53 directions on $\Phi^\perp$ (at $N = 64$) are strictly positive ($\mu_1 \approx +0.3483 > 0$).
2. **The Decisive Restoration:** Along the solitary negative direction $w_{\mathrm{bad}}$, the deficit is $R_{\mathrm{comp}} \approx -0.4870$. But the off-diagonal Archimedean remainder supplies $R_{\mathrm{arch}} \approx +0.4398$, and the pole projector supplies $R_{\mathrm{pole}} \approx +0.2121$.
3. **The Net Margin:**
   $$R_{\mathrm{net}}(w_{\mathrm{bad}}) = -0.4870 + 0.4398 + 0.2121 = +0.1649 > 0.$$
   The restoring efficiency is $\rho_{\mathrm{restore}} \approx 1.3387 > 1$, leaving an excess margin of $+0.165$.
4. **Asymptotic Freeze:** As $N \to 64$, the net margin stabilizes at $\approx 0.165$, and the Fourier profile of $v_{\mathrm{bad}} = U_{\mathrm{cont}} w_{\mathrm{bad}}$ freezes with its dominant node at $m^* = 26$.

The central mathematical question of **Cell 133** is:
> *Why does the sum of the restoring terms ($+0.652$) systematically overpower the solitary competition deficit ($|-0.487|$) by a stable macroscopic margin ($+0.165$)? What physical mechanics in coordinate and frequency space govern this balance?*

---

## 2. Continuous Coordinate Wavepacket Representation of $v_{\mathrm{bad}}$

Every discrete coefficient vector $v \in \mathbb{R}^{N+1}$ in the canonical even basis corresponds to a unique continuous wavefunction $T_v(t)$ on the interval $[0, L]$:
$$T_v(t) = v_0 \phi_0(t) + \sum_{m=1}^N v_m \phi_m(t),$$
where the orthonormal basis under the normalized measure $d\mu = \frac{dt}{L}$ is:
$$\phi_0(t) \equiv 1, \qquad \phi_m(t) \equiv \sqrt{2} \cos\left(\frac{2\pi m t}{L}\right) \quad (m \ge 1).$$

### 2.1 Orthonormality and Isometric Parseval Identity
Under the inner product $\langle f, g \rangle_{L^2(d\mu)} \equiv \frac{1}{L} \int_0^L f(t) g(t) dt$:
$$\langle \phi_m, \phi_n \rangle_{L^2(d\mu)} = \delta_{mn}.$$
Consequently, the $L^2(d\mu)$ norm of $T_v(t)$ coincides identically with the $\ell^2$ Euclidean norm of $v$:
$$\|T_v\|_{L^2(d\mu)}^2 = \frac{1}{L} \int_0^L |T_v(t)|^2 dt = v_0^2 + \sum_{m=1}^N v_m^2 = \|v\|_2^2 = 1.$$

### 2.2 Boundary Values and Midpoint Symmetries
Because $\cos(0) = \cos(2\pi m) = 1$ for all integers $m$:
$$T_v(0) = T_v(L) = v_0 + \sqrt{2} \sum_{m=1}^N v_m.$$
At the midpoint $t = L/2$, $\cos(\pi m) = (-1)^m$:
$$T_v(L/2) = v_0 + \sqrt{2} \sum_{m=1}^N (-1)^m v_m.$$

### 2.3 Spatial Kinetic Semi-Norm
The spatial derivative is:
$$T'_v(t) = -\frac{2\pi \sqrt{2}}{L} \sum_{m=1}^N m v_m \sin\left(\frac{2\pi m t}{L}\right).$$
Its $L^2(d\mu)$ semi-norm represents the physical kinetic energy:
$$\mathcal{K}_{\mathrm{kin}}[T_v] \equiv \frac{1}{L} \int_0^L |T'_v(t)|^2 dt = \left(\frac{2\pi}{L}\right)^2 \sum_{m=1}^N m^2 v_m^2.$$

For $v_{\mathrm{bad}}$, the modal energy is concentrated around the frozen peak mode $m^* = 26$. With $L = \log 13 \approx 2.56495$, the scale factor is:
$$\left(\frac{2\pi}{L}\right)^2 \approx 5.9996 \approx 6.00.$$
A modal concentration at $m^* = 26$ carries a large kinetic energy density:
$$\mathcal{K}_{\mathrm{kin}} \approx 6.00 \times 26^2 \times 0.120 \approx 487 \gg 1.$$
This demonstrates that $T_{v_{\mathrm{bad}}}(t)$ is an oscillatory wavepacket with characteristic spatial wavelength:
$$\lambda^* = \frac{L}{m^*} = \frac{2.56495}{26} \approx 0.09865.$$

---

## 3. The Step-Potential Well and the Bound-State Exclusion Obstruction

### 3.1 The Physical Prime Step Potential $W(t)$
Recall from Cell 129 and Cell 130 that the prime operator decomposes as:
$$Q_{\mathrm{prime}}^{\mathrm{even}} \equiv -\widetilde{W} + D_{\mathrm{per}} + \Delta \widetilde{\mathcal{D}},$$
where $\widetilde{W}$ is the Galerkin projection of the physical step-potential operator defined in Cell 129 Definition 5.1:
$$\boxed{W(t) \equiv 2 \sum_{\substack{q \le c \\ \log q \le t}} w_q = 2 \sum_{q \le c} \frac{\Lambda(q)}{\sqrt{q}} \cdot \mathbf{1}_{[\log q, L]}(t) \ge 0.}$$

The factor of 2 originates from the symmetrized autocorrelation in Weil's explicit formula and the Connes–CvS Galerkin matrix:
$$\langle v, Q_{\mathrm{prime}} v \rangle = - \frac{2}{L} \sum_{q \le c} w_q \int_{\log q}^L T_v(t) T_v(t - \log q) \, dt.$$
Under the algebraic polarization $2 T_v(t) T_v(t - \log q) = T_v(t)^2 + T_v(t - \log q)^2 - |T_v(t) - T_v(t - \log q)|^2$, the local component contributes $2 \int_{\log q}^L T_v(t)^2 dt$ by midpoint reflection symmetry. Hence the physical step potential multiplying $T_v(t)^2 \frac{dt}{L}$ has weight $2 w_q$.

For $c = 13$, the primes and prime powers $q = p^k \le 13$ are:
- $q = 2$: $\log 2 \approx 0.69315$, weight $w_2 = \frac{\log 2}{\sqrt{2}} \approx 0.49013$, step increment $2 w_2 \approx 0.98026$
- $q = 3$: $\log 3 \approx 1.09861$, weight $w_3 = \frac{\log 3}{\sqrt{3}} \approx 0.63424$, step increment $2 w_3 \approx 1.26848$
- $q = 4$: $\log 4 \approx 1.38629$, weight $w_4 = \frac{\log 2}{2} \approx 0.34657$, step increment $2 w_4 \approx 0.69315$
- $q = 5$: $\log 5 \approx 1.60944$, weight $w_5 = \frac{\log 5}{\sqrt{5}} \approx 0.71977$, step increment $2 w_5 \approx 1.43954$
- $q = 7$: $\log 7 \approx 1.94591$, weight $w_7 = \frac{\log 7}{\sqrt{7}} \approx 0.73546$, step increment $2 w_7 \approx 1.47092$
- $q = 8$: $\log 8 \approx 2.07944$, weight $w_8 = \frac{\log 2}{\sqrt{8}} \approx 0.24506$, step increment $2 w_8 \approx 0.49013$
- $q = 9$: $\log 9 \approx 2.19722$, weight $w_9 = \frac{\log 3}{3} \approx 0.36620$, step increment $2 w_9 \approx 0.73239$
- $q = 11$: $\log 11 \approx 2.39790$, weight $w_{11} = \frac{\log 11}{\sqrt{11}} \approx 0.72300$, step increment $2 w_{11} \approx 1.44600$
- $q = 13$: $\log 13 \approx 2.56495 = L$, weight $w_{13} = \frac{\log 13}{\sqrt{13}} \approx 0.71138$, step increment $2 w_{13} \approx 1.42277$

The profile of $W(t)$ features:
1. **A Completely Flat Plateau on $[0, \log 2)$:** Because the smallest prime is $p = 2$, no prime power satisfies $\log q < \log 2$. Hence:
   $$W(t) \equiv 0 \qquad \forall t \in [0, \log 2) \approx [0, 0.69315).$$
   This flat zone occupies $\frac{\log 2}{\log 13} \approx 27.02\%$ of the total interval $[0, L]$.
2. **A Deepening Well on $[\log 2, L]$:** At each $\log q$, $W(t)$ steps up by $2 w_q$. At the right endpoint:
   $$W(L) = 2 \sum_{q \le 13} w_q \approx 8.9212.$$
   The average potential value is $-\psi_{\mathrm{prime}}'(0) = \frac{1}{L}\int_0^L W(t) dt = 2 \sum_{q \le c} w_q (1 - \log q / L) \approx 5.316$.

### 3.2 Exact Continuous Integration Identity
**Theorem 133.1 (Coordinate-Matrix Equivalence):**
For any $v \in \mathbb{R}^{N+1}$, the matrix quadratic form evaluates identically to the spatial continuous integral:
$$\boxed{v^T \widetilde{W} v = \frac{1}{L} \int_0^L W(t) |T_v(t)|^2 dt.}$$

*Proof:* By definition of the matrix elements $\widetilde{W}_{mn} = \frac{1}{L} \int_0^L W(t) \phi_m(t) \phi_n(t) dt$ (derived in Cell 129 via exact Fourier moments), summing over $m, n \in \{0, \dots, N\}$ yields:
$$\sum_{m, n=0}^N v_m v_n \widetilde{W}_{mn} = \frac{1}{L} \int_0^L W(t) \left( \sum_{m=0}^N v_m \phi_m(t) \right) \left( \sum_{n=0}^N v_n \phi_n(t) \right) dt = \frac{1}{L} \int_0^L W(t) |T_v(t)|^2 dt. \quad \blacksquare$$

### 3.3 The Bound-State Exclusion Obstruction
Why does the potential energy $\langle v, \widetilde{W} v \rangle$ fail to overpower the kinetic/translation stiffness $\Omega$?
- If an unconstrained wavefunction $f \in L^2$ were placed on $[0, L]$, it could concentrate $100\%$ of its mass inside the deep well region $[\log 11, L]$, achieving an energy of $\approx 8.92$.
- However, $v_{\mathrm{bad}} \in \Phi^\perp = \operatorname{Ran}(I - P_{\mathrm{bound}})$, meaning $v_{\mathrm{bad}}$ is **strictly orthogonal to the 11 lowest bound states** $u_0, u_1, \dots, u_{10}$ of $Q_{\mathrm{even}}$.
- Orthogonality to the bound-state spectral subspace forces the state out of the low-energy eigenspace of $Q_{\mathrm{even}}$. Numerically, the principal competition direction $w_{\mathrm{bad}}$ is predominantly represented by high Fourier modes ($59.3\%$ of energy in modes $m \ge 11$, with peak mode at $m^* = 26$).
- While low modes do not vanish entirely ($|v_0|^2 \approx 0.0286$, modes $1:3 \approx 17.0\%$, modes $4:10 \approx 20.8\%$), this high-frequency distribution spreads mass across the entire interval $[0, L]$.
- Consequently, $T_{v_{\mathrm{bad}}}(t)$ cannot concentrate purely in the deep well $[\log 2, L]$:
  $$\mathcal{M}_{\mathrm{flat}}[v] \equiv \frac{1}{L} \int_0^{\log 2} |T_v(t)|^2 dt \approx 0.3938 \quad (39.4\%), \qquad \mathcal{M}_{\mathrm{well}}[v] \equiv \frac{1}{L} \int_{\log 2}^L |T_v(t)|^2 dt \approx 0.6062 \quad (60.6\%).$$
  Because nearly $40\%$ of the probability mass remains in the zero-potential plateau $[0, \log 2)$ (where $W(t) \equiv 0$), the potential energy harvested by $T_{v_{\mathrm{bad}}}$ is substantially diluted:
  $$v_{\mathrm{bad}}^T \widetilde{W} v_{\mathrm{bad}} \approx 3.710 \ll W(L) \approx 8.921.$$

---

## 4. Archimedean Spectral Partitioning and Sidelobe Transfer

### 4.1 Continuous Fourier Transform $\Phi_v(r)$
In continuous Fourier space, the wavepacket $T_v(t) \cdot \mathbf{1}_{[0, L]}(t)$ has Fourier transform:
$$\Phi_v(r) \equiv \sum_{m=0}^N v_m \phi_m(r),$$
where $\phi_m(r)$ are the continuous Fourier basis amplitudes:
$$\phi_0(r) = \frac{2}{\sqrt{L}} \frac{\sin(r L / 2)}{r}, \qquad \phi_m(r) = \frac{2\sqrt{2}}{\sqrt{L}} \frac{r \sin(r L / 2)}{r^2 - a_m^2} \cdot (-1)^m \quad (m \ge 1).$$

The Archimedean quadratic form evaluates to:
$$\langle v, Q_{\mathrm{even}}^{\mathrm{arch}} v \rangle = \frac{1}{\pi} \int_0^\infty |\Phi_v(r)|^2 h_+(r) dr.$$

### 4.2 The Archimedean Sign Problem on $\Phi^\perp$
The multiplier $h_+(r) = \operatorname{Re}\psi(1/4 + ir/2) - \log \pi$ is strictly negative on $[0, r_*]$ ($r_* \approx 6.28984$, with $h_+(0) \approx -5.37218$).
We partition the continuous integral into:
$$I_{\mathrm{arch}}^{(-)} \equiv \frac{1}{\pi} \int_0^{r_*} |\Phi_{v_{\mathrm{bad}}}(r)|^2 h_+(r) dr < 0,$$
$$I_{\mathrm{arch}}^{(+)} \equiv \frac{1}{\pi} \int_{r_*}^\infty |\Phi_{v_{\mathrm{bad}}}(r)|^2 h_+(r) dr > 0.$$

**Proposition 133.2 (Suppression of the Negative Archimedean Contribution):**
On the continuum subspace $\Phi^\perp$, the negative Archimedean contribution evaluates to:
$$I_{\mathrm{arch}}^{(-)} \approx -0.0935 \quad (N = 64).$$
While not negligible on the scale of the net margin $+0.165$ (it represents roughly half the surplus), it is **significantly suppressed relative to the total Archimedean energy**:
$$\frac{|I_{\mathrm{arch}}^{(-)}|}{I_{\mathrm{arch}}^{\mathrm{tot}}} = \frac{0.0935}{1.6136} \approx 5.8\%.$$
The dominant contribution to the Archimedean form comes from the positive spectrum ($I_{\mathrm{arch}}^{(+)} \approx 1.7071$).

### 4.3 Origin of the Positive Off-Diagonal Excess $R_{\mathrm{arch}} \approx +0.440$
The discrete diagonal Archimedean energy is:
$$D_{\mathrm{arch}} \equiv \sum_{m=0}^N v_m^2 h_+\left(\frac{2\pi m}{L}\right) \approx 1.1738 \quad (N = 64).$$
The off-diagonal divided-difference contribution is:
$$R_{\mathrm{arch}} \equiv \langle v, \Delta Q_{\mathrm{arch}} v \rangle = \langle v, Q_{\mathrm{arch}} v \rangle - D_{\mathrm{arch}} = I_{\mathrm{arch}}^{\mathrm{tot}} - D_{\mathrm{arch}} = 1.6136 - 1.1738 = +0.4398.$$
The numerical decomposition rigorously verifies this value. The proposed physical hypothesis is that windowing to $[0, L]$ generates broad continuous Fourier sinc sidelobes that sample the logarithmic growth of $h_+(r) \sim \log(r/2)$ at high frequencies, but an analytical proof of this off-diagonal positivity remains an active Gate 1 research question.

---

## 5. The Zeta Pole Form $\mathcal{Q}_{\mathrm{pole}}[T]$

The pole contribution in the Weil explicit formula represents the contribution of the trivial poles of $\zeta(s)$ at $s = 0, 1$:
$$\mathcal{W}_{\mathrm{pole}}[g] = 2 \int_{-\infty}^\infty g(x) \cosh(x/2) dx.$$
In the truncated Galerkin setting on $[0, L]$, the autocorrelation $g(x) = (T_v \star T_v^*)(x)$ generates a positive semidefinite operator $Q_{\mathrm{even}}^{\mathrm{pole}} \succeq 0$.

In physical coordinate space, the overlap functional is:
$$J_{\mathrm{pole}}[T_v] \equiv \frac{1}{L} \int_0^L T_v(t) \cosh\left(\frac{1}{2}\left(t - \frac{L}{2}\right)\right) dt.$$
Because $\cosh(x/2) \ge 1$ is positive and symmetric about the midpoint $t = L/2$, any wavefunction with non-zero average value or even-parity envelope couples strongly to this functional.
In Cell 132, we computed:
$$R_{\mathrm{pole}}(w_{\mathrm{bad}}) \approx +0.212100 > 0.$$
This provides an autonomous, strictly positive stabilizing contribution that is completely decoupled from the competition deficit.

---

## 6. Synthesis: The Variational Origin of the Macroscopic Margin

Combining all four pieces, the net energy functional along any unit state $v \in \Phi^\perp$ is:
$$\mathcal{E}_{\mathrm{net}}[T_v] = \underbrace{\big(\mathcal{E}_{\Omega}[T_v] - \mathcal{E}_{\mathrm{neg}}[T_v]\big)}_{R_{\mathrm{comp}}} + \underbrace{\mathcal{E}_{\mathrm{arch}}^{\mathrm{off}}[T_v]}_{R_{\mathrm{arch}}} + \underbrace{\mathcal{Q}_{\mathrm{pole}}[T_v]}_{R_{\mathrm{pole}}}.$$

Along the solitary negative eigenvector $v_{\mathrm{bad}}$ that minimizes the competition deficit:
1. **The Deficit Sinks to a Floor:** To maximize $\mathcal{E}_{\mathrm{neg}}[T_v] = \frac{1}{L}\int_0^L W(t) |T_v|^2 dt - \langle v, \Delta \widetilde{\mathcal{D}} v \rangle$, $T_v$ attempts to concentrate mass in $[\log 2, L]$. But orthogonality to $u_0, \dots, u_{10}$ forces a high-frequency profile ($m^* = 26$), diluting the harvested potential to $\approx 3.71$ (with $39.4\%$ of mass trapped in the flat plateau $[0, \log 2)$). The kinetic/translation stiffness $\Omega$ stabilizes the minimum at:
   $$R_{\mathrm{comp}} \approx -0.4870.$$
2. **The High-Frequency State Activates the Restoring Channels:** The modal localization produces:
   $$R_{\mathrm{arch}} \approx +0.4398, \qquad R_{\mathrm{pole}} \approx +0.2121.$$
3. **The Macroscopic Surplus:**
   $$\Delta_{\mathrm{surplus}} = R_{\mathrm{arch}} + R_{\mathrm{pole}} - |R_{\mathrm{comp}}| \approx 0.4398 + 0.2121 - 0.4870 = +0.1649 > 0.$$

### Epistemic Assessment & Path to Cell 134
The computational audit rigorously certifies exact finite-$N$ identities (including coordinate-matrix potential equivalence to $< 10^{-40}$) and demonstrates that along the solitary negative direction $w_{\mathrm{bad}}$, the deficit $R_{\mathrm{comp}} \approx -0.487$ is overcome by $R_{\mathrm{arch}} \approx +0.440$ and $R_{\mathrm{pole}} \approx +0.212$, leaving a positive surplus margin $\Delta_{\mathrm{surplus}} = R_{\mathrm{net}} \approx +0.165 > 0$.

The data strongly support the hypothesis that finite-$N$ coupled positivity has a stable coordinate and spectral mechanism governed by a localized spatial wavepacket ($M_{\mathrm{well}} \approx 60.6\%$, $t_{\mathrm{peak}} \approx 2.539$, $m^* = 26$). However, neither the limiting profile nor a uniform positive margin as $N \to \infty$ has yet been proved analytically.

This establishes the analytical agenda for **Cell 134**:
$$\boxed{\textbf{Can the competition minimizer } v_{\mathrm{bad}} \textbf{ be characterized as a constrained variational state for a 1D step potential?}}$$
