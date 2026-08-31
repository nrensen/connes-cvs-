# Connes–CvS exploratory cell history

**Repository:** `nrensen/connes-cvs-`  
**Historical snapshot audited:** commit `2275abc0bacbc804a4197ef5c5e65d03aa82312f`  
**Purpose of this document:** preserve a compact map of what each exploratory `cell<n>.py` was intended to investigate, what it established, and where discrepancies or historical errors were encountered.

This document is a **research-history map**, not a claim that every historical cell is mathematically correct. In particular, known errors are intentionally recorded rather than silently corrected.

---

## How to read the map

The cells evolved as an investigation. Later cells often supersede, clarify, or localise questions raised by earlier cells.

The statuses used below are:

- **Established** — the intended question was satisfactorily answered.
- **Diagnostic / superseded** — useful exploratory work, but a later cell provides the cleaner result.
- **Discrepancy found** — the cell exposed a genuine mismatch that required further investigation.
- **Historical error** — the cell contains a mathematical mistake that is intentionally retained as part of the research record.
- **Open / next step** — the investigation was not yet complete.

A particularly important distinction is between **historical calculations** and **current canonical definitions**. Historical cells should remain reproducible even when we now know that a calculation was wrong.

---

# Cells 0–4 — initial reconstruction and dictionary work

## Cell 0 — initial ground-state / zero-side reconstruction

### Intended purpose

Cell 0 was the initial exploratory reconstruction of the repository calculation at the small test case

- `c = 13`
- `N = 8`
- `T = 60`
- high precision.

It independently reconstructed the even-sector Fourier response, checked the first Riemann zero, differentiated the response, reconstructed the canonical coefficient representation, and explored the constraint geometry associated with zero and pole conditions.

### What it established

It established the basic numerical setting for the later investigation:

- construction of the ground state;
- canonical/full coefficient relationship;
- evaluation of `F(gamma_j)` and `F'(gamma_j)`;
- pole functional;
- constrained maximisation/projection calculations.

The cell became an exploratory foundation rather than the final canonical implementation.

### Status

**Diagnostic / superseded.**

---

## Cell 1 — canonical basis and response construction

### Intended purpose

Cell 1 developed the canonical basis-level constructions:

- `F_basis(k, tau)`;
- `Fprime_basis(k, tau)`;
- coefficient-vector response;
- relation between canonical and full Fourier coefficients.

It was an early attempt to make the finite-dimensional Fourier dictionary explicit rather than treating the repository implementation as a black box.

### What it established

The basis-response machinery became important later and was progressively distilled into `cell.py`.

### Status

**Diagnostic / superseded.**

---

## Cell 2 — zero-side Guinand–Weil reconstruction

### Intended purpose

Cell 2 explicitly reconstructed the zero-side finite Guinand–Weil dictionary for the constrained maximising vector.

The later section is explicitly labelled:

> CELL_2 — ZERO-SIDE GUINAND--WEIL CHECK

and independently evaluates the full-space quadratic form after embedding the canonical vector.

### What it established

It connected:

`canonical coefficient vector`
→ `Fourier/test-function representation`
→ `zero constraints`
→ `quadratic form`.

It also exposed the importance of keeping canonical and full coordinates distinct.

### Status

**Established / diagnostic precursor to later dictionary audits.**

---

## Cell 3 — full finite Weil dictionary

### Intended purpose

Cell 3 extended the dictionary work substantially. It reconstructed:

- the complex spectral response;
- `G_complex`;
- the Galerkin quadratic form;
- prime contribution;
- pole contribution;
- Archimedean contribution;
- the Fourier-side `K`/`ghat` machinery.

This was one of the first places where the distinction between the **linear spectral response** and the **quadratic Weil functional** became important.

### Important historical feature

Cell 3 contains its own local implementation of `G_complex`, rather than using the later canonical `cell.py` definition.

### Status

**Diagnostic / superseded.**

---

## Cell 4 — prime functional audit

### Intended purpose

Cell 4 independently audited the prime contribution.

It compares:

1. the direct divided-difference prime matrix;
2. the explicit Guinand–Weil Fourier-side prime expression.

The intended equality is the finite prime-side dictionary.

### What it established

The prime contribution could be independently reconstructed from the prime-power source and checked against the Fourier-side `ghat` representation.

### Status

**Established.**

---

# Cell 5 — the major historical discrepancy

## Cell 5 — Archimedean source / Weil-form comparison

### Intended purpose

Cell 5 was intended to extend the independent dictionary audit to the Archimedean contribution.

It:

- reconstructed the Archimedean source;
- evaluated finite-`T` source values;
- constructed the corresponding source matrix;
- compared the repository Archimedean quadratic form against an explicit Weil-side expression;
- examined the result as `T` increased.

### What went wrong

The historical explicit Archimedean calculation used the quantity

$$\frac1\pi\int_0^T h_+(r) \thinspace \mathop{\mathrm{Re}}G_v(r) \thinspace dr \thinspace$$

where

$$G_v(r)=\sum_k v_kG_k(r).$$

That quantity is a coefficient-weighted **linear** construction, not the quadratic Weil functional represented by the Galerkin matrix.

This produced a substantial and misleading discrepancy and caused a long subsequent investigation.

### Historical significance

This is an important research-history event. The cell is **not to be corrected or rewritten** merely because the error is now understood.

The error motivated:

- the corrected Cell-5 experiments;
- the source/dictionary audits;
- the canonical/full coordinate investigation;
- the eventual `K_v` quadratic construction;
- and, ultimately, the refactoring of `cell.py`.

### Status

**Historical error — deliberately retained.**

---

## Cell 5 corrected variants

### `cell5_corrected.py`

An abandoned attempt to correct the Cell-5 Archimedean calculation by replacing the historical linear `G` expression with the appropriate quadratic construction.

### `cell5_corrected2.py`

A subsequent refinement of that correction.

These files are part of the investigation's history and should not be confused with the original Cell 5.

### Status

**Historical experimental branches.**

---

# Cell 6 — independent Archimedean source/dictionary audit

## Intended purpose

Cell 6 was created to independently reconstruct the Archimedean matrix from

$$\psi_{R,T}(x)=\frac{1}{2\pi^2}\int_{-T}^{T}h_+(r)S(r,x,L) \thinspace dr$$

then form the divided-difference matrix and compare:

1. source-derived Archimedean matrix;
2. repository Archimedean matrix;
3. explicit Weil-side Archimedean quantity.

It also independently reconstructs the prime and pole matrices.

### What it established

This became a much more systematic source-level audit.

It separated:

- source construction;
- divided differences;
- repository matrix decomposition;
- explicit Weil-side calculation.

### Historical discrepancy

The cell originally described the `G_complex` integral as an “explicit Weil-side Archimedean quadratic form”. We now understand that the `G_complex`/`sum_v_G` construction itself is a coefficient-weighted linear sum.

That historical semantic error is **retained**.

The recent API refactor merely changed the executable name from `G_complex` to `sum_v_G`; it did not change the calculation.

### Status

**Diagnostic / discrepancy localisation; historical semantic error retained.**

---

# Cells 6a–6c — follow-up Archimedean investigations

These are progressively narrower investigations following Cell 6.

They explore the Archimedean source, its matrix representation, finite-`T` behaviour, and related discrepancies.

### Status

**Diagnostic / supporting investigations.**

They should be read together with Cells 6 and 7–8 rather than as independent final results.

---

# Cell 7 — symbolic / analytic Archimedean dictionary

### Intended purpose

Cell 7 explicitly states its goal as deriving and numerically verifying

$$
\text{completed-zeta Archimedean factor}
\to h_+(\tau)
\to \text{basis Fourier response}
\to S_x(\tau)
\to \psi_{\rm arch}(x).
$$

It deliberately avoids constructing a large Galerkin matrix.

### What it established

It checked the Archimedean `h_+` convention against the logarithmic derivative of the completed Archimedean factor and exposed possible factor-of-two convention errors.

### Status

**Established / diagnostic.**

---

# Cell 8 — explicit-formula / CvS Archimedean dictionary

### Intended purpose

Cell 8 independently evaluates

$$S(r,x,L)=\int_0^L\sin(2\pi x(1-y/L))\cos(ry) \thinspace dy$$

in two ways:

1. direct quadrature;
2. closed-form trigonometric expression.

It then compares both with the repository kernel and finally compares the resulting source with `operator.py::psi_arch`.

### What it established

It provided an independent check of the Archimedean source kernel and its relationship to the repository operator implementation.

### Status

**Established / important supporting audit.**

---

# Cell 9 — finite-dimensional Weil dictionary audit

### Intended purpose

Cell 9 explicitly tests what a coefficient vector in the trigonometric Galerkin basis represents.

The intended chain is:

$$
v
\to f_v(t)
\to F_v(\tau)
\to \text{translated test function}
\to \text{Weil quadratic form}.
$$

It deliberately uses arbitrary test vectors rather than relying solely on the ground state.

### What it found

It established the direct finite-dimensional Fourier representation but exposed a discrepancy in the spectral/centering convention that required further symbolic analysis.

### Status

**Discrepancy found; resolved by later Cells 10–12.**

---

# Cell 10 — structural spectral dictionary audit

### Intended purpose

Cell 10 resolves the Cell-9 spectral-function discrepancy symbolically.

It determines the exact relationship between:

- the direct Fourier transform;
- the repository `g_k` representation;
- the centering phase;
- the positive/negative frequency conventions.

### What it established

The discrepancy was a **centering/sign convention issue**, not a numerical quadrature failure.

### Status

**Established.**

---

# Cell 11 — exact Fourier / extraction dictionary audit

### Intended purpose

Cell 11 determines exactly which Fourier transform is implemented by `extract_zeros()`.

It distinguishes $F_+(\tau)=\int f(t)e^{+i\tau t} \thinspace dt$
from $F_-(\tau)=\int f(t)e^{-i\tau t} \thinspace dt$

and their centred versions.

It deliberately does not perform zero finding.

### What it established

It pinned down the precise transform/sign/centering convention used by the extraction machinery.

### Status

**Established.**

---

# Cell 12 — closed Weil / spectral quadratic-form audit

### Intended purpose

Cell 12 closes the loop between:

$$
v
\to f_v
\to H_v
\to \text{Weil quadratic form}
$$

and $v^TQv$.

It explicitly builds on the conclusions of Cells 9–11.

### What it established

The spectral representation was shown to reproduce the finite quadratic form encoded by `Q`, subject to the conventions established by the preceding cells.

### Status

**Established / major milestone.**

---

# Cell 13 — controlled Parseval / Fourier-norm audit

### Intended purpose

Cell 13 investigates an apparent Parseval discrepancy.

The numerical frequency integral appeared to approach approximately `0.686` rather than the expected unit norm.

It tests four possibilities:

- quadrature failure;
- finite-frequency tail;
- normalization error;
- incorrect closed transform.

### What it found

The discrepancy was ultimately traced to the representation/normalisation issue rather than a failure of the underlying Fourier dictionary.

### Status

**Discrepancy found; resolved by Cell 14.**

---

# Cell 14 — corrected full/canonical Parseval audit

### Intended purpose

Cell 14 explicitly identifies the Cell-13 problem:

> treating the full `(2N+1)` vector as though it were the canonical `(N+1)` vector.

It checks:

- full → canonical;
- canonical → full;
- round-trip accuracy;
- normalization;
- direct reconstruction of `f(t)`;
- direct Fourier coefficients;
- Parseval.

### What it established

The apparent Parseval discrepancy was a **coordinate-representation error**, not a failure of Parseval or of the Fourier transform.

### Status

**Established.**

---

# Cell 15 — Cell-5 discrepancy closure audit

### Intended purpose

Cell 15 deliberately revisits the historical Cell-5 discrepancy.

It constructs the canonical/full representations correctly and deliberately reproduces the old Cell-5 coordinate mistake.

The central question is whether the historical discrepancy was caused by feeding a full-space vector into a canonical-space functional.

### What it established

It isolated a canonical/full coordinate error that was also present in the historical Cell-5 investigation.

Importantly, this did **not** erase the separate Archimedean `G` category error from Cell 5. The two issues must remain conceptually distinct.

### Status

**Established as a historical-coordinate audit.**

---

# Cell 16 — Archimedean discrepancy localisation

### Intended purpose

Cell 16 localises the remaining Archimedean discrepancy progressively through:
$S
\to dS
\to \text{divided differences}
\to \text{basis kernel}
\to \text{quadratic-form integrand}
\to r\text{-integration}$.

It is explicitly diagnostic rather than production code.

### What it contributed

It narrowed the discrepancy to the Archimedean construction and provided detailed comparisons of the basis-level and vector-level responses.

It also helped make visible the distinction between:

- `G_k`;
- the coefficient-weighted `sum_v_G`;
- genuinely quadratic constructions.

### Status

**Diagnostic / important precursor to Cell 17.**

---

# Cell 17 — quadratic Archimedean kernel audit

### Intended purpose

Cell 17 is the crucial transition to the correct quadratic construction.

It establishes the relationship between:

1. the Archimedean quadratic form obtained from the source/divided differences;
2. the quadratic Volterra kernel $K_v(\omega)=2\int_0^\omega T_v(t)T_v(\omega-t) \thinspace dt$;
3. the direct Archimedean integral constructed from `K_v`.

The key identity under test is $D_v(\omega)=\pi K_v(\omega)$.

### Particularly important design decision

Cell 17 explicitly says:

> This cell deliberately does NOT use `G_complex()`.

That was exactly the correct conceptual separation.

### What it established

It supplied the genuinely quadratic construction that avoids the historical `G`-based category error and provides the route for the subsequent Archimedean dictionary audit.

### Status

**Major established result / direct precursor to the planned Cell 21 work.**

---

# Cell 18 — historical `G_complex` equivalence audit

### Intended purpose

Cell 18 determines whether the historical Cell-5 `G_complex` and the current `cell.py` implementation are mathematically equivalent.

It reconstructs Cell 5's historical `N` and `L` environment rather than guessing them, then compares the two implementations.

### What it established

The historical Cell-5 construction and the current canonical implementation are mathematically equivalent.

The recent refactor changed the current name from `G_complex` to `sum_v_G` without changing its mathematics.

### Important historical distinction

Cell 18 deliberately retains the historical local function name `G_complex` when extracting and inspecting Cell 5. That is correct: the historical source is part of the record.

### Status

**Established.**

---

# Cell 19 — linear-vs-quadratic homogeneity audit

### Intended purpose

Cell 19 is a deliberately extensive numerical audit of the distinction between
coefficient-weighted linear constructions and genuinely quadratic forms.

In particular, it is intended to make the different scaling laws explicit:

$$G(av,r)=aG(v,r)$$

whereas

$$
(av)^*Q(av)=|a|^2v^*Qv.
$$

The purpose is not merely to demonstrate a general mathematical fact, but to
provide a numerical regression/forensic test against the specific objects used
in the Connes–CvS implementation.

### Execution status

Cell 19 is implemented and has been launched, but is computationally very
large. Its expected runtime is measured in days.

Consequently, its final output is not yet available and its conclusions should
not be treated as established until execution completes.

### Historical significance

Cell 19 grew directly out of the `G_complex` investigations. It is intended to
provide an explicit numerical guard against the class of category error that
occurred in Cell 5 and subsequently resurfaced during the Cell 20 investigation.

### Status

**Running — results pending.**

---

# Cell 20 — corrected Archimedean quadratic audit

### Intended purpose

Cell 20 is the major current corrected calculation.

Its stated mathematical correction is to replace the historical Cell-5 expression

$$\frac1\pi\int_0^T h_+(r)\mathop{\mathrm{Re}}G_v(r) \thinspace dr$$

with the required quadratic functional

$$\frac1\pi\int_0^T h_+(r)\int_0^L K_v(1-y/L)\cos(ry) \thinspace dy \thinspace dr.$$

It uses the closed finite-Fourier representation of `K` developed in Cell 17.

### What it is intended to establish

Whether the corrected explicit Archimedean calculation agrees with the repository Archimedean quadratic form.

### Historical treatment

Cell 20 is **not a correction to Cell 5 itself**. It is a new, corrected investigation derived from Cell 5.

The historical Cell 5 remains unchanged.

### Status

**Current substantive audit / important precursor to Cell 21.**

---

# Cell 20a — pole sanity-check forensics

### Intended purpose

Cell 20a investigates an apparent discrepancy between
$\langle u,Q_{\rm pole}u\rangle$ and
$2\mathop{\mathrm{Re}}G_v(i/2)$ where $G_v$ is implemented by `sum_v_G`.

It deliberately performs:

- no Archimedean integration;
- no new Galerkin matrix construction.

Instead it tests the more basic distinction between $P(v)$
as a pole **linear functional**, and $\langle u,Q_{\rm pole}u\rangle$
as a pole **quadratic form**.

It also tests whether $2\mathop{\mathrm{Re}}\sum_vG(v,i/2)$
is proportional to the pole functional.

### What it established

This was a focused forensic investigation into another place where a linear object and a quadratic object could easily be conflated.

The updated version uses `sum_v_G` and explicitly labels its linear character.

### Status

**Diagnostic / confirmed as useful semantic clarification.**

---

# Cell 21 — modern Cell-5 reimplementation

## Intended purpose

Cell 21 was the planned clean reimplementation of what historical Cell 5 was intended to calculate, using the modern `cell.py` machinery and explicitly avoiding the historical `G_complex` / `sum_v_G` category error.

The principal calculation is the genuinely quadratic Archimedean functional

$$
A_{\rm arch}=\frac{1}{\pi}
\int_0^T h_+(r)
\int_0^L
K_v(1-y/L)\cos(ry)\thinspace dy\thinspace dr,
$$

where $K_v$ is the quadratic Volterra kernel established in Cell 17.

The principal path therefore uses

$$
v
\to
K_v
\to
\text{Fourier representation}
\to
A_{\rm arch},
$$

rather than treating

$$\sum_k v_k G_k(r)$$

as a surrogate for the quadratic Weil functional.

The historical Cell-5 expression involving `sum_v_G` is retained in Cell 21 only as a separately labelled forensic comparison.

## What it established

Cell 21 successfully reproduced the corrected Archimedean quadratic calculation and compared it with the Archimedean quadratic form obtained from the Galerkin matrix.

At 20 dps the direct nested numerical calculation produced a small discrepancy relative to the subsequently derived analytic reduction. The Archimedean calculation was therefore repeated at 40 dps.

At 40 dps the independently evaluated nested numerical integral converged to the same value as the analytic reduction to essentially the full available precision:

$$A_{\rm arch}=-1.659033087490935669112988625892145556527\ldots$$

The 20-dps result differed from the high-precision value by approximately $3.13\times10^{-21}$,

whereas the 40-dps result agreed to approximately its full working precision.

This established that the earlier discrepancy was a numerical-precision / nested-quadrature issue rather than a mathematical discrepancy in the quadratic construction.

## Computational significance

The result also exposed the enormous computational cost of evaluating the corrected expression directly.

At 20 dps the Archimedean calculation took approximately 10,874 s (about 3.0 hours).

At 40 dps it took approximately 44,081 s (about 12.2 hours).

The calculation is therefore suitable as an independent validation method, but not as the preferred production implementation.

## Status

**Established — major validation result.**

Cell 21 is now the independent brute-force control against which the analytic Cell-22/23 implementations can be checked.

---

# Cell 22 — analytic Archimedean reduction

## Intended purpose

Cell 22 was created after Cell 21 exposed the prohibitive cost of the nested numerical integration.

The objective was to preserve the same mathematics while analytically evaluating the inner $y$-integral.

The Archimedean calculation was reduced from

$$
\frac{1}{\pi}
\int_0^T
h_+(r)
\int_0^L
K_v(1-y/L)\cos(ry)\thinspace dy\thinspace dr
$$

to

$$
\frac{1}{\pi}
\int_0^T
h_+(r)J_v(r)\thinspace dr
$$

where $J_v(r)$ is evaluated by a finite analytic Fourier sum.

There is therefore only one remaining numerical quadrature: the outer $r$-integral.

## What it established

Cell 22 agreed with the independent Cell-21 result and converged systematically as working precision was increased.

Using the Cell-23 120-dps result as the reference, the approximate absolute discrepancies were:

| dps |            discrepancy |
| --: | ---------------------: |
|  20 | $3.14\times10^{-14}$ |
|  40 | $2.87\times10^{-25}$ |
|  60 | $3.95\times10^{-39}$ |
|  80 | $9.62\times10^{-59}$ |
| 100 | $1.80\times10^{-79}$ |
| 120 |       $\sim10^{-91}$ |

The relatively poor 20-dps result therefore proved to be a precision/conditioning issue. Increasing precision caused the result to converge rapidly to the same value obtained independently by Cell 21 and Cell 23.

## Computational significance

Cell 22 reduced the approximately three-hour Cell-21 calculation at 20 dps to approximately 23 seconds.

At 40 dps it required approximately 58 seconds.

This demonstrated that the expensive nested numerical integration was unnecessary once the inner integral was analytically reduced.

## Status

**Established — validated analytic reduction.**

Cell 22 is retained as the first-generation analytic implementation and as an independent computational route to the Cell-23 result.

---

# Cell 23 — optimised analytic Archimedean calculation

## Intended purpose

Cell 23 takes the analytic reduction established in Cell 22 and improves its numerical and computational efficiency without changing the mathematics.

The principal optimisations are:

* exploit the symmetry between the $(m,n)$ and $(n,m)$ terms;
* evaluate each Fourier-mode $S_m(r)$ only once for a given $r$;
* use analytically stable forms for expressions such as $1-\cos(x)$;
* express the kernel integral using a sinc-based representation with the removable $k=0$ limit handled explicitly.

The resulting calculation still evaluates

$$A_{\rm arch}=\frac{1}{\pi}\int_0^T h_+(r)J_v(r)\thinspace dr$$

but does substantially less repeated arithmetic than Cell 22.

## What it established

Cell 23 converges extremely rapidly with working precision.

Using the 120-dps Cell-23 result as the reference, the approximate discrepancies are:

| dps | discrepancy |
| --: | ----------: |
| 20 | $1.44\times10^{-21}$ |
| 40 | $2.15\times10^{-41}$ |
| 60 | $4.16\times10^{-62}$ |
| 80 | $2.01\times10^{-82}$ |
| 120 | reference |

The 60-, 80-, 100- and 120-dps calculations demonstrate very strong stability of the resulting value.

The limiting value is

$$A_{\rm arch}=-1.6590330874909356691129886258921455565271140176152095515115701580412723893977269298\ldots$$

The agreement with Cell 21 at 40 dps is particularly important because Cell 21 obtains the result through the original nested numerical integration rather than through the analytic reduction.

Cell 22 independently converges to the same value, although more slowly with respect to working precision.

Thus Cell 23 is supported by two independent computational routes:

$$\text{Cell 21}\longrightarrow A_{\rm arch}$$

and

$$\text{Cell 22}\longrightarrow A_{\rm arch}$$

with both converging to the Cell-23 result.

## Computational significance

The analytic optimisation is dramatic.

At 40 dps, Cell 23 requires only a few seconds compared with approximately 12 hours for Cell 21.

At 120 dps, Cell 23 still completes in roughly a minute.

The brute-force nested integral is therefore no longer an appropriate production method for this calculation. It is best regarded as an independent validation method.

## Status

**Established — current preferred Archimedean computational implementation.**

Cell 23 is now the natural basis for future high-precision Archimedean calculations.

---

# Cells 19 and 5_corrected — implications of the analytic reduction

## Cell 19

Cell 19 was designed as an extensive numerical audit of the distinction between coefficient-weighted linear constructions and genuinely quadratic forms.

Its mathematical purpose remains legitimate, but its implementation performs very expensive nested numerical integrations.

At the current recorded run, Cell 19 at 50 dps had already consumed approximately

$$6325\ {\rm minutes}\approx105.4\ {\rm hours},$$

with no useful final result yet available.

This is now understood to be a computationally obsolete route for the quadratic Archimedean calculations.

The completion of Cell 21, together with the successful Cell 22/23 reductions, means that the mathematical question motivating Cell 19 has already been addressed through substantially cheaper and independently validated machinery.

The running calculation may be retained temporarily as a historical experiment, but it should not be regarded as a prerequisite for the current research programme.

### Status

**Running / computationally superseded.**

---

## `cell5_corrected.py`

`cell5_corrected.py` is another historical branch that evaluates the corrected Archimedean quantity using expensive nested numerical integrations, including repeated calculations at different values of $T$.

At the current recorded run, the 80-dps calculation had consumed approximately

$$8551\ {\rm minutes}\approx142.5\ {\rm hours}.$$

It had completed the $T=20$ calculation and was still working toward the next $T$ value.

The successful Cell-21 result has now provided an independent brute-force validation of the corrected quadratic construction, while Cells 22 and 23 have demonstrated that the same mathematical quantity can be evaluated vastly more efficiently.

Consequently, there is no longer a mathematical need to complete the entire historical `cell5_corrected.py` sequence merely to establish the correctness of the quadratic Archimedean construction.

The historical file should remain untouched. Any decision to terminate its current long-running calculation is a computational-resource decision, not a change to the historical record.

### Status

**Historical experimental branch / computationally superseded.**

---

# Current research state after Cells 21–23

The Archimedean investigation has now passed an important transition.

The historical sequence

$$
\text{Cell 5}
\to
\text{discrepancy}
\to
\text{Cells 6–20}
$$

has established why the original calculation failed and has produced the correct quadratic construction.

Cell 21 then provided an independent brute-force implementation of that construction.

Cells 22 and 23 subsequently transformed the calculation from an impractical nested numerical integral into an efficient analytic computation.

The current picture is therefore:

> $$
\begin{array}{c}
\text{Cell 21}\\
\text{independent nested numerical validation}
\end{array}
$$

$$\Downarrow$$

> $$
\begin{array}{c}
\text{Cell 22}\\
\text{analytic reduction} \Downarrow
\end{array}
$$

$$\Downarrow$$

> $$
\begin{array}{c}
\text{Cell 23}\\
\text{optimised analytic implementation}
\end{array}
$$

The three routes converge to the same Archimedean value.

This substantially increases confidence in the quadratic Archimedean dictionary and removes the need to use the expensive nested integral as the normal computational route.

# Efficient calculation of the Archimedean integral

## Derivation of the analytic Fourier representation `K_fourier`

The purpose of `K_fourier` is to evaluate the Fourier-side representation of the quadratic kernel

$$
K_v(\omega) = 2\int_0^\omega T_v(t)\thinspace T_v(\omega-t)\thinspace dt,
\qquad 0\leq\omega\leq1
$$

without performing the inner convolution integral numerically.

The derivation below starts from this convolution definition and reduces it analytically to a finite sum in the canonical coefficient vector

$$
v=(v_0,\ldots,v_N).
$$

### 1. Full symmetric Fourier representation

For the derivation, temporarily introduce the corresponding full symmetric coefficients $u_m$:

$$
u_0=v_0,
\qquad
u_{+m}=u_{-m}=\frac{v_m}{\sqrt{2}},
\qquad
1\leq m\leq N.
$$

Thus

$$
T_v(t) = \sum_{m=-N}^{N} u_m e^{2\pi i m t}.
$$

The vector $u$ is used only as an intermediate mathematical representation in this derivation. The computational interface of `cell.py` is canonical $v$.

Define

$$
a_m=\frac{2\pi m}{L}.
$$

The Fourier-side quantity evaluated by `K_fourier` is

$$
J_v(r) = \int_0^L K_v\left(1-\frac{y}{L}\right) \cos(ry)\thinspace dy.
$$

The substitution

$$
\omega=1-\frac{y}{L}
$$

is useful because the integer Fourier modes satisfy

$$
e^{2\pi i m\omega} = e^{2\pi i m(1-y/L)} = e^{-ia_m y}.
$$

### 2. Expand the convolution

Substitute the Fourier expansion of $T_v$ into the convolution:

$$
K_v(\omega) = 2 \sum_{m,n=-N}^{N} u_m u_n e^{2\pi i n\omega}
\int_0^\omega e^{2\pi i(m-n)t}\thinspace dt.
$$

For $m=n$,

$$
2e^{2\pi i m\omega} \int_0^\omega 1\thinspace dt = 2\omega e^{2\pi i m\omega}.
$$

For $m\ne n$,

$$
2e^{2\pi i n\omega} \int_0^\omega e^{2\pi i(m-n)t}\thinspace dt =
\frac{1}{\pi i(m-n)} \left( e^{2\pi i m\omega} - e^{2\pi i n\omega} \right).
$$

Therefore

$$
K_v(\omega) = 2\sum_{m=-N}^{N} u_m^2\thinspace\omega e^{2\pi i m\omega} +
\sum_{m,n=-N}^{N} \frac{u_m u_n}{\pi i(m-n)}
$$

$$
K_v(\omega) = 2\sum_{m=-N}^{N} u_m^2\thinspace\omega e^{2\pi i m\omega} +
\sum_{m,n=-N}^{N} \frac{u_m u_n}{\pi i(m-n)}
\left( e^{2\pi i m\omega} - e^{2\pi i n\omega} \right).
$$

K_v(\omega) = 2\sum_{m=-N}^{N} u_m^2\thinspace\omega e^{2\pi i m\omega}
+ \sum_{\substack{m,n=-N \\ m \ne n}}^{N} \frac{u_m u_n}{\pi i(m-n)}
\left( e^{2\pi i m\omega} - e^{2\pi i n\omega} \right).
$$

### 3. Fourier transform of the diagonal terms

Set

$$
\omega=1-\frac{y}{L}.
$$

For the diagonal terms,

$$
\omega e^{2\pi i m\omega} = \left(1-\frac{y}{L}\right)e^{-ia_m y}.
$$

Taking the real part after multiplication by $\cos(ry)$ gives

$$
\int_0^L \left(1-\frac{y}{L}\right) \cos(a_m y)\cos(ry)\thinspace dy.
$$

Define

$$
C_m(r) = \int_0^L \left(1-\frac{y}{L}\right) \cos(a_m y)\cos(ry)\thinspace dy.
$$

Hence the complete diagonal contribution is

$$
J_{\mathrm{diag}}(r) = 2\sum_{m=-N}^{N} u_m^2 C_m(r).
$$

Using the canonical symmetry,

$$
u_0=v_0, \qquad u_{\pm m}^2=\frac{v_m^2}{2},
$$

and $C_{-m}=C_m$, this becomes

> $$
J_{\mathrm{diag}}(r) = 2\sum_{m=0}^{N} v_m^2 C_m(r).
$$

This is the `diag` term in `K_fourier`.

### 4. Fourier transform of an off-diagonal ordered pair

For $m\ne n$, after setting

$$
\omega=1-\frac{y}{L},
$$

the exponential difference becomes

$$
e^{-ia_m y}-e^{-ia_n y}.
$$

Define

$$
S_m(r) = \int_0^L \sin(a_m y)\cos(ry)\thinspace dy.
$$

Since

$$
e^{-ia_m y} = \cos(a_m y)-i\sin(a_m y),
$$

the real part of the Fourier transform of the off-diagonal term is

$$
\frac{u_m u_n}{\pi} \frac{S_n(r)-S_m(r)}{m-n}.
$$

Thus the full off-diagonal contribution is

$$
J_{\mathrm{off}}(r) = \sum_{\substack{m,n=-N\\m\ne n}}^N
\frac{u_m u_n}{\pi} \frac{S_n(r)-S_m(r)}{m-n}.
$$

The summand is symmetric under interchange of $m$ and $n$:

$$
u_m u_n \frac{S_n-S_m}{m-n} = u_n u_m \frac{S_m-S_n}{n-m}.
$$

Therefore the two ordered terms $(m,n)$ and $(n,m)$ are identical, and

> $$
J_{\mathrm{off}}(r) = \frac{2}{\pi}\sum_{m<n} u_m u_n\frac{S_n(r)-S_m(r)}{m-n}
.
$$

This is the important multiplicity point: restricting the ordered sum to $m<n$ introduces exactly the displayed factor of $2$. There is no additional factor of $2$ merely because the sum is now triangular.

### 5. Reduce the off-diagonal sum to canonical coordinates

We now substitute

$$
u_0=v_0, \qquad u_{\pm m}=\frac{v_m}{\sqrt2}.
$$

The off-diagonal terms naturally split into four classes.

#### 5.1 Terms involving the zero mode

For each $m>0$, the pairs $(0,m)$ and $(-m,0)$ give equal contributions.

Since $S_0=0$,

$$
J_{0,m}+J_{-m,0} = -\frac{2\sqrt2}{\pi} v_0v_m \frac{S_m(r)}{m}.
$$

Summing over $m$ gives

> $$
J_{0}(r) = -\frac{2\sqrt2\thinspace v_0}{\pi}\sum_{m=1}^{N}\frac{v_m S_m(r)}{m}.
$$

This is the `off_zero` term in `K_fourier`.

#### 5.2 The opposite-sign pair $(m,-m)$

For each $m>0$, the pair $(-m,m)$ gives

$$
\frac{2}{\pi} \frac{v_m^2}{2} \frac{S_m-S_{-m}}{-2m}.
$$

Because

$$
S_{-m}=-S_m,
$$

this reduces to

$$
-\frac{v_m^2}{\pi} \frac{S_m(r)}{m}.
$$

Therefore

> $$
J_{\pm m}(r) = -\frac{1}{\pi} \sum_{m=1}^{N} \frac{v_m^2 S_m(r)}{m}.
$$

This is the `off_diag` term in `K_fourier`.

#### 5.3 Same-sign positive and negative pairs

For $1\le m<n\le N$, the positive pair $(m,n)$ contributes

$$
\frac{v_m v_n}{\pi} \frac{S_n-S_m}{m-n}.
$$

The negative pair $(-n,-m)$ contributes the same amount.

Together,

$$
\frac{2v_m v_n}{\pi} \frac{S_n-S_m}{m-n}.
$$

#### 5.4 Mixed-sign pairs

The two mixed-sign pairs $(-n,m)$ and $(-m,n)$ together contribute

$$
-\frac{2v_m v_n}{\pi} \frac{S_m+S_n}{m+n}.
$$

Combining the same-sign and mixed-sign contributions gives

$$
\frac{2v_m v_n}{\pi} \left[ \frac{S_n-S_m}{m-n} - \frac{S_m+S_n}{m+n} \right].
$$

The bracket simplifies algebraically:

$$
\begin{aligned}
\frac{S_n-S_m}{m-n} - \frac{S_m+S_n}{m+n}
&=
\frac{(S_n-S_m)(m+n)-(S_m+S_n)(m-n)} {(m-n)(m+n)} \\[4pt]
&= \frac{2(mS_m-nS_n)} {n^2-m^2}.
\end{aligned}
$$

Hence the complete positive-mode pair contribution is

> $$
J_{m,n}(r) = \frac{4v_m v_n}{\pi} \frac{mS_m(r)-nS_n(r)} {n^2-m^2},
\qquad 1\le m<n\le N.
$$

This is the $O(N^2)$ triangular term evaluated by `K_fourier`.

### 6. Final canonical formula

Combining the diagonal, zero-mode, opposite-sign, and positive-mode pair contributions gives

$$
\begin{aligned}
J_v(r) ={}&
2\sum_{m=0}^{N} v_m^2 C_m(r)
\\ &- \frac{2\sqrt2\thinspace v_0}{\pi} \sum_{m=1}^{N} \frac{v_mS_m(r)}{m}
\\ &- \frac{1}{\pi} \sum_{m=1}^{N}
\frac{v_m^2S_m(r)}{m} \\
&+ \frac{4}{\pi}
\sum_{1\le m<n\le N} v_m v_n
\frac{mS_m(r)-nS_n(r)} {n^2-m^2}.
\end{aligned}
$$

This is exactly the decomposition implemented by `K_fourier`.

The calculation therefore performs no approximation to the quadratic kernel. It is an algebraic reduction of the original convolution integral to a finite Fourier sum. The numerical gain comes from evaluating the Fourier transform analytically and computing each $S_m(r)$ only once for a given $(r,L)$.

### 7. Relationship to the computational implementation

The implementation corresponds term-by-term to the derivation:

* `diag` evaluates

  $$2\sum_{m=0}^{N}v_m^2C_m(r);$$

* `off_diag` evaluates

  $$-\frac{1}{\pi}\sum_{m=1}^{N}\frac{v_m^2S_m(r)}{m};$$

* `off_zero` evaluates

  $$-\frac{2\sqrt2\thinspace v_0}{\pi}\sum_{m=1}^{N}\frac{v_mS_m(r)}{m};$$

* the triangular loop evaluates

  $$\frac{4}{\pi}\sum_{1\le m<n\le N}v_m v_n\frac{mS_m(r)-nS_n(r)}{n^2-m^2}.$$

Thus `K_fourier` is mathematically equivalent to the original convolution construction `K_canonical`, while avoiding the numerical inner integration over the convolution variable.


---

# Future directions

The purpose of this section is deliberately forward-looking. As the investigation progresses, these directions will themselves become part of the historical record.

## 1. Adopt Cell 23 as the production Archimedean implementation

Future calculations requiring the Archimedean quadratic functional should preferentially use the Cell-23 analytic machinery.

The nested Cell-21 calculation should be reserved for occasional independent validation.

The historical Cell-5 and Cell-19 implementations should remain available for forensic purposes but should not be used as normal computational infrastructure.

## 2. Add CPU timing to future calculations

The previous cells generally recorded wall-clock time. This is useful for knowing how long a calculation took on the machine, but it does not distinguish computation from waiting for CPU availability.

Future substantial cells should preferably record both:

* wall-clock time;
* process CPU time.

This will allow computational cost to be compared meaningfully when multiple long-running experiments share the same machine.

Historical cells should not be rewritten merely to add this information.

## 3. Revisit finite-$T$ behaviour using the efficient implementation

The analytical Cell-23 route makes it practical to investigate the dependence on $T$ at substantially higher precision and with many more samples than was practical using the nested integral.

This should allow the finite-$T$ convergence of the Archimedean contribution to be studied without spending days on each calculation.

## 4. Consolidate the validated Archimedean dictionary

The work of Cells 17, 20, 21, 22 and 23 should eventually be distilled into a clean statement of the finite-dimensional Archimedean dictionary:

$$
v
\to
T_v
\to
K_v
\to
\widehat g_v
\to
Q_{\rm arch}.
$$

The goal is not merely a fast numerical routine, but a transparent mathematical chain in which the semantic type of every object is clear.

## 5. Resume the broader Connes–CvS investigation

The project should now be able to move beyond the historical Cell-5 discrepancy.

The central computational infrastructure is considerably better understood:

* the finite Fourier dictionary has been audited;
* canonical/full coordinates are understood;
* the prime contribution has been independently checked;
* the Archimedean source has been independently checked;
* the genuinely quadratic $K_v$ construction has been established;
* the historical `G_complex` / `sum_v_G` distinction is understood;
* the Archimedean quadratic calculation now has an efficient, independently validated implementation.

The next work should therefore focus on what mathematical consequences can be extracted from the validated finite-dimensional construction, rather than continuing to reproduce the historical numerical calculations.

---

# Updated major historical arc

```text
Cells 0–4
    Initial reconstruction
    ↓
Cell 5
    Archimedean discrepancy discovered
    ↓
Cells 5_corrected*
    Attempts to repair / understand discrepancy
    ↓
Cell 6
    Independent Archimedean source/matrix audit
    ↓
Cells 7–8
    Archimedean source dictionary
    ↓
Cells 9–12
    Fourier / centering / Weil dictionary
    ↓
Cells 13–15
    Parseval and canonical/full-coordinate discrepancy
    ↓
Cell 16
    Localise remaining Archimedean discrepancy
    ↓
Cell 17
    Correct quadratic K_v construction established
    ↓
Cell 18
    Historical G_complex equivalence confirmed
    ↓
Cell 19
    Extensive linear-vs-quadratic numerical audit
    [running; computationally superseded]
    ↓
Cell 20
    Corrected Archimedean quadratic audit
    ↓
Cell 20a
    Pole linear-vs-quadratic forensic audit
    ↓
Cell 21
    Clean modern Cell-5 reimplementation
    Independent brute-force validation
    ↓
Cell 22
    Analytic elimination of inner y-integral
    ↓
Cell 23
    Optimised analytic Archimedean implementation
    ↓
Current direction
    Efficient high-precision Archimedean calculations
    + broader Connes–CvS mathematical investigation
```

## Current status summary

At the current stage:

* The finite Fourier/zero-side dictionary is substantially audited.
* The canonical/full coordinate distinction is understood.
* The prime-side dictionary has been independently audited.
* The Archimedean source itself has been independently audited.
* The genuinely quadratic $K_v$ construction has been established and tested.
* The historical `G_complex` construction has been confirmed mathematically equivalent to the current `sum_v_G`.
* The historical Cell-5/Cell-6 misuse of that linear construction remains preserved as historical record.
* Cell 21 has independently validated the corrected quadratic Archimedean calculation through expensive nested numerical integration.
* Cell 22 has analytically eliminated the inner numerical integral.
* Cell 23 has optimised that analytic calculation and demonstrated stable high-precision convergence.
* The three computational routes converge to the same Archimedean value.
* The expensive nested-integral implementations are now best regarded as validation/forensic calculations rather than production machinery.
* The investigation is ready to move beyond the Cell-5 discrepancy and back toward the broader mathematical objectives of the project.

