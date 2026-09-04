# Connes–CvS exploratory cell history

**Repository:** `nrensen/connes-cvs-`  
**Historical snapshot audited:** commit `150fa5fe3788018d7582d67d488d3c95a314a155`  
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

# Cells 19 and 5_corrected — retrospective implications of the analytic reduction

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

# Cell 24 — finite-T Archimedean convergence map

### Intended purpose

Cell 24 begins the post-validation investigation of $\text{finite-}T$ behaviour using the now-validated analytic Archimedean implementation.

The forensic ground state is held fixed:

* $c=13$
* $N=8$
* Galerkin $T=400$
* generation precision 150 dps.

Only the upper limit of the Archimedean $r$-integral is varied:

$$
A_{\rm arch}(T) = \frac1\pi\int_0^T h_+(r)J_v(r)\thinspace dr.
$$

No Galerkin matrix is rebuilt and no numerical $y$-integration is performed.

### What it established

It provided a systematic $\text{finite-}T$ convergence map using the efficient analytic representation and confirmed that the fixed forensic state could be reused while extending the Archimedean cutoff.

### Status

**Established / transition to long-range tail investigation.**

---

# Cell 25 — finite-T Archimedean convergence and historical cross-check

### Intended purpose

Cell 25 uses the analytic `K_fourier` reduction to reproduce key $\text{finite-}T$ values from the historical `cell5_corrected` calculation and then extends the cutoff substantially beyond the historical $T=200$ range.

The ground state remains fixed at

$$
c=13,\qquad N=8,\qquad T_{\rm ground}=60,
$$

with generation precision 80 dps.

### What it established

The analytic calculation reproduced the historical $\text{finite-}T$ Archimedean values while making much larger cutoff experiments computationally practical.

This established an important continuity between the historical source-level calculation and the modern analytic implementation.

### Status

**Established — historical cross-check and extended $\text{finite-}T$ investigation.**

---

# Cell 26 — long-range forensic Archimedean tail

### Intended purpose

Cell 26 extends the $\text{finite-}T$ investigation to $T=10,000$ for the fixed forensic ground state

$$
c=13,\qquad N=8,\qquad T_{\rm ground}=400,
$$

using the analytic Archimedean integrand

$$
I(r)=h_+(r)K_{\rm fourier}(v_\star,r,L).
$$

No ground-state regeneration occurs as $T$ changes.

### What it found

The $\text{finite-}T$ Archimedean contribution continued to exhibit a small but persistent $\text{large-}T$ contribution, motivating a direct investigation of the pointwise tail rather than immediately assuming a particular asymptotic law.

### Historical significance

Cell 26 is the point at which the investigation changed from ordinary finite-cutoff validation to an explicit study of the asymptotic Archimedean tail.

### Status

**Diagnostic / precursor to Cells 27–32.**

---

# Cell 27 — Archimedean tail anatomy

### Intended purpose

Cell 27 examines the $\text{large-}r$ pointwise structure of

$$
J(r)=K_{\rm fourier}(v_\star,r,L)
$$

and

$$
I(r)=h_+(r)J(r)
$$

for the fixed forensic ground state.

It deliberately performs no long-range integration and makes no asymptotic assumption.

It also introduces phase-locked samples satisfying

$$
rL=k\pi
$$

and

$$
rL=(k+\tfrac12)\pi
$$

in order to distinguish ordinary decay from oscillatory structure.

### What it contributed

It established that the $\text{large-}r$ integrand is strongly oscillatory and that ordinary pointwise ratios are not a reliable way to infer its asymptotic decay.

The phase-locked samples showed that the oscillation is strongly tied to the phase $rL$.

### Status

**Diagnostic / structural precursor to the integrated-tail investigation.**

---

# Cell 28 — direct integrated Archimedean tail

### Intended purpose

Cell 28 abandons pointwise extrapolation and integrates the actual analytic integrand over successive finite intervals:

$$
A(a,b)=\int_a^b h_+(r)K_{\rm fourier}(v_\star,r,L)\thinspace dr.
$$

No asymptotic power law is assumed.

### What it contributed

It showed that the successive signed interval contributions remained positive over the investigated range and provided the first direct evidence that the tail was not simply disappearing through obvious local cancellation.

It also supplied interval-level data suitable for subsequent logarithmic/dyadic analysis.

### Status

**Diagnostic / precursor to Cells 29–30.**

---

# Cell 29 — log-scale integrated-tail scaling

### Intended purpose

Cell 29 introduces the dyadic interval quantity

$$
D(T)=\int_T^{2T}I(r)\thinspace dr
$$

and examines the empirical ratio

$$
\frac{D(2T)}{D(T)}
$$

without assuming a value for the power-law exponent.

It also tracks the cumulative tail and diagnostic quantities such as $T^pD(T)$.

### What it contributed

The dyadic contributions remained positive and decreased with $T$. The effective scaling appeared broadly compatible with a decay somewhat slower than $1/T$ over the explored range.

However, the data did not establish a particular asymptotic law.

### Status

**Diagnostic / exploratory asymptotic analysis.**

---

# Cell 30 — forensic asymptotic tail test

### Intended purpose

Cell 30 extends the dyadic analysis from $T\leq20,480$ to

$$
T=20,971,520
$$

and tests the empirical hypothesis

$$
D(T)\sim\frac{C}{T}.
$$

It reports

$$
C_T=TD(T)
$$

and an empirical local exponent $p_{\rm eff}$, together with several extrapolation diagnostics.

### What it appeared to show

The computed $D(T)$ values remained positive over the extended range and continued to decrease approximately on the scale of $1/T$, while $C_T$ continued to drift upward rather than reaching an obvious plateau.

This suggested that a logarithmic correction such as

$$
D(T)\sim\frac{\log T}{T}
$$

might be compatible with the observed trend.

### Critical qualification

The Cell-30 integrals were evaluated using unsubdivided `mp.quad` over extremely large oscillatory intervals. Cell 31 subsequently demonstrated that changing interval subdivision changes these values by percent-level amounts even when working precision is increased substantially.

Consequently, the numerical values and asymptotic interpretation of Cell 30 are **not established**.

Cell 30 should therefore be preserved as an important historical hypothesis-generating experiment, not as evidence for a $1/T$ or $(\log T)/T$ asymptotic law.

### Status

**Diagnostic / superseded as a quantitative tail estimate by Cell 31.**

---

# Cell 31 — large-T quadrature forensic

### Intended purpose

Cell 31 was created to determine whether the $\text{large-}T$ results of Cell 30 were actually resolving the highly oscillatory integrand.

It independently varies:

1. working precision;
2. interval subdivision.

The mathematical integrand and fixed forensic ground state are unchanged.

### What it established

The precision sweep showed essentially identical results at 80, 100 and 120 dps for the large test intervals. Increasing numerical precision therefore did **not** resolve the discrepancy.

In contrast, subdividing the same intervals changed the computed integrals by percent-level amounts.

This establishes that the principal numerical problem in Cell 30 was **quadrature resolution of the oscillatory interval**, rather than insufficient arithmetic precision.

The result is a critical methodological warning:

> Agreement across working precision is not sufficient evidence of convergence when an oscillatory integral is being evaluated over an enormous interval.

### Consequence for the preceding tail analysis

The numerical values of $D(T)$, $C_T$, and $p_{\rm eff}$ reported by Cell 30 cannot presently be treated as quantitatively converged.

In particular, the apparent $1/T$ behaviour, upward drift of $TD(T)$, and possible logarithmic correction remain unresolved.

The next step should therefore be to understand the analytic frequency structure of $K_{\rm fourier}$ and construct a quadrature method adapted to that structure, rather than simply increasing precision or blindly increasing the number of interval subdivisions.

### Status

**Established — critical numerical-methodology result.**

Cell 31 supersedes the quantitative conclusions of Cell 30 while preserving Cell 30's role as a hypothesis-generating experiment.

---

# Cell 32 — exact analytical structure of the Archimedean tail

### Intended purpose

Cell 32 responds directly to the quadrature problem exposed by Cell 31.

Rather than attempting further $\text{large-}T$ integration, it analytically reduces the exact $\text{finite-}N$ `K_fourier` expression using

$$
a_m=\frac{2\pi m}{L},
\qquad
a_mL=2\pi m.
$$

This reveals the common oscillatory structure of the Fourier modes.

### What it established

For the $\text{finite-}N$ kernel,

$$
K_{\rm fourier}(v,r,L) = (1-\cos rL)R_v(r),
$$

where $R_v(r)$ is a purely rational function of $r$.

The factorisation was independently checked numerically against the existing `K_fourier` implementation at high precision.

The $\text{large-}r$ behaviour was then found to be

$$
R_v(r) = \frac{A(v)}{r^2} + O(r^{-4}),
$$

with

$$
A(v) = \frac{2}{L} \left( v_0+\sqrt2\sum_{m=1}^{N}v_m \right)^2.
$$

Equivalently, writing

$$
T_v(0) = v_0+\sqrt2\sum_{m=1}^{N}v_m,
$$

the leading coefficient is

$$
A(v)=\frac{2T_v(0)^2}{L}.
$$

For the $c=13,N=8$ forensic ground state, $T_v(0)$ is extraordinarily small, so the nominal $r^{-2}$ tail is strongly suppressed.

Cell 32 also numerically probes the next $r^{-4}$ coefficient, providing a target for subsequent exact symbolic derivation.

### Mathematical significance

This is the first point in the tail investigation at which the $\text{large-}r$ structure is explained analytically rather than inferred from numerical integration.

It also changes the numerical problem fundamentally: the tail is no longer an opaque highly oscillatory function. Its dominant oscillatory factor and leading rational decay are explicitly known.

### Status

**Major established analytical result / foundation for subsequent tail bounds and N-dependence analysis.**

---

# Cell 33 — N-dependence of the leading Archimedean tail coefficient

### Intended purpose

Cell 33 begins the complementary $N$-dependence investigation suggested by Cell 32.

For each finite Galerkin dimension $N$, Cell 32 gives the exact leading coefficient

$$
A_N=\frac{2}{L}T_{v_N}(0)^2,
$$

where $v_N$ is the corresponding ground-state vector and

$$
T_v(0)=v_0+\sqrt2\sum_{m=1}^{N}v_m.
$$

Cell 33 therefore asks whether the extraordinary suppression of $T_v(0)$ observed for the forensic $N=8$ state is an isolated numerical feature or a systematic property of the Galerkin ground state.

The cell:

* constructs ground states for a sequence of even $N$;
* evaluates $T_v(0)$;
* computes the exact leading tail coefficient $A_N$;
* numerically extracts the next coefficient $B$ from the exact reduced rational kernel;
* checks the stability of the extracted $B$ as $r$ increases;
* reports $\|v\|$ and $\lambda_{\min}$ as numerical sanity checks.

No Archimedean tail integration is performed.

### What it established

The survey showed that the strong suppression of $T_v(0)$, and hence of the leading coefficient $A_N$, is not peculiar to a single $N$. The leading tail coefficient decreases extremely rapidly as the Galerkin dimension increases.

The experiment also supplied stable numerical estimates of the next $r^{-4}$ coefficient, providing a concrete target for an exact asymptotic derivation.

At this stage, however, the $B$ coefficient remained numerical and no $N\to\infty$ law was claimed.

### Status

Diagnostic / structural $N$-dependence survey; precursor to the exact coefficient analysis of Cells 36–38.

---

# Cell 34 — systematic $N$-scan of ground-state tail structure

### Intended purpose

Cell 34 extends Cell 33 from a sparse $\text{even-}N$ survey to every integer

$$
N=1,\ldots,24,
$$

deliberately including both odd and even dimensions.

The primary question remains whether the rapid suppression of

$$
T_v(0)=v_0+\sqrt2\sum_{m=1}^{N}v_m
$$

is systematic as $N$ increases.

The cell also records:

* $T_v(L)$;
* the low-order spectral moments $M_2$ and $M_4$;
* the exact leading coefficient

$$A_N=\frac{2T_v(0)^2}{L};$$

* numerical extraction of the next coefficient $B$;
* $\lambda_{\min}$;
* $\text{successive-}N$ scaling ratios.

No Archimedean integration is performed. Existing ground-state cache entries are reused where available.

### What it contributed

The broader $N$-scan confirmed that the suppression of the endpoint quantity $T_v(0)$ persists across successive Galerkin dimensions and is not an artefact of selecting only even $N$.

The scan also showed that the suppression extends beyond the leading coefficient and is accompanied by strong structure in the low-order spectral moments. This motivated the more targeted endpoint-jet investigation in Cell 35.

The numerical $B$ extraction remained exploratory: it was useful for identifying the next coefficient but was not yet an analytical result.

### Status

Diagnostic / structural $N$-dependence survey; superseded quantitatively by the exact $\text{finite-}N$ coefficient derivation in Cells 36–38.

---

# Cells 35–38 — endpoint jets and the exact Archimedean tail hierarchy

*Updated 4 September 2026.*

This sequence marks the transition from empirical $\text{finite-}N$ tail observations to an exact analytical description of the complete inverse-power tail of the $\text{finite-}N$ Archimedean kernel.

The starting point is the exact factorisation established in Cell 32,

$$
K_{\mathrm{Fourier}}(v,r,L) = (1-\cos(rL))R_v(r),
$$

with

$$
R_v(r)\sim \frac{A_0}{r^2} + \frac{A_1}{r^4} + \frac{A_2}{r^6} + \cdots.
$$

Cell 32 established

$$
A_0=\frac{2}{L}T_v(0)^2.
$$

Cells 35–38 show that this is the first member of an exact hierarchy governed by the even endpoint jet of the finite-band test function.

## Cell 35 — endpoint jets and spectral moments

### Intended purpose

Cell 35 investigates the endpoint data of

$$
T_v(t) = v_0+\sqrt2\sum_{m=1}^{N} v_m\cos\left(\frac{2\pi mt}{L}\right).
$$

For $k\ge1$,

$$
T_v^{(2k)}(0) = \sqrt2\thinspace(-1)^k \left(\frac{2\pi}{L}\right)^{2k} M_{2k},
$$

where

$$
M_{2k} = \sum_{m=1}^{N}m^{2k}v_m.
$$

The cell surveys $T_v(0)$, $T_v(L)$, and the even endpoint derivatives through $T_v^{(8)}(0)$, together with the corresponding dimensionless spectral moments.

### What it established

The ground-state vectors exhibit strong suppression of endpoint quantities as $N$ increases. In particular, the small value of $T_v(0)$ identified in Cell 32 is part of a broader pattern involving higher even endpoint derivatives.

This suggested that the increasingly small Archimedean tail may be connected to increasing suppression of the endpoint jet.

At this stage the connection remained a structural observation. No asymptotic coefficient was yet expressed analytically in terms of the endpoint data.

### Status

Established structural observation / precursor to Cells 36–38.

---

## Cell 36 — exact finite $N$ tail coefficients

### Intended purpose

Cell 36 converts the $\text{large-}r$ expansion of the exact reduced rational kernel into an explicit algebraic calculation.

For

$$
\kappa=\frac{2\pi}{L},
$$

and

$$
H_k(m,n) = \frac{n^{2k+2}-m^{2k+2}}{n^2-m^2} = \sum_{j=0}^{k}n^{2(k-j)}m^{2j},
$$

the coefficient for $k\ge1$ is

$$
A_k = \frac{4(k+1)}{L}\kappa^{2k} \sum_m m^{2k}v_m^2 +
\frac{4\sqrt{2}}{L}\kappa^{2k}v_0M_{2k} +
\frac{4}{\pi}\kappa^{2k+1} \sum_{m < n}v_mv_nH_k(m,n).
$$

### What it established

The coefficients obtained from this expression agree with numerical extraction from the exact rational kernel. Successive subtraction of the asymptotic terms exposes the predicted next inverse power of $r$.

Thus the $\text{finite-}N$ inverse-power hierarchy is established algebraically rather than inferred from numerical fitting.

The result also identifies the precise polynomial structure that must be reorganised to obtain an endpoint formulation.

### Status

Established — exact $\text{finite-}N$ asymptotic coefficient formula.

---

## Cell 37 — moment-convolution identity

### Intended purpose

Cell 37 examines the polynomial $H_k(m,n)$ appearing in Cell 36 and asks whether the pairwise spectral interaction can be expressed entirely through ordinary spectral moments.

Using the diagonal continuation

$$
H_k(m,m)=(k+1)m^{2k},
$$

define

$$
Q_k = \sum_{m,n\ge1} v_mv_nH_k(m,n).
$$

The polynomial identity gives

$$
Q_k = \sum_{j=0}^{k} M_{2j}M_{2(k-j)}.
$$

### What it established

The numerical residuals between the two sides are at working-precision noise throughout the tested $N$-range, confirming the exact algebraic identity.

This is the key simplification in the tail calculation: the apparently complicated pairwise spectral interaction is exactly a convolution of the even spectral moments.

Cell 37 therefore supplies the algebraic bridge from the coefficient representation of Cell 36 to an endpoint-jet representation.

### Status

Established — exact moment-convolution identity.

---

## Cell 38 — closed form for the complete tail hierarchy

### Intended purpose

Cell 38 performs the remaining algebraic cancellation and rewrites the exact Cell-36 coefficients in terms of endpoint derivatives.

For $k\ge1$, the Cell-36 coefficient reduces to

$$
A_k = \frac{4}{L}\kappa^{2k} \left[ \sum_{j=0}^{k} M_{2j}M_{2(k-j)} + \sqrt2\thinspace v_0M_{2k} \right].
$$

Define

$$
D_0:=T_v(0), \qquad D_k:=T_v^{(2k)}(0) \quad(k\ge1).
$$

Then

$$
A_k = \frac{2}{L}(-1)^k \sum_{j=0}^{k}D_jD_{k-j}, \qquad k\ge1.
$$

Together with

$$
A_0=\frac{2}{L}D_0^2,
$$

this gives the complete $\text{finite-}N$ inverse-power hierarchy in terms of the even endpoint jet.

Equivalently,

$$
A_k = \frac{4}{L}(-1)^k \left[ T_v(0)T_v^{(2k)}(0) + \frac12 \sum_{j=1}^{k-1} T_v^{(2j)}(0) T_v^{(2k-2j)}(0) \right].
$$

### What it established

Cell 38 establishes an exact $\text{finite-}N$ algebraic identity: every inverse-power coefficient of the reduced Archimedean tail is a quadratic convolution of the even endpoint derivatives of $T_v$.

This result is independent of the ground-state property. It holds for any finite coefficient vector $v$ for which the finite Fourier representation is defined.

Consequently, the endpoint suppression observed numerically in Cell 35 has an exact analytical interpretation: small endpoint derivatives produce small coefficients in the $\text{large-}r$ expansion through the quadratic convolution above.

No $N\to\infty$ limit is assumed, and no numerical quadrature is involved.

### Status

Major established analytical result — exact closed form for the complete $\text{finite-}N$ Archimedean tail hierarchy.

---

# Current research state after Cells 24–38

*Updated 4 September 2026.*

The investigation has now passed through three distinct stages:

$$
\text{dictionary validation}
\quad\longrightarrow\quad
\text{numerical tail forensics}
\quad\longrightarrow\quad
\text{exact finite-}N\text{ tail analysis}.
$$

The important distinction between mathematical structure and numerical evaluation remains.

For finite $N$,

$$
K_{\mathrm{Fourier}}(v,r,L) = (1-\cos(rL))R_v(r),
$$

where $R_v(r)$ is rational and

$$
R_v(r) \sim \sum_{k\ge0}\frac{A_k}{r^{2k+2}}.
$$

The coefficients now have the exact representation

$$
A_0=\frac{2}{L}D_0^2,
$$

and, for $k\ge1$,

$$
A_k = \frac{2}{L}(-1)^k \sum_{j=0}^{k}D_jD_{k-j},
$$

where

$$
D_k=T_v^{(2k)}(0).
$$

Thus the $\text{finite-}N$ tail is no longer an opaque numerical phenomenon. Its complete asymptotic coefficient structure is known exactly.

The earlier extreme-range numerical experiments remain useful historically, but they are no longer the appropriate primary route for understanding the tail.

---

# Current priorities

## 1. Understand the endpoint jet of the ground state as $N\to\infty$

The exact hierarchy reduces the remaining problem to the behaviour of

$$
D_k(N)=T_{v_N}^{(2k)}(0)
$$

for the sequence of Galerkin ground states $v_N$.

The key question is now whether the observed suppression of

$$
D_0(N),D_1(N),D_2(N),\ldots
$$

has a systematic asymptotic structure.

In particular, we should distinguish:

* suppression of each fixed $D_k(N)$;
* the rate of that suppression with $N$;
* possible relationships between different $D_k(N)$;
* whether the endpoint jet approaches zero in a controlled functional sense.

No $N\to\infty$ law has yet been established.

## 2. Determine whether the endpoint hierarchy yields a usable tail bound

The exact coefficient formula suggests a route to a controlled remainder estimate.

For finite $N$, one ultimately wants a bound of the form

$$
\left| R_v(r) - \sum_{k=0}^{K} \frac{A_k}{r^{2k+2}} \right| \le \frac{C_{K,N}}{r^{2K+4}},
$$

with an explicit and preferably tractable $C_{K,N}$.

Combined with bounds on $h_+(r)$, this would turn the asymptotic expansion into a rigorous estimate for the Archimedean tail integral.

## 3. Look for a generating-function or boundary-operator interpretation

The convolution

$$
\sum_{j=0}^{k}D_jD_{k-j}
$$

is characteristic of a product of two formal power series.

This suggests that the entire hierarchy may admit a compact generating-function representation rather than being treated coefficient-by-coefficient.

Such a representation could expose positivity, cancellation, or remainder structure that is difficult to see from the individual $A_k$.

## 4. Keep $N$ and $T$ separate

The two limiting questions remain distinct:

$$
T\to\infty
$$

for the Archimedean integral at fixed finite $N$, and

$$
N\to\infty
$$

for the Galerkin sequence itself.

The exact $\text{finite-}N$ tail hierarchy substantially clarifies the first problem, but it does not by itself justify exchanging these limits.

## 5. Return to numerical integration only when analytically justified

Cell 31 demonstrated that extreme-range unsubdivided oscillatory quadrature can give apparently stable but incorrect results.

The exact factorisation and coefficient hierarchy now provide the appropriate analytical basis for any future numerical tail computation. Brute-force extension of the old Cell-30 calculation should therefore not be treated as the primary research route.

---

# Updated major historical arc

```
Cells 0–4
    Initial reconstruction
    ↓
Cell 5
    Archimedean discrepancy discovered
    ↓
Cells 6–20
    Source, coordinate, category and quadratic-form forensics
    ↓
Cell 21
    Independent brute-force quadratic validation
    ↓
Cells 22–23
    Analytic elimination of inner integral
    + optimised Archimedean implementation
    ↓
Cell 24
    Finite-T convergence map
    ↓
Cell 25
    Historical finite-T cross-check + extension
    ↓
Cell 26
    Long-range forensic tail
    ↓
Cell 27
    Pointwise tail anatomy / phase structure
    ↓
Cell 28
    Direct signed interval integration
    ↓
Cell 29
    Dyadic integrated-tail scaling
    ↓
Cell 30
    Extreme-range asymptotic hypothesis
    [quantitative interpretation later invalidated]
    ↓
Cell 31
    Quadrature forensic
    [precision exonerated; interval resolution identified]
    ↓
Cell 32
    Exact common oscillatory factor
    + analytical r^-2 tail
    + suppressed leading coefficient
    ↓
Cell 33
    Initial N-dependence survey
    ↓
Cell 34
    Systematic N-scan
    ↓
Cell 35
    Endpoint jets and spectral moments
    ↓
Cell 36
    Exact finite-N tail coefficients
    ↓
Cell 37
    Exact moment-convolution identity
    ↓
Cell 38
    Exact endpoint-jet hierarchy
    ↓
Current
    Ground-state endpoint-jet asymptotics
    + rigorous tail bounds
    + generating-function structure
```

# Current status summary

*Updated 4 September 2026.*

At the current stage:

* The finite Fourier/zero-side dictionary is substantially audited.
* The canonical/full coordinate distinction is understood.
* The prime-side dictionary has been independently audited.
* The Archimedean source has been independently audited.
* The genuinely quadratic $K_v$ construction has been established.
* The historical linear `G_complex` / current `sum_v_G` distinction is understood and preserved.
* Cell 21 provides an independent brute-force validation of the corrected Archimedean quadratic calculation.
* Cells 22–23 establish an efficient analytic implementation.
* Cells 24–26 establish the long-range $\text{finite-}T$ investigation.
* Cell 27 identifies strong $rL$-dependent oscillatory structure.
* Cells 28–30 investigate the integrated tail, while Cell 31 establishes that the extreme-range `mp.quad` results are not converged with respect to interval resolution.
* Cell 32 analytically identifies the exact common factor $1-\cos(rL)$ and the leading $r^{-2}$ coefficient.
* Cells 33–34 establish that the leading endpoint suppression is a systematic $N$-dependent phenomenon worth investigating further.
* Cell 35 identifies corresponding suppression in higher even endpoint derivatives and spectral moments.
* Cell 36 derives the exact $\text{finite-}N$ inverse-power coefficients algebraically.
* Cell 37 proves the moment-convolution identity that reorganises the pairwise spectral terms.
* Cell 38 reduces the entire tail hierarchy to an exact quadratic convolution of the even endpoint jet.

The central research question has therefore shifted again:

> **What endpoint-jet structure is imposed on the Galerkin ground state as $N\to\infty$, and can that structure provide a uniform or $N$-dependent bound on the Archimedean tail?**

This is now the natural bridge from the exact $\text{finite-}N$ analysis to the $\text{large-}N$ problem.

