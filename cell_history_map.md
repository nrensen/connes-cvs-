# Connes–CvS exploratory cell history

**Repository:** `nrensen/connes-cvs-`  
**Historical snapshot audited:** commit `9aa29074b05a457c012dba7ec33cf2f5ee04f38c`  
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

$$G_v(r)=\sum_k v_kG_k(r)$$.

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

$$\text{completed-zeta Archimedean factor}
\to h_+(\tau)
\to \text{basis Fourier response}
\to S_x(\tau)
\to \psi_{\rm arch}(x)$$.

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

$$v
\to f_v(t)
\to F_v(\tau)
\to \text{translated test function}
\to \text{Weil quadratic form}$$.

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

$$v
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

$$\frac1\pi\int_0^T h_+(r)\int_0^L K_v(1-y/L)\cos(ry) \thinspace dy \thinspace dr$$.

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

# Cell 21 — planned next substantive cell

## Intended purpose

Cell 21 has not yet been committed.

Its intended role is to move from the historical discrepancy investigations to a clean, independent treatment of the **Archimedean quadratic dictionary**.

The principal path should use:
$T_v\to K_v\to\widehat g_v\to v^*Q_{\rm arch}v$,

rather than `sum_v_G` as a surrogate quadratic object.

The intended question is essentially:

> Does the independently constructed quadratic Archimedean functional agree with the Archimedean quadratic form represented by the Galerkin matrix?

Cell 17 and Cell 20 provide the immediate groundwork.

### Design constraint

Cell 21 should be deliberately resistant to the historical Cell-5/Cell-6 category error.

In particular, `sum_v_G` may appear as a separately labelled diagnostic if useful, but it should **not be part of the principal quadratic calculation**.

### Status

**Planned / open.**

---

# Major historical arc

The overall development can be summarised as:

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
Cell 20
    Corrected Archimedean quadratic audit
    ↓
Cell 20a
    Pole linear-vs-quadratic forensic audit
    ↓
Cell 21
    Correct Archimedean quadratic dictionary
    [planned]
```
## The two distinct historical errors

It is important that the repository history does not collapse these into one issue.

1. Canonical/full coordinate error

At various stages the full $(2N+1)$ representation and canonical
$(N+1)$ representation were confused.

This was investigated and resolved by Cells 13–15.

2. Linear  $G_v$ versus quadratic Weil functional

The historical Archimedean calculation treated $G_v(r)=\sum_k v_kG_k(r)$
as though it represented the quadratic Archimedean functional.

This was the deeper semantic error exposed through Cell 5 and revisited in Cell 6.

Cell 17 established the appropriate genuinely quadratic $K_v$ route.

These are **different errors** and should remain separately documented.

## `cell.py` refactor — historical significance

The `cell.py` refactor at commit `8b317297bb1361d37d64cdc747c656ee50918add` was not a mathematical correction.

It changed the vocabulary:

```text
F_vector          → sum_v_F
Fprime_vector     → sum_v_Fprime
G_complex_basis   → G_basis_complex
G_complex         → sum_v_G
```

The pre/post regression test confirmed that the numerical fingerprint was unchanged.

The subsequent migration commit `9aa29074b05a457c012dba7ec33cf2f5ee04f38c` mechanically updated the few cells that import these objects.

The purpose of this refactor is therefore historical as well as practical:

> **Make the mathematical category of an object more visible without changing what the historical calculations actually do.**

In particular, `sum_v_G` makes explicit that the operation is the coefficient-weighted sum
$\sum_k v_kG_k$, which is precisely the distinction that was obscured by the old name `G_complex`.

# Research-state summary

At the current stage:

* The finite Fourier/zero-side dictionary is substantially audited.
* The canonical/full coordinate distinction is understood.
* The prime-side dictionary has been independently audited.
* The Archimedean source itself has been independently audited.
* The genuinely quadratic `K_v` construction has been established and tested.
* The historical `G_complex` construction has been confirmed mathematically equivalent to the current `sum_v_G`.
* The historical Cell-5 and Cell-6 misuse of that linear construction remains preserved as historical record.
* Cell 20 provides the current corrected Archimedean calculation.
* Cell 20a separately investigates the pole linear/quadratic distinction.
* Cell 21 remains the natural next substantive synthesis of the Archimedean quadratic dictionary.
