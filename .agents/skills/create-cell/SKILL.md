---
name: create-cell
description: Procedural authoring guide, code templates, numerical recipes, and quality checklists for creating and updating computational cell scripts (cell*.py) in the Connes–CvS Galerkin project.
---

# Computational Cell Authoring Skill (`create-cell`)

This skill defines the canonical architectural standards, code skeletons, numerical recipes, and pre-flight quality checklists for authoring computational scripts (`cell*.py`) in the Connes–CvS Galerkin project.

---

## 1. Core Operating Constraints (Non-Negotiable)

1. **Zero Local Execution:**
   The agent must **NEVER** run Python scripts locally on the agent workspace. All computational cells must be written as self-contained, high-precision Python scripts ready for external execution on user compute nodes.
2. **Dispassionate, Dry Output:**
   Script outputs (`cell*.out`) must report computed values, matrix norms, eigenvalues, condition numbers, residuals, and convergence tables with absolute objectivity. No narrative hype, speculative analytical claims, or victory declarations in output strings.
3. **Deterministic Termination Sentinel:**
   Every cell script must terminate cleanly with the exact three-line sentinel:
   ```python
   print("=" * 80)
   print("CELL <ID> EXECUTION COMPLETE")
   print("=" * 80)
   ```

---

## 2. Strategic Alignment & Documentation Protocol

Before writing any new cell script, agents must ensure it satisfies the **Strategic Roadmap Traceability Principle**:
1. **Strategic Gate Alignment:** Every cell script must explicitly declare in its header docstring which of the five **Mathematical Gates** in [ROADMAP.md](file:///c:/data/github/connes-cvs-/ROADMAP.md) it serves (e.g. *Gate 1: Finite-N Spectral Mechanism & Asymptotic Tail Extinction*).
2. **Target Statement & Falsification Criterion:** The docstring must state the exact theorem, lemma, or hypothesis under test, along with the precise quantitative threshold for falsification or confirmation.
3. **Subordinate Diagnostic Role:** Cells are diagnostic and empirical probes, **never** standalone milestones. They do not replace mathematical theorems.
4. **Companion Analytical Note Protocol (`cell<n>.md`):**
   All theoretical derivations, operator identities, tested hypotheses, and intermediate reductions supporting a computational cell must be placed directly in a dedicated companion note named `cell<n>.md` (e.g. `cell101.md`). Ad-hoc notes with disparate names are prohibited.
   - $\text{cell note } (\text{cell}\langle n\rangle\text{.md}) = \text{what we are thinking}$ (exploratory, working reductions).
   - $\text{manuscripts } (\text{Paper-NR1, NR2}) = \text{what we believe is mathematically worth claiming}$ (settled theorems promoted after review).
5. **Notebook Logging in `cell_history_map.md`:** Once a cell is executed and its `.out` file analyzed, its mathematical rationale, computed findings, and refuted hypotheses must be logged in [cell_history_map.md](file:///c:/data/github/connes-cvs-/cell_history_map.md), with explicit links to `cell<n>.md`, **NOT** added as narrative clutter to [ROADMAP.md](file:///c:/data/github/connes-cvs-/ROADMAP.md).

---

## 3. Canonical Script Layout

Every cell script must follow this modular structure:

```
1. Header Docstring (Cell ID, Purpose, Target Gate, Tested Proposition, Falsification Criterion)
2. Precision & Parameter Configuration (mpmath dps, c, L, T, N_LIST)
3. Mathematical Utilities & Numerical Recipes (quadratures, basis functions)
4. Canonical Parity Projection Matrix (E: R^{N+1} -> R^{2N+1})
5. Problem-Specific Operator Construction & Eigensolvers
6. Dimension Loop & Structured Diagnostics (Guarded against ZeroDivisionError)
7. Multi-Dimension Synthesis Table
8. Clean Completion Sentinel
```

---

## 3. Standard Code Skeletons & Numerical Recipes

### 3.1 Precision & Global Configuration
```python
import time
import mpmath as mp
from connes_cvs import build_galerkin_matrix
from cell import get_galerkin_matrix

# Canonical precision baseline
mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

# Standard discrete dimension sweep
N_LIST = [4, 8, 12, 16, 20, 24]
```

### 3.1.0 Cache-Accelerated Matrix Retrieval (`cell.py`)
To avoid redundant, expensive numerical quadratures across sweeps and multiple diagnostic cells, construct the Galerkin matrix via `get_galerkin_matrix` from [`cell.py`](file:///c:/data/github/connes-cvs-/cell.py):
```python
# Returns (Q, metadata) with persistent content-addressed JSON caching and sub-dimension reuse
Q_full, _ = get_galerkin_matrix(
    c=C_PARAM,
    N=N,
    T=T_PARAM,
    dps=GROUND_DPS,
    verbose=False,
)
```
- **Persistent Caching:** Cached entries are stored in `.cell_cache/galerkin_matrix/` and load in < 10 ms.
- **Sub-Dimension Reuse:** If a cached entry for $N' \ge N$ exists for the same parameters, $Q$ is assembled instantly without evaluating any quadratures.
- **Incremental Evaluation:** When extending a dimension sweep to higher $N$, only missing basis points $n \in [N_{\mathrm{prev}} + 1, N]$ are evaluated.

### 3.1.1 `mp.mpf` Precision, Literal Construction & String Formatting
1. **Explicit String Wrapping for Constants:**
   Never pass raw floating-point numbers to `mp.mpf` if precision loss can occur. Always wrap constants in strings:
   - ❌ Incorrect: `mp.mpf(0.5)`, `mp.mpf(1e-35)` (causes binary float rounding prior to mpmath conversion)
   - ✅ Correct: `mp.mpf('0.5')`, `mp.mpf('1e-35')`

2. **No Python Float Format Specifiers on `mp.mpf` Objects:**
   `mpmath.mpf` objects do **NOT** support standard Python float format specifiers like `f"{val:.2f}"` or `f"{val:.6e}"` directly. Doing so raises `TypeError: unsupported format string passed to mpf.__format__`.
   - ❌ Incorrect: `shift_pct = f"{rec['rel_diff_pct']:.2f}%"` (crashes with TypeError)
   - ✅ Correct (via `mp.nstr`): `shift_pct = mp.nstr(rec['rel_diff_pct'], 4) + "%"`
   - ✅ Correct (via explicit float cast): `shift_pct = f"{float(rec['rel_diff_pct']):.2f}%"`
   - Always use `mp.nstr(val, digits)` for table entries and scientific outputs.

### 3.2 Canonical Parity Projection Matrix ($v$-Basis)
Always ensure calculations on the even sector operate in the canonical $(N+1)$-dimensional $v$-basis rather than the full $(2N+1)$-dimensional exponential basis:
```python
def canonical_even_basis(N: int) -> mp.matrix:
    """
    Construct the (2N+1) x (N+1) orthonormal matrix E mapping canonical v-basis
    v in R^{N+1} to the full exponential basis c in R^{2N+1}: c = E v.
    """
    dim_full = 2 * N + 1
    dim_can = N + 1
    E = mp.matrix(dim_full, dim_can)

    E[N, 0] = mp.mpf(1)
    inv_sqrt2 = 1 / mp.sqrt(2)
    for m in range(1, dim_can):
        E[N + m, m] = inv_sqrt2
        E[N - m, m] = inv_sqrt2

    return E

# Usage:
# Q_weil = E.T * Q_full * E
# Q_weil = mp.mpf('0.5') * (Q_weil + Q_weil.T)
```

### 3.3 Golub–Welsch Gauss–Legendre Quadrature
When continuous integration over an interval $[a, b]$ is required (e.g. against $h_+(r)$):
```python
def gauss_legendre_nodes_weights(order: int, a: mp.mpf, b: mp.mpf) -> tuple[list[mp.mpf], list[mp.mpf]]:
    """
    Compute Gauss-Legendre quadrature nodes and weights on [a, b] using
    the Golub-Welsch tridiagonal eigenvalue method at current mpmath dps.
    """
    J = mp.matrix(order, order)
    for i in range(order - 1):
        k = i + 1
        b_k = mp.mpf(k) / mp.sqrt(4 * k * k - 1)
        J[i, i + 1] = b_k
        J[i + 1, i] = b_k

    nodes_std, V = mp.eigsy(J)

    mid = (b + a) / 2
    half_width = (b - a) / 2

    nodes = []
    weights = []
    for i in range(order):
        x_i = nodes_std[i]
        w_i = 2 * (V[0, i] ** 2)
        nodes.append(mid + half_width * x_i)
        weights.append(half_width * w_i)

    return nodes, weights
```

### 3.4 Fourier Basis Amplitudes & Removable Singularities
Fourier basis functions $\Phi_m(r)$ possess removable singularities at $r = \pm a_m$ ($a_m = 2\pi m / L$). Always expand them locally via Taylor series to avoid catastrophic division-by-zero or loss of precision:
```python
def phi_basis(m: int, r: mp.mpf, L: mp.mpf) -> mp.mpf:
    """Canonical Fourier basis amplitude phi_m(r) on R."""
    half_L = L / 2
    if m == 0:
        if abs(r) < mp.mpf('1e-25'):
            return mp.sqrt(L)
        return (2 / mp.sqrt(L)) * (mp.sin(r * half_L) / r)

    a_m = 2 * mp.pi * m / L
    diff_pos = r - a_m
    diff_neg = r + a_m

    # Removable singularity at r = a_m
    if abs(diff_pos) < mp.mpf('1e-12'):
        u = diff_pos * half_L
        sinc_u = 1 - (u**2) / 6 + (u**4) / 120 - (u**6) / 5040
        geom = (a_m + diff_pos) / (2 * a_m + diff_pos)
        sgn = (-1) ** m
        return (2 * mp.sqrt(2) / mp.sqrt(L)) * (sgn * half_L * geom * sinc_u)

    # Removable singularity at r = -a_m
    if abs(diff_neg) < mp.mpf('1e-12'):
        u = diff_neg * half_L
        sinc_u = 1 - (u**2) / 6 + (u**4) / 120 - (u**6) / 5040
        geom = (-a_m + diff_neg) / (-2 * a_m + diff_neg)
        sgn = (-1) ** m
        return (2 * mp.sqrt(2) / mp.sqrt(L)) * (sgn * half_L * geom * sinc_u)

    denom = r**2 - a_m**2
    return (2 * mp.sqrt(2) / mp.sqrt(L)) * (r * mp.sin(r * half_L) / denom)
```

### 3.5 Generalized Eigensolver with Nullspace Audit
Before solving $A x = \lambda B x$ where $B \succeq 0$, always audit the singular spectrum and numerical nullspace of $B$:
```python
def solve_generalized_eigenproblem(
    A: mp.matrix,
    B: mp.matrix,
    tol_rel: mp.mpf = mp.mpf('1e-35'),
) -> tuple[list[mp.mpf], mp.matrix, list[mp.mpf], int]:
    dim = A.rows
    b_vals, V_b = mp.eigsy(B)

    # Sort descending
    idx_b = sorted(range(dim), key=lambda i: b_vals[i], reverse=True)
    B_evals = [b_vals[i] for i in idx_b]
    V_sorted = mp.matrix(dim, dim)
    for j, i in enumerate(idx_b):
        for r in range(dim):
            V_sorted[r, j] = V_b[r, i]

    sigma_max = B_evals[0]
    rank = 0
    for s in B_evals:
        if s > tol_rel * sigma_max and s > 0:
            rank += 1
        else:
            break

    null_dim = dim - rank

    # Whitening transformation S = V_r * diag(1 / sqrt(sigma_i))
    S = mp.matrix(dim, rank)
    for j in range(rank):
        inv_sqrt = 1 / mp.sqrt(B_evals[j])
        for r in range(dim):
            S[r, j] = V_sorted[r, j] * inv_sqrt

    # Projected standard symmetric problem M = S^T * A * S
    M = S.T * A * S
    M = mp.mpf('0.5') * (M + M.T)
    m_vals, Y = mp.eigsy(M)

    idx_m = sorted(range(rank), key=lambda i: m_vals[i])
    evals = [m_vals[i] for i in idx_m]

    # Pull back to original coordinates x = S * y and normalize ||x||_2 = 1
    evecs = mp.matrix(dim, rank)
    for j, i in enumerate(idx_m):
        y_col = Y[:, i]
        x_col = S * y_col
        norm_x = mp.sqrt(sum(x_col[r, 0] ** 2 for r in range(dim)))
        for r in range(dim):
            evecs[r, j] = x_col[r, 0] / norm_x

    return evals, evecs, B_evals, null_dim
```

---

## 5. Pre-Flight Quality Checklist

Before committing any cell script, audit against this checklist:

- [ ] **Strategic Gate Alignment:** Does the header docstring explicitly cite which of the 5 Mathematical Gates in `ROADMAP.md` is targeted, along with the precise mathematical statement and falsification criterion?
- [ ] **Singularity & Division-by-Zero Protection:** Are all denominator differences guarded? (e.g. In spectral product or tail analyses, are local/focus modes such as $\ell = j$ and $\ell = j+1$ strictly guarded or excluded before evaluating differences like $|E_{j+1} - E_\ell|$?)
- [ ] **F-String Bracket Escaping:** Are mathematical sets in f-strings escaped with double braces?
  - ❌ Incorrect: `print(f"Energy in {e_0, e_1} = {val}")` (causes `NameError: name 'e_0' is not defined`)
  - ✅ Correct: `print(f"Energy in {{e_0, e_1}} = {val}")`
- [ ] **`mp.mpf` String Formatting:** Are all `mp.mpf` objects formatted without Python float specifiers (like `:.2f`)?
  - ❌ Incorrect: `f"{rec['shift_pct']:.2f}%"` (crashes with `TypeError: unsupported format string passed to mpf.__format__`)
  - ✅ Correct: `mp.nstr(val, digits)` or `f"{float(val):.2f}%"`
- [ ] **Basis Consistency:** Is the script operating strictly in the intended basis ($v$-basis of size $N+1$ vs full-space exponential basis of size $2N+1$)?
- [ ] **Removable Singularities:** Are all continuous denominators protected against algebraic zeroes via Taylor series expansions?
- [ ] **Eigenpair Residual Checking:** Are computed eigenvalues and eigenvectors checked for numerical residuals (e.g. $\|A x - \lambda B x\|_2$)?
- [ ] **Epistemic Labeling:** Does the script avoid describing numerical dominance tests as "rigorous certificates"?
- [ ] **Dispassionate Output:** Are all printed labels and messages purely technical and descriptive?
- [ ] **Completion Sentinel:** Does the script conclude with the exact three-line sentinel?
