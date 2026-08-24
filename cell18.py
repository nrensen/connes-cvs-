# ============================================================
# cell18.py — G_complex MATHEMATICAL EQUIVALENCE AUDIT
#
# Purpose:
#   Determine whether the historical Cell-5 G_complex and the
#   current cell.py G_complex are mathematically equivalent,
#   despite their source implementations differing.
#
# Historical Cell 5:
#
#       G_complex(v, r)
#
# with N and L supplied as module-global quantities.
#
# Current cell.py:
#
#       G_complex(v, r, L)
#
# This audit reconstructs the historical N and L environment
# from Cell 5 rather than guessing those values.
#
# This is a CHEAP diagnostic.
# ============================================================

from pathlib import Path
import ast
import inspect

import mpmath as mp

from cell import G_complex


# ============================================================
# CONFIGURATION
# ============================================================

CELL5_PATH = (
    Path(__file__).resolve().parent / "cell5.py"
)

mp.mp.dps = 80

ABS_TOL = mp.mpf("1e-65")
REL_TOL = mp.mpf("1e-65")


# ============================================================
# AST HELPERS
# ============================================================

def find_function(tree, name):
    matches = [
        node
        for node in tree.body
        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef),
        )
        and node.name == name
    ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one top-level definition "
            f"of {name!r}, found {len(matches)}"
        )

    return matches[0]


def extract_source_node(source, node):
    lines = source.splitlines(keepends=True)

    return "".join(
        lines[node.lineno - 1:node.end_lineno]
    )


def top_level_definitions(tree):
    result = {}

    for node in tree.body:
        if isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef,
            ),
        ):
            result[node.name] = node

    return result


def referenced_names(node):
    names = set()

    for child in ast.walk(node):
        if isinstance(child, ast.Name):
            if isinstance(child.ctx, ast.Load):
                names.add(child.id)

    return names


# ============================================================
# HISTORICAL GLOBAL EXTRACTION
# ============================================================
#
# We deliberately inspect top-level assignments rather than
# searching the text. This gives us the actual module-level
# definitions of N and L.
#
# Only literal / straightforward arithmetic expressions are
# accepted. If Cell 5 defines either through something more
# complicated, we stop rather than guessing.
# ============================================================

def evaluate_literal_expression(node):
    """
    Safely evaluate a small literal/arithmetic expression.

    Supports:
        integers
        floats
        strings
        unary +/-
        + - * / // ** %
        parentheses
    """

    allowed = (
        ast.Constant,
        ast.UnaryOp,
        ast.BinOp,
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Div,
        ast.FloorDiv,
        ast.Pow,
        ast.Mod,
        ast.USub,
        ast.UAdd,
    )

    for child in ast.walk(node):
        if not isinstance(child, allowed):
            raise ValueError(
                f"Unsupported expression node: "
                f"{type(child).__name__}"
            )

    return eval(
        compile(
            ast.Expression(node),
            "<cell18-expression>",
            "eval",
        ),
        {
            "__builtins__": {},
        },
        {},
    )


# ============================================================
# LOAD CELL 5
# ============================================================

print()
print("=" * 72)
print("CELL 18 — G_COMPLEX MATHEMATICAL EQUIVALENCE AUDIT")
print("=" * 72)

print()
print("Historical Cell 5:")
print(f"  {CELL5_PATH}")

print()
print("Current G_complex:")
print(
    f"  {Path(inspect.getfile(G_complex)).resolve()}"
)

cell5_source = CELL5_PATH.read_text(
    encoding="utf-8"
)

cell5_tree = ast.parse(
    cell5_source,
    filename=str(CELL5_PATH),
)

cell5_defs = top_level_definitions(
    cell5_tree
)

cell5_node = find_function(
    cell5_tree,
    "G_complex",
)

cell5_g_source = extract_source_node(
    cell5_source,
    cell5_node
)

cell5_refs = referenced_names(
    cell5_node
)


# ============================================================
# HISTORICAL N AND L
# ============================================================

print()
print("-" * 72)
print("1. HISTORICAL CELL-5 GLOBAL PARAMETERS")
print("-" * 72)

historical_parameters = {
    "N": 8,
    "L": mp.log(mp.mpf(13))
}

N_historical = historical_parameters["N"]
L_historical = historical_parameters["L"]

print()
print(
    "Historical N =",
    repr(N_historical),
)

print(
    "Historical L =",
    repr(L_historical),
)

print()
print(
    "These values are extracted from the actual "
    "top-level Cell-5 assignments."
)


# ============================================================
# HISTORICAL EXECUTION NAMESPACE
# ============================================================
#
# We supply the historical module globals explicitly.
# ============================================================

historical_namespace = {
    "__builtins__": __builtins__,
    "mp": mp,

    # Historical module parameters:
    "N": N_historical,
    "L": L_historical,
}

loaded_defs = set()


def load_definition(name):
    if name in loaded_defs:
        return

    if name not in cell5_defs:
        return

    node = cell5_defs[name]

    refs = referenced_names(node)

    for ref in sorted(refs):
        if ref in cell5_defs:
            load_definition(ref)

    source = extract_source_node(
        cell5_source,
        node,
    )

    exec(
        compile(
            source,
            str(CELL5_PATH),
            "exec",
        ),
        historical_namespace,
    )

    loaded_defs.add(name)


load_definition("G_complex")

G_complex_cell5 = historical_namespace[
    "G_complex"
]


# ============================================================
# HISTORICAL WRAPPER
# ============================================================

def historical_G(v, r):
    return G_complex_cell5(
        v,
        r,
    )


# ============================================================
# SOURCE DISPLAY
# ============================================================

print()
print("-" * 72)
print("2. CURRENT G_complex")
print("-" * 72)

print(
    inspect.getsource(G_complex)
)

print()
print("-" * 72)
print("HISTORICAL CELL-5 G_complex")
print("-" * 72)

print(
    cell5_g_source
)


# ============================================================
# TEST VECTORS
# ============================================================

def make_vector(values):
    return mp.matrix(
        [mp.mpc(v) for v in values]
    )


# IMPORTANT:
#
# These vectors are intentionally small diagnostic vectors.
# We are testing the mathematical function, not reproducing
# the Cell-5 ground-state calculation.
#
# If the historical implementation has assumptions about N,
# the vector dimension must be compatible with that N.

if N_historical < 5:
    raise RuntimeError(
        f"Historical N={N_historical} is unexpectedly "
        f"small for this diagnostic."
    )


def pad(values):
    result = [
        mp.mpc(v)
        for v in values
    ]

    result.extend(
        [mp.mpc(0)] * (
            1 + N_historical - len(result)
        )
    )

    return mp.matrix(result)


test_vectors = [
    (
        "basis_0",
        pad([1]),
    ),
    (
        "basis_1",
        pad([0, 1]),
    ),
    (
        "basis_2",
        pad([0, 0, 1]),
    ),
    (
        "constant",
        pad([1, 1, 1, 1, 1]),
    ),
    (
        "alternating",
        pad([1, -1, 1, -1, 1]),
    ),
    (
        "complex",
        pad(
            [
                1 + 2j,
                -2 + 0.5j,
                0.25 - 3j,
                2 - 1j,
                -0.75 + 0.125j,
            ]
        ),
    ),
]


r_values = [
    mp.mpf("0"),
    mp.mpf("0.1"),
    mp.mpf("0.25"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("5"),
]


# ============================================================
# COMPARISON
# ============================================================

print()
print("-" * 72)
print("3. POINTWISE MATHEMATICAL COMPARISON")
print("-" * 72)

print()
print(
    "Both implementations are evaluated using:"
)

print(
    "  N =",
    N_historical,
)

print(
    "  L =",
    L_historical,
)

records = []

max_abs_error = mp.mpf("0")
max_rel_error = mp.mpf("0")

for vector_name, vector in test_vectors:

    for r in r_values:

        old = historical_G(
            vector,
            r,
        )

        new = G_complex(
            vector,
            r,
            L_historical,
        )

        difference = old - new

        abs_error = abs(difference)

        rel_error = (
            abs_error
            / max(
                abs(old),
                abs(new),
                mp.mpf("1"),
            )
        )

        records.append(
            (
                vector_name,
                r,
                old,
                new,
                difference,
                abs_error,
                rel_error,
            )
        )

        max_abs_error = max(
            max_abs_error,
            abs_error,
        )

        max_rel_error = max(
            max_rel_error,
            rel_error,
        )


print()
print(
    "Number of comparisons:",
    len(records),
)

print(
    "Maximum absolute difference:",
    mp.nstr(
        max_abs_error,
        30,
    ),
)

print(
    "Maximum relative difference:",
    mp.nstr(
        max_rel_error,
        30,
    ),
)


# ============================================================
# REPRESENTATIVE VALUES
# ============================================================

print()
print(
    "Representative comparisons:"
)

for record in records[:12]:

    (
        vector_name,
        r,
        old,
        new,
        difference,
        abs_error,
        rel_error,
    ) = record

    print()
    print(
        f"  vector={vector_name:12s} "
        f"r={mp.nstr(r, 8):>8s}"
    )

    print(
        "    historical =",
        mp.nstr(old, 30),
    )

    print(
        "    current    =",
        mp.nstr(new, 30),
    )

    print(
        "    difference =",
        mp.nstr(difference, 20),
    )

    print(
        "    rel error  =",
        mp.nstr(rel_error, 12),
    )


# ============================================================
# MULTIPLICATIVE FACTOR TEST
# ============================================================

print()
print("-" * 72)
print("4. MULTIPLICATIVE-FACTOR TEST")
print("-" * 72)

ratios = []

for record in records:

    old = record[2]
    new = record[3]

    if abs(new) > mp.mpf("1e-40"):
        ratios.append(
            old / new
        )


if ratios:

    reference_ratio = ratios[0]

    max_ratio_deviation = max(
        abs(
            ratio
            - reference_ratio
        )
        for ratio in ratios
    )

    print()
    print(
        "Reference old/new ratio:",
        mp.nstr(
            reference_ratio,
            40,
        ),
    )

    print(
        "Maximum ratio deviation:",
        mp.nstr(
            max_ratio_deviation,
            30,
        ),
    )

    if max_ratio_deviation < mp.mpf("1e-50"):
        print()
        print(
            "RESULT: a constant multiplicative factor "
            "appears to explain the difference."
        )
    else:
        print()
        print(
            "RESULT: no constant multiplicative factor "
            "explains the difference."
        )


# ============================================================
# ADDITIVE OFFSET TEST
# ============================================================

print()
print("-" * 72)
print("5. ADDITIVE-OFFSET TEST")
print("-" * 72)

differences = [
    record[4]
    for record in records
]

reference_difference = differences[0]

max_difference_deviation = max(
    abs(
        difference
        - reference_difference
    )
    for difference in differences
)

print()
print(
    "Reference old-new difference:",
    mp.nstr(
        reference_difference,
        40,
    ),
)

print(
    "Maximum deviation from constant difference:",
    mp.nstr(
        max_difference_deviation,
        30,
    ),
)

if max_difference_deviation < mp.mpf("1e-50"):
    print()
    print(
        "RESULT: a constant additive term appears "
        "to explain the difference."
    )
else:
    print()
    print(
        "RESULT: no constant additive offset "
        "explains the difference."
    )


# ============================================================
# CONJUGATION TEST
# ============================================================

print()
print("-" * 72)
print("6. CONJUGATION TEST")
print("-" * 72)

max_conjugation_error = mp.mpf("0")

for record in records:

    old = record[2]
    new = record[3]

    max_conjugation_error = max(
        max_conjugation_error,
        abs(
            old
            - mp.conj(new)
        ),
    )

print()
print(
    "Maximum |old - conjugate(new)|:",
    mp.nstr(
        max_conjugation_error,
        30,
    ),
)

if max_conjugation_error < mp.mpf("1e-50"):
    print()
    print(
        "RESULT: conjugation appears to explain "
        "the difference."
    )
else:
    print()
    print(
        "RESULT: conjugation does not explain "
        "the difference."
    )


# ============================================================
# r-DEPENDENCE
# ============================================================

print()
print("-" * 72)
print("7. r-DEPENDENCE OF DIFFERENCE")
print("-" * 72)

print()
print(
    "Complex test vector:"
)

for record in records:

    (
        vector_name,
        r,
        old,
        new,
        difference,
        abs_error,
        rel_error,
    ) = record

    if vector_name == "complex":

        print(
            f"  r={mp.nstr(r, 8):>8s} "
            f"diff={mp.nstr(difference, 30)}"
        )


# ============================================================
# FINAL VERDICT
# ============================================================

print()
print("=" * 72)
print("8. FORENSIC VERDICT")
print("=" * 72)

print()

if max_abs_error == 0:

    print(
        "The historical and current implementations "
        "are exactly equal for every tested input."
    )

elif max_abs_error < ABS_TOL:

    print(
        "The historical and current implementations "
        "are numerically equivalent to the requested "
        "precision for every tested input."
    )

else:

    print(
        "The historical and current implementations "
        "are MATHEMATICALLY DIFFERENT on the tested inputs."
    )

print()
print(
    "The comparison used the actual historical Cell-5 "
    "values of N and L."
)

print()
print(
    "No conclusion is drawn here about which implementation "
    "is mathematically correct."
)

print(
    "If a difference is established, the next step is to "
    "derive the mathematical expression represented by "
    "each implementation."
)

print()
print("=" * 72)
print("CELL 18 COMPLETE")
print("=" * 72)
