#!/usr/bin/env python3
"""Exact rational arithmetic for the whole-fiber trace-constant bound.

Analytic premises are derived in the accompanying full-fiber gauge lemma and constants remark in sections/21_trace_proof.tex.
This certifies their scalar inequalities, not the operator proof itself.
No quadrature, eigensolver, or floating-point input is used. Outputs are
protected by exclusive creation and must be supplied explicitly.
"""
import argparse
from datetime import datetime
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
from zoneinfo import ZoneInfo


def encoded(x):
    with localcontext() as ctx:
        ctx.prec = 40
        decimal = str(Decimal(x.numerator) / Decimal(x.denominator))
    return {"numerator": x.numerator, "denominator": x.denominator,
            "decimal_display_only": decimal}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    source = root / "sections/21_trace_proof.tex"
    pi_lower, pi_upper = F(31, 10), F(22, 7)
    sqrt2_upper = F(99, 70)
    kappa, log_trace = F(67, 10), F(15)
    kappa_squared_upper = 4 * (4 * pi_upper**2 / 3 - 2)
    logarithm_trace_upper = pi_upper * (kappa + 2 * sqrt2_upper) / 2
    gauge_upper = (kappa * (F(13, 12) + 79*pi_upper/48)
                   + log_trace * (F(7, 12) + 13*pi_upper/48)
                   + pi_upper/3 + pi_upper**2/48)

    r = pi_upper/2
    strip_upper = sum((r**j / (factorial(j)*(2*j+1)) for j in range(5)), F(0))
    strip_upper += r**5 / (factorial(5)*11) / (1-r/6)
    scalar_defect_upper = pi_upper*(3+sqrt2_upper)/24
    b1_upper = (3+F(7, 10))/pi_lower**2 + F(14, 3)/pi_lower
    b2_upper = 2/pi_lower**2 + F(15, 4)/pi_lower**4
    scalar_portion_upper = (4*strip_upper + 2 + 2/pi_lower**2
                            + scalar_defect_upper + b1_upper + 2*b2_upper)
    coarse_scalar_upper = 4*F(191,100)+2+2/pi_lower**2+F(3,5)+F(19,10)+F(1,2)
    checks = {
        "sqrt2_bound": sqrt2_upper**2 > 2,
        "kappa_below_67_over_10": kappa_squared_upper < kappa**2,
        "logarithm_trace_below_15": logarithm_trace_upper < 15,
        "gauge_cost_below_65": gauge_upper < 65,
        "strip_below_191_over_100": strip_upper < F(191,100),
        "scalar_defect_below_3_over_5": scalar_defect_upper < F(3,5),
        "B1_below_19_over_10": b1_upper < F(19,10),
        "B2_below_1_over_4": b2_upper < F(1,4),
        "scalar_portion_below_13": scalar_portion_upper < coarse_scalar_upper < 13,
        "assembled_trace_cost_below_143": 2*gauge_upper+scalar_portion_upper < 143,
        "exp_2_above_7": sum((F(2)**j/factorial(j) for j in range(6)), F(0)) > 7,
        "exp_7_over_10_above_2": sum((F(7,10)**j/factorial(j) for j in range(5)), F(0)) > 2,
        "exp_1_above_5_over_2": sum((F(1)/factorial(j) for j in range(4)), F(0)) > F(5,2),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    values = {
        "kappa_squared_upper": kappa_squared_upper,
        "logarithm_trace_upper": logarithm_trace_upper,
        "gauge_cost_upper": gauge_upper,
        "strip_upper": strip_upper,
        "scalar_defect_upper": scalar_defect_upper,
        "B1_upper": b1_upper,
        "B2_upper": b2_upper,
        "scalar_portion_upper": scalar_portion_upper,
        "coarse_scalar_portion_upper": coarse_scalar_upper,
        "assembled_trace_cost_upper": 2*gauge_upper+scalar_portion_upper,
    }
    result = {
        "status": "PASS_EXACT_RATIONAL_CONSTANT_BOUNDS",
        "recorded_at": datetime.now(ZoneInfo("Asia/Jerusalem")).strftime("%Y-%m-%d %H:%M:%S %Z"),
        "scope": "Exact arithmetic checks conditional on the analytic operator estimates in the full-fiber gauge lemma and constants remark in sections/21_trace_proof.tex; no numerical spectral claim.",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_path": str(source.relative_to(root)),
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "checks": checks,
        "values": {name: encoded(value) for name, value in values.items()},
    }
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
