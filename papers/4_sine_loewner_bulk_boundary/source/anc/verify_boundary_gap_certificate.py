#!/usr/bin/env python3
"""Verify a physical half-line gap enclosure with exact-input ball arithmetic.

This independently uses direct double-cell integrals, not the jump-sum
implementation used to discover the witness. No numerical quadrature,
eigenvalue solver, or saved floating-point result enters this check.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import flint
from flint import arb, acb, fmpq, ctx


def stamp():
    return subprocess.check_output(
        ['env', 'TZ=Asia/Jerusalem', 'date', '+%Y-%m-%d %H:%M:%S %Z'],
        text=True).strip()


def primitive(z):
    """Continuous primitive of log|z| on the real line."""
    return arb(0) if z.is_zero() else z*(abs(z).log()-1)


def certify(witness, bits):
    ctx.prec = bits
    ctx.threads = 1
    t = [fmpq(x) for x in witness['positive_mesh_exact_dyadics']]
    assert t[0] == 0 and t[-1] == 1
    assert all(a < b for a, b in zip(t, t[1:]))
    n = len(t)-1
    edges = [arb(x) for x in [x-1 for x in t[:-1]]+t]
    widths = [b-a for a, b in zip(edges, edges[1:])]
    assert all(h > 0 for h in widths)
    coeff = [acb(arb(fmpq(a)), arb(fmpq(b))) for a, b in zip(
        witness['coeff_real_exact_dyadics'], witness['coeff_imag_exact_dyadics'])]
    assert len(coeff) == 2*n
    values = [c/h.sqrt() for c, h in zip(coeff, widths)]
    lam = arb(fmpq(witness['lambda_exact_dyadic']))

    def projected_hilbert(a, b):
        # Integrate log|x-left_j| - log|x-right_j| directly over (a,b).
        integral = acb(0)
        for left, right, value in zip(edges, edges[1:], values):
            exact_integral = (primitive(b-left)-primitive(a-left)
                              -primitive(b-right)+primitive(a-right))
            integral += value*exact_integral
        return integral/(arb.pi()*(b-a).sqrt())

    h_inside = [projected_hilbert(a, b) for a, b in zip(edges, edges[1:])]
    q_inside = [(c+acb(0, 1)*h)/2 for c, h in zip(coeff, h_inside)]
    j_coeff = coeff[n:]+coeff[:n]
    norm = sum((c.real*c.real+c.imag*c.imag for c in coeff), arb(0))
    q_energy = sum(((c.conjugate()*q).real for c, q in zip(coeff, q_inside)), arb(0))
    cross = sum(((q.conjugate()*j).real for q, j in zip(q_inside, j_coeff)), arb(0))
    assert q_energy > 0

    outside = [fmpq(x) for x in witness['outside_positive_exact_edges']]
    assert outside[0] >= 1 and all(a < b for a, b in zip(outside, outside[1:]))
    projected_outside_energy = arb(0)
    for left, right in zip(outside, outside[1:]):
        a, b = arb(left), arb(right)
        for hp in (projected_hilbert(a, b), projected_hilbert(-b, -a)):
            projected_outside_energy += hp.real*hp.real+hp.imag*hp.imag

    # Bessel gives an outside-energy LOWER bound. Hilbert isometry then
    # gives a frequency-residual UPPER bound; no omitted-tail estimate is needed.
    upper_squared = (q_energy-projected_outside_energy/4
                     -2*lam*cross+lam*lam*norm)
    assert upper_squared > 0
    physical_upper = (upper_squared/q_energy).sqrt()
    lo, hi = [arb(fmpq(x)) for x in witness['target_interval']]
    checks = {
        'inside_positive_central_gap': bool(lo > 0 and hi < 1/arb(2).sqrt()),
        'strict_lower_margin': bool(physical_upper < lam-lo),
        'strict_upper_margin': bool(physical_upper < hi-lam),
        'physical_residual_below_stated_threshold': bool(
            physical_upper < arb(fmpq(witness.get('residual_threshold', '591/10000')))),
    }
    if witness['target_interval'] == ['81/200', '107/250']:
        checks.update({
            'printed_trial_lower': bool(lam > arb(fmpq(41635, 100000))),
            'printed_trial_upper': bool(lam < arb(fmpq(41637, 100000))),
            'printed_synthesis_lower': bool(q_energy > arb(fmpq(32783, 100000))),
            'printed_residual_squared_upper': bool(upper_squared < arb(fmpq(38641, 1000000000))),
        })
    assert all(checks.values()), checks
    return {'precision_bits': bits, 'norm_squared': str(norm),
            'synthesis_norm_squared': str(q_energy), 'mixed_inner_product': str(cross),
            'outside_projected_energy': str(projected_outside_energy),
            'frequency_residual_squared_upper': str(upper_squared),
            'physical_residual_upper': str(physical_upper), 'checks': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--witness', type=Path)
    args = parser.parse_args()
    if args.output:
        assert not args.output.exists(), 'Refusing to overwrite an output file'
    source = Path(__file__).resolve()
    fixture = args.witness or source.with_name('boundary_gap_witness.json')
    raw = fixture.read_bytes()
    witness = json.loads(raw)
    record = {'started': stamp(), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'witness_sha256': hashlib.sha256(raw).hexdigest(), 'python_flint': flint.__version__,
              'method': 'Direct double-cell Hilbert integrals, rigorous outside Bessel bound',
              'scope': 'At least one phase-zero physical half-line eigenvalue in the stated interval; no uniqueness or precise point value',
              'target_interval': witness['target_interval'],
              'runs': []}
    for bits in (192, 256):
        record['runs'].append(certify(witness, bits))
        print(stamp(), bits, 'bits: interval certificate PASS', flush=True)
    record['finished'] = stamp()
    record['status'] = 'PASS'
    text = json.dumps(record, indent=2)+'\n'
    if args.output:
        with args.output.open('x') as stream:
            stream.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
