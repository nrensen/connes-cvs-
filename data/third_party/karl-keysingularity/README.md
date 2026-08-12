# karl-keysingularity c=13 reproduction

This directory preserves, verbatim, the machine-readable artifact submitted by
GitHub user [karl-keysingularity](https://github.com/karl-keysingularity) in
[issue #3](https://github.com/akivag613/connes-cvs-/issues/3). The contributor
explicitly permitted its inclusion in the package validation set.

The artifact reports a `connes-cvs` **v0.2.2** run at
`c=13, N=100, T=400, dps=80` under native Windows and Python 3.11.9. The
contributor reported a Ryzen 9 9950X and python-flint assembly, followed by
zero extraction using the corrected full-precision `mpmath.log(13)` input.

Maintainer verification against `tests/reference_values.json` finds:

- `lambda_min`: same sign and decimal exponent, with **22 leading significant
  digits** matching the 25-digit reference;
- `gamma1_error`: same sign and decimal exponent, with **all 20 significant
  digits carried by the reference** matching.

The comparison stops at the shorter supplied mantissa and never invents digits
by zero-padding. The raw artifact SHA-256 is
`5dddabbce30b4a4fa1a88f8ce34d82bb0d6b07e78801f5c36e54d7f7c428c79c`.

Limitations: this is a v0.2.2 result, not a v0.3.0 rerun. The artifact does not
record its exact dependency versions, command line, or source hash; timings are
specific to the contributor's environment. Digits beyond the 22/20 comparison
above have not been independently validated here. Credit is to the public
GitHub handle only; no real-world identity is inferred.
