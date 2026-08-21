from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros
import mpmath as mp

# Seconds-level API smoke cell (not a paper benchmark)
Q = build_galerkin_matrix(c=13, N=8, T=60, dps=30)

# Diagonalize
lam_min, eigvec = compute_ground_state(Q)
print(f"\xce\xbb_min(c=13) = {mp.nstr(lam_min, 6)}")
# \xce\xbb_min(c=13) \xe2\x89\x88 4.43043e-23

# Extract the first detected Riemann zero.
# Preferred form (v0.3.0): pass the cutoff c and the package computes
# L = log(c) internally at full working precision.
zeros = extract_zeros(eigvec, c=13, n_zeros=1, dps=30)
# Equivalent, still supported: extract_zeros(eigvec, L=mp.log(13), ...).
# Passing L as a Python float carries only ~16 digits and caps the
# extraction accuracy near 1e-16; v0.3.0 emits a UserWarning.
print(f"\xce\xb3\xe2\x82\x81 detected = {mp.nstr(zeros[0]['gamma_detected'], 12)}")
print(f"|\xce\xb3\xe2\x82\x81 error|  = {mp.nstr(zeros[0]['error'], 4)}")
# \xce\xb3\xe2\x82\x81 detected \xe2\x89\x88 14.1347251417
# |\xce\xb3\xe2\x82\x81 error|  \xe2\x89\x88 2.52738e-17

