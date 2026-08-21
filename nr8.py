import inspect
import connes_cvs.operator as op

for name in [
    "psi_prime",
    "psi_pole",
    "psi_arch",
    "psi_prime_deriv",
    "psi_pole_deriv",
    "psi_arch_deriv",
]:
    fn = getattr(op, name)
    print(f"\n=== {name} ===")
    print(inspect.signature(fn))
    print(inspect.getsource(fn))
