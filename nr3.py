import inspect
import connes_cvs.operator as op
import connes_cvs.kernels as ker

print("=== operator ===")
print([x for x in dir(op) if not x.startswith("_")])

print("\n=== kernels ===")
print([x for x in dir(ker) if not x.startswith("_")])
