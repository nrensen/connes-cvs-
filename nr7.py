import inspect
import connes_cvs.operator as op

print(inspect.signature(op._compute_psi_pair))
print(inspect.getsource(op._compute_psi_pair))
