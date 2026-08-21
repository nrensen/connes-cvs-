import inspect
import connes_cvs.operator as op

print(inspect.signature(op.build_galerkin_matrix))
print(inspect.getsource(op.build_galerkin_matrix))
