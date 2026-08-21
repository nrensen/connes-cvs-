import inspect
import connes_cvs.operator as op

print(inspect.signature(op.extract_zeros))
print(inspect.getsource(op.extract_zeros))
