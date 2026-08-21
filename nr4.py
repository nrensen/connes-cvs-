import inspect
import connes_cvs.kernels as ker

print("=== S_hat_x ===")
print(inspect.signature(ker.S_hat_x))
print(inspect.getsource(ker.S_hat_x))

print("\n=== dS_hat_x_dx ===")
print(inspect.signature(ker.dS_hat_x_dx))
print(inspect.getsource(ker.dS_hat_x_dx))
