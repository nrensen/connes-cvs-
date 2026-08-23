import mpmath as mp

from cell import FORENSIC_GROUND_STATE, get_ground_state


mp.mp.dps = FORENSIC_GROUND_STATE["dps"]

print("=" * 70)
print("FORCING CANONICAL FORENSIC GROUND STATE INTO CACHE")
print("=" * 70)
print()

print("Parameters:")
for key, value in FORENSIC_GROUND_STATE.items():
    print(f"  {key} = {value}")
print()

lambda_min, u_star, meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
    verbose=True,
)

print()
print("=" * 70)
print("CACHE OPERATION COMPLETE")
print("=" * 70)
print()
print(f"lambda_min = {mp.nstr(lambda_min, 50)}")
print()
print("Cache statistics:")
print(meta)
