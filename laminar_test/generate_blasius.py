"""
Generate the Blasius self-similar solution file for the BL code.

Solves the Blasius ODE: f''' + f*f'' = 0
with f(0)=0, f'(0)=0, f'(inf)=1

Output format (text):
  n_source
  eta(1) eta(2) ... eta(n_source)
  f(1)   f(2)   ... f(n_source)
  df(1)  df(2)  ... df(n_source)
"""
import numpy as np
from scipy.integrate import solve_ivp

def blasius_ode(eta, y):
    """Blasius ODE: f''' + f*f'' = 0"""
    f, fp, fpp = y
    return [fp, fpp, -f * fpp]

# Shooting: f''(0) = 0.332057 is the known value
eta_max = 20.0
n_points = 1000
eta_span = (0, eta_max)
eta_eval = np.linspace(0, eta_max, n_points)

y0 = [0.0, 0.0, 0.46960]  # f(0), f'(0), f''(0) for f''' + ff'' = 0

sol = solve_ivp(blasius_ode, eta_span, y0, t_eval=eta_eval, rtol=1e-12, atol=1e-12)

eta = sol.t
f = sol.y[0]
df = sol.y[1]  # f' = U/U_inf

print(f"Blasius solution: eta_max={eta_max}, n={len(eta)}")
print(f"  f'(inf) = {df[-1]:.10f} (should be 1.0)")
print(f"  f''(0)  = {sol.y[2][0]:.10f} (should be 0.33206)")

# Write in the format the code expects
fname = 'blasius_solution.dat'
with open(fname, 'w') as fout:
    fout.write(f'{len(eta)}\n')
    fout.write(' '.join(f'{e:.15e}' for e in eta) + '\n')
    fout.write(' '.join(f'{v:.15e}' for v in f) + '\n')
    fout.write(' '.join(f'{v:.15e}' for v in df) + '\n')

print(f'Wrote {fname}')
