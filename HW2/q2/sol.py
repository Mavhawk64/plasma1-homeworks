import sympy as sp
from sympy import Function, Matrix, symbols

# Given values
E = Matrix([10, 0, 0])  # Electric field in V/m
B = Matrix([0, 0, 1e-4])  # Magnetic field in T
q = 1.602e-19  # Proton charge in C
m = 1.672e-27  # Proton mass in kg
v0 = Matrix([2.414e5, -1.414e5, 0])  # Initial velocity in m/s
x0 = Matrix([0, 0, 0])  # Initial position in m

# Part (a) Cyclotron frequency
Omega_c = q * B.norm() / m

# Part (b) E x B drift velocity
drift_velocity = E.cross(B) / B.norm() ** 2

# Part (c) Solve x(t) and v(t)
# done analytically
# Part (d) Gyration speed
v_parallel = drift_velocity
v_perp = (v0 - v_parallel).norm()

# Part (e) Gyration energy in eV
energy_J = 0.5 * m * v_perp**2
energy_eV = energy_J / q

# Part (f) Magnetic moment
mu = (m * v_perp**2) / (2 * B.norm())

# Output results
results = {
    "Cyclotron frequency (rad/s)": Omega_c,
    "E x B drift velocity (m/s)": drift_velocity,
    "Gyration speed (m/s)": v_perp,
    "Gyration energy (eV)": energy_eV,
    "Magnetic moment (A*m^2)": mu,
}

for key, value in results.items():
    print(f"{key}: {value.evalf()}")
