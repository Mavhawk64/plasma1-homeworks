import sympy as sp

phi, x = sp.symbols("\\phi x", real=True)
A = sp.symbols("2\\frac{e^2n_0}{\\varepsilon_0k_BT}", real=True, positive=True)

# d^2/dx^2 phi = A phi
lhs = sp.diff(phi, x, x) - A * phi

# solve the differential equation
sol = sp.dsolve(lhs, x)

print(sol)
