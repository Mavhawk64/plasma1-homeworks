import sympy as sp

Q, e0, r = sp.symbols("Q \\varepsilon_0 r", real=True)
l = sp.symbols("\\lambda_D", positive=True)

phi = Q / (4 * sp.pi * e0 * r) * sp.exp(-sp.sqrt(2) * r / l)

# 1/r^2 d/dr (r^2 d/dr phi)
lhs = -4 * sp.pi * e0 * sp.diff(r**2 * sp.diff(phi, r), r)

lhs = sp.simplify(lhs)

# integrate wrt r (0->R)
R = sp.symbols("R", positive=True)

Q_p = sp.integrate(lhs, (r, 0, R))

Q_p = sp.simplify(Q_p)

print(sp.latex(Q_p))

print("\n\n")

# subs R = l, 2l, 10l and print each result:
for R_val in [l, 2 * l, 10 * l]:
    print(
        "Q_p(" + str(sp.latex(R_val)) + ") =",
        sp.latex(sp.simplify(Q_p.subs(R, R_val))),
        "\\\\",
    )

print("\n\n")
# same thing but output as a decimal
for R_val in [l, 2 * l, 10 * l]:
    print(
        "Q_p(" + str(R_val) + ") =",
        sp.N(sp.simplify(Q_p.subs(R, R_val))),
    )
