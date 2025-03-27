import sympy as sp

v0x, v0y, v0z = sp.symbols("v_{0x} v_{0y} v_{0z}", real=True)
v0 = sp.Matrix([v0x, v0y, v0z])
vx, vy, vz = sp.symbols("v_x v_y v_z", real=True)
n_0 = sp.symbols("n_0", real=True, positive=True)
n_s, T_s = sp.symbols("n_s T_s", real=True, positive=True)
m = sp.symbols("m", real=True, positive=True)
v = sp.Matrix([vx, vy, vz])
# f_\text{beam}(\mathbf v,\mathbf v, t)=n_0\delta(v_x-v_0)\delta(v_y)\delta(v_z)
# f = n_0 * sp.DiracDelta(vx - v0x) * sp.DiracDelta(vy) * sp.DiracDelta(vz)
f = (
    n_s
    * (m / (2 * sp.pi * T_s)) ** (sp.Rational(3, 2))
    * sp.exp(-sp.Rational(1, 2) * m / T_s * (v - v0).dot(v - v0))
)

n = sp.integrate(f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo))
# n = n_s

n = sp.simplify(n)
print("n =", sp.latex(n), "\\\\")

# fmt: off
u = 1 / n * sp.integrate(v * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo))
# fmt: on

u = sp.simplify(u)
print("\\mathbf u =", sp.latex(u), "\\\\")

# fmt: off
xx=sp.simplify(1/n * sp.integrate(vx**2 * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)))
# xy=sp.simplify(1/n * sp.integrate(vx*vy * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)))
# xz=sp.simplify(1/n * sp.integrate(vx*vz * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)))
yy=sp.simplify(1/n * sp.integrate(vy**2 * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)))
# yz=sp.simplify(1/n * sp.integrate(vy*vz * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)))
zz=sp.simplify(1/n * sp.integrate(vz**2 * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)))
# fmt: on

# fmt: off
# vv_avg = 1 / n * sp.integrate((v * v.T) * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo))
# fmt: on
tr_vv_avg = sp.simplify(
    sp.Matrix([[xx, 0, 0], [0, yy, 0], [0, 0, zz]])
)  # significanyly faster than the above line

print("\\Tr(\\langle\\mathbf v\\mathbf v\\rangle) =", sp.latex(tr_vv_avg), "\\\\")

uu = u * u.T
uu = sp.simplify(uu)

print("\\mathbf u\\mathbf u =", sp.latex(uu), "\\\\")

p = 1 / 3 * (tr_vv_avg - sp.trace(uu)) * m * n
p = sp.simplify(p)

print("p =", sp.latex(p), "\\\\")

ked = (
    1
    / 2
    * m
    * sp.integrate(
        v.dot(v) * f, (vx, -sp.oo, sp.oo), (vy, -sp.oo, sp.oo), (vz, -sp.oo, sp.oo)
    )
)
ked = sp.simplify(ked)
print(
    "\\text{KED} = \\langle nmv^2/2\\rangle = \\frac{3}{2}nT =", sp.latex(ked), "\\\\"
)
