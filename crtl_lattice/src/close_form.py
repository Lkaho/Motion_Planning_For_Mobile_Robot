# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

from sympy import *
init_printing()
px0, py0, pz0, vx0, vy0, vz0, pxf, pyf, pzf, vxf, vyf, vzf, T= symbols('px0 py0 pz0 vx0 vy0 vz0 pxf pyf pzf vxf vyf vzf T')

delta_px = pxf - vx0 * T - px0
delta_py = pyf - vy0 * T - py0
delta_pz = pzf - vz0 * T - pz0
delta_vx = vxf - vx0
delta_vy = vyf - vy0
delta_vz = vzf - vz0

b = Matrix([delta_px, delta_py, delta_pz, delta_vx, delta_vy, delta_vz])

A = Matrix([
    [1/6 * T ** 3, 0, 0, 1/2 * T** 2, 0, 0],
    [0,  1/6 * T ** 3, 0, 0, 1/2 * T ** 2,  0],
    [0, 0, 1/6 * T ** 3, 0, 0, 1/2 * T ** 2],
    [1/2  * T ** 2, 0, 0, T, 0, 0],
    [0, 1/2 * T ** 2, 0, 0, T, 0],
    [0, 0, 1/2 * T ** 2, 0, 0, T]
])

x = simplify(A ** (-1) * b)

M = Matrix(
    [
        [(T**3)/3,        0,        0, (T**2)/2,        0,        0],
        [       0, (T**3)/3,        0,        0, (T**2)/2,        0],
        [       0,        0, (T**3)/3,        0,        0, (T**2)/2],
        [(T**2)/2,        0,        0,        T,        0,        0],
        [       0, (T**2)/2,        0,        0,        T,        0],
        [       0,        0, (T**2)/2,        0,        0,        T],
    ]
)

J = collect(
    expand(
        simplify(
            Transpose(x) * M * x
        )[0] + T
    ),
    syms=T
)

print(J)

dotJ = collect(
    expand(
        simplify(diff(J, T))
    ),
    syms = T
)
print(dotJ)

dotJ_new = collect(
    expand(
        simplify(dotJ * T**4)
    ),
    syms = T
)






