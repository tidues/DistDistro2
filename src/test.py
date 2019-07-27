from sympy.abc import p, q
from sympy import Piecewise, integrate, simplify, Max
from sympy.plotting import plot3d, plot

#phi_pq = Piecewise(
#        (180 * p * q * (p + q - 1) ** 2, ((p >= 0) & (q >= 0) & (p + q <= 1))),
#        (180 * (p - 1) * (q - 1) * (p + q - 1) ** 2, ((p <= 1) & (q <= 1) & (p + q >= 1))),
#        (0, True)
#        )
phi_pq = Piecewise(
        (180 * p * q * (p + q - 1) ** 2, (p + q <= 1)),
        (180 * (p - 1) * (q - 1) * (p + q - 1) ** 2, (p + q >= 1))
        )

plot3d(phi_pq, (p, 0, 1), (q, 0 , 1))

intp = simplify(integrate(phi_pq, (p, 0, 1)))
intq = simplify(integrate(phi_pq, (q, 0, 1)))
print(intp)
print(intq)

plot(intp, (q, 0, 1))
plot(intq, (p, 0, 1))


