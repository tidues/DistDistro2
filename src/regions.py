import region as rg
from sympy import *
from sympy.abc import x, p, q
from commonFuncs import theta

def get_R(e, f, g, le, lf, duv, d, p1, p2, q1, q2):
    # case e == f
    if e == f:
        # case ee00
        for (i, j) in g.two2:
            a = g.a[e, f, i, j]
            b = g.b[e, f, i, j]
            xp = theta(a, b, x)
            if (i, j) == (0, 0):
                r1 = rg.RegionBase(0, xp/le, 0, p)
                r2 = rg.RegionBase(xp/le, 1, p - xp/le, p)
                R = rg.Region((r1,r2))
            elif (i, j) == (0, 1):
                r1 = rg.RegionBase(0, 1 - xp/le, p, p + xp/le)
                r2 = rg.RegionBase(1 - xp/le, 1, p, 1)
                R = rg.Region((r1,r2))
            elif (i, j) == (1, 0):
                r1 = rg.RegionBase(1 - (xp - duv)/le, 1, 0, p - 1 + (xp - duv)/le)
                R = rg.Region((r1,))
            else:
                r1 = rg.RegionBase(0, (xp - duv)/le,  p + 1 - (xp - duv)/le, 1)
                R = rg.Region((r1,))
            g.Rx[e, f, i, j] = R
            g.R[e, f, i, j] = R.R_const(g.b[e, f, i, j])
    else:
        for (i, j) in g.two2:
            a = g.a[e, f, i, j]
            b = g.b[e, f, i, j]
            xp = theta(a, b, x)
            if (i, j) == (0, 0):
                r1 = rg.RegionBase(0, (xp - d[i, j])/le, 0, (xp - le * p - d[i, j])/lf)
                r2 = rg.RegionBase(p1, (xp - d[i, j])/le, 0, (xp - le * p - d[i, j])/lf, sgn=-1, eta=True, a=le * p1 + d[i, j])
                r3 = rg.RegionBase(0, (xp - lf * q1 - d[i, j])/le, q1, (xp - le * p - d[i, j])/lf, sgn=-1, eta=True, a=lf * q1 + d[i, j])
                R = rg.Region((r1, r2, r3))
            elif (i, j) == (0, 1):
                r1 = rg.RegionBase(0, (xp - d[i, j])/le, (- xp + le * p + d[i, j] + lf)/lf, 1)
                r2 = rg.RegionBase(p2, (xp - d[i, j])/le, (- xp + le * p + d[i, j] + lf)/lf, 1, sgn=-1, eta=True, a=le * p2 + d[i, j])
                r3 = rg.RegionBase(0, (xp + lf * q1 - d[i, j] - lf)/le, (- xp + le * p + d[i, j] + lf)/lf, q1, sgn=-1, eta=True, a=(1 - q1) * lf + d[i,j])
                R = rg.Region((r1, r2, r3))
            elif (i, j) == (1, 0):
                r1 = rg.RegionBase((- xp + d[i, j] + le)/le, 1, 0, (xp + le * p - d[i, j] - le)/lf)
                r2 = rg.RegionBase((- xp + d[i, j] + le)/le, p1, 0, (xp + le * p - d[i, j] - le)/lf, sgn=-1, eta=True, a=(1 - p1) * le + d[i,j])
                r3 = rg.RegionBase((- xp + lf * q2 + d[i, j] + le)/le, 1, q2, (xp + le * p - d[i, j] - le)/lf, sgn=-1, eta=True, a=lf * q2 + d[i, j])
                R = rg.Region((r1, r2, r3))
            else:
                r1 = rg.RegionBase((- xp + d[i, j] + le)/le, 1, (- xp - le * p + d[i, j] + le + lf)/lf, 1)
                r2 = rg.RegionBase((- xp + d[i, j] + le)/le, p2, (- xp - le * p + d[i, j] + le + lf)/lf, 1, sgn=-1, eta=True, a=(1 - p2) * le + d[i,j])
                r3 = rg.RegionBase((- xp - lf * q2 + d[i, j] + le + lf)/le, 1, (- xp - le * p + d[i, j] + le + lf)/lf, q2, sgn=-1, eta=True, a=(1 - q2) * lf + d[i,j])
                R = rg.Region((r1, r2, r3))
            g.Rx[e, f, i, j] = R
            g.R[e, f, i, j] = R.R_const(g.b[e, f, i, j])


# get L
def get_L(e, f, g, le, lf, duv, d, p1, p2, q1, q2):
    if e == f:
        for (i, j) in g.two2:
            a = g.a[e, f, i, j]
            b = g.b[e, f, i, j]
            xp = theta(a, b, x)
            if (i, j) == (0, 0):
                r = rg.RegionBase(xp/le, 1, None, None)
                R = rg.Region((r,))
            elif (i, j) == (0, 1):
                r = rg.RegionBase(0, 1 - xp/le, None, None)
                R = rg.Region((r,))
            elif (i, j) == (1, 0):
                r = rg.RegionBase(1 - (xp - duv)/le, 1, None, None)
                R = rg.Region((r,))
            elif (i, j) == (1, 1):
                r = rg.RegionBase(0, (xp - duv)/le, None, None)
                R = rg.Region((r,))
            g.L[e, f, i, j] = R
    else:
        for (i, j) in g.two2:
            a = g.a[e, f, i, j]
            b = g.b[e, f, i, j]
            xp = theta(a, b, x)
            
            if (i, j) == (0, 0):
                r = rg.RegionBase(
                        Max(0,(2 * xp - d[0, 0] - d[0, 1] - lf)/(2 * le)), 
                        Min(p1, (xp - d[i, j])/le), 
                        None, None)
                R = rg.Region((r,))
            elif (i, j) == (0, 1):
                r = rg.RegionBase(
                        Max(0, (2 * xp - d[0, 0] - d[0, 1] - lf)/(2 * le)), 
                        Min(p2, (xp - d[i, j])/le), 
                        None, None)
                R = rg.Region((r,))
            elif (i, j) == (1, 0):
                r = rg.RegionBase(
                        Max(p1, (- xp + d[i, j] + le)/le), 
                        Min(1, (- 2 * xp + d[1, 0] + d[1, 1] + 2* le + lf)/(2 * le)), 
                        None, None)
                R = rg.Region((r,))
            elif (i, j) == (1, 1):
                r = rg.RegionBase(
                        Max(p2, (- xp + d[i, j] + le)/le), 
                        Min(1, (- 2 * xp + d[1, 0] + d[1, 1] + 2* le + lf)/(2 * le))
                        , None, None)
                R = rg.Region((r,))

            g.L[e, f, i, j] = R

