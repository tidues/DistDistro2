import Nregion as rg
from numericFuncs import theta

def get_R(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val=None):
    # case e == f
    if e == f:
        if i == 0:
            a = 0
        else:
            a = d
        b = (d + le)/2
        if x_val is None:
            x_val = b
        xp = theta(a, b, x_val)
        if (i, j) == (0, 0):
            r = rg.RegionBase(0, 1, lambda p: max(0, p - xp/le), lambda p: p, x_val)
        elif (i, j) == (0, 1):
            r = rg.RegionBase(0, 1, lambda p: p, lambda p: min(1, p + xp/le), x_val)
        elif (i, j) == (1, 0):
            r = rg.RegionBase(1 - (xp - d)/le, 1, lambda p: 0, lambda p: p - 1 + (xp - d)/le, x_val)
        else:
            r = rg.RegionBase(0, (xp - d)/le,  lambda p: p + 1 - (xp - d)/le, lambda p: 1, x_val)
        return r
    else:
        a = d[i, j]
        b = (le + lf + 
                min(d[0, 0] + d[1, 1], 
                    d[0, 1] + d[1, 0])) / 2
        if x_val is None:
            x_val = b
        xp = theta(a, b, x_val)
        if (i, j) == (0, 0):
            r = rg.RegionBase(0, min(p1, (xp - d[i, j])/le), lambda p: 0, lambda p: min(q1, (xp - le * p - d[i, j])/lf), x_val, bd=(True, False))
        elif (i, j) == (0, 1):
            r = rg.RegionBase(0, min(p2, (xp - d[i, j])/le), lambda p: max(q1, (- xp + le * p + d[i, j] + lf)/lf), lambda p: 1, x_val, bd=(True, False))
        elif (i, j) == (1, 0):
            r = rg.RegionBase(max(p1, (- xp + d[i, j] + le)/le), 1, lambda p: 0, lambda p: min(q2, (xp + le * p - d[i, j] - le)/lf), x_val)
        else:
            r = rg.RegionBase(max(p2, (- xp + d[i, j] + le)/le), 1, lambda p: max(q2, (- xp - le * p + d[i, j] + le + lf)/lf), lambda p: 1, x_val)
        return r


