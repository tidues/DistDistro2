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
            r1 = rg.RegionBase(0, xp/le, lambda p: 0, lambda p: p, x_val, bd=(True, False))
            r2 = rg.RegionBase(xp/le, 1, lambda p: p - xp/le, lambda p: p, x_val)
            R = rg.Region((r1,r2))
        elif (i, j) == (0, 1):
            r1 = rg.RegionBase(0, 1 - xp/le, lambda p: p, lambda p: p + xp/le, x_val, bd=(True, False))
            r2 = rg.RegionBase(1 - xp/le, 1, lambda p: p, lambda p: 1, x_val)
            R = rg.Region((r1,r2))
        elif (i, j) == (1, 0):
            r1 = rg.RegionBase(1 - (xp - d)/le, 1, lambda p: 0, lambda p: p - 1 + (xp - d)/le, x_val)
            R = rg.Region((r1,))
        else:
            r1 = rg.RegionBase(0, (xp - d)/le,  lambda p: p + 1 - (xp - d)/le, lambda p: 1, x_val)
            R = rg.Region((r1,))
        return R
    else:
        a = d[i, j]
        b = (le + lf + 
                min(d[0, 0] + d[1, 1], 
                    d[0, 1] + d[1, 0])) / 2
        if x_val is None:
            x_val = b
        xp = theta(a, b, x_val)
        if (i, j) == (0, 0):
            r1 = rg.RegionBase(0, (xp - d[i, j])/le, lambda p: 0, lambda p: (xp - le * p - d[i, j])/lf, x_val)
            r2 = rg.RegionBase(p1, (xp - d[i, j])/le, lambda p: 0, lambda p: (xp - le * p - d[i, j])/lf, x_val, sgn=-1, eta=True, a=le * p1 + d[i, j])
            r3 = rg.RegionBase(0, (xp - lf * q1 - d[i, j])/le, lambda p: q1, lambda p: (xp - le * p - d[i, j])/lf, x_val, sgn=-1, eta=True, a=lf * q1 + d[i, j])
            R = rg.Region((r1, r2, r3))
        elif (i, j) == (0, 1):
            r1 = rg.RegionBase(0, (xp - d[i, j])/le, lambda p: (- xp + le * p + d[i, j] + lf)/lf, lambda p: 1, x_val)
            r2 = rg.RegionBase(p2, (xp - d[i, j])/le, lambda p: (- xp + le * p + d[i, j] + lf)/lf, lambda p: 1, x_val, sgn=-1, eta=True, a=le * p2 + d[i, j])
            r3 = rg.RegionBase(0, (xp + lf * q1 - d[i, j] - lf)/le, lambda p: (- xp + le * p + d[i, j] + lf)/lf, lambda p: q1, x_val, sgn=-1, eta=True, a=(1 - q1) * lf + d[i,j])
            R = rg.Region((r1, r2, r3))
        elif (i, j) == (1, 0):
            r1 = rg.RegionBase((- xp + d[i, j] + le)/le, 1, lambda p: 0, lambda p: (xp + le * p - d[i, j] - le)/lf, x_val)
            r2 = rg.RegionBase((- xp + d[i, j] + le)/le, p1, lambda p: 0, lambda p: (xp + le * p - d[i, j] - le)/lf, x_val, sgn=-1, eta=True, a=(1 - p1) * le + d[i,j], bd=(True, False))
            r3 = rg.RegionBase((- xp + lf * q2 + d[i, j] + le)/le, 1, lambda p: q2, lambda p: (xp + le * p - d[i, j] - le)/lf, x_val, sgn=-1, eta=True, a=lf * q2 + d[i, j])
            R = rg.Region((r1, r2, r3))
        else:
            r1 = rg.RegionBase((- xp + d[i, j] + le)/le, 1, lambda p: (- xp - le * p + d[i, j] + le + lf)/lf, lambda p: 1, x_val)
            r2 = rg.RegionBase((- xp + d[i, j] + le)/le, p2, lambda p: (- xp - le * p + d[i, j] + le + lf)/lf, lambda p: 1, x_val, sgn=-1, eta=True, a=(1 - p2) * le + d[i,j], bd=(True, False))
            r3 = rg.RegionBase((- xp - lf * q2 + d[i, j] + le + lf)/le, 1, lambda p: (- xp - le * p + d[i, j] + le + lf)/lf, lambda p: q2, x_val, sgn=-1, eta=True, a=(1 - q2) * lf + d[i,j])
            R = rg.Region((r1, r2, r3))
        return R

## get L
#def get_L(e, f, g, le, lf, duv, d, p1, p2, q1, q2):
#    if e == f:
#        for (i, j) in g.two2:
#            a = g.a[e, f, i, j]
#            b = g.b[e, f, i, j]
#            xp = theta(a, b, x)
#            if (i, j) == (0, 0):
#                r = rg.RegionBase(xp/le, 1, None, None)
#                R = rg.Region((r,))
#            elif (i, j) == (0, 1):
#                r = rg.RegionBase(0, 1 - xp/le, None, None)
#                R = rg.Region((r,))
#            elif (i, j) == (1, 0):
#                r = rg.RegionBase(1 - (xp - duv)/le, 1, None, None)
#                R = rg.Region((r,))
#            elif (i, j) == (1, 1):
#                r = rg.RegionBase(0, (xp - duv)/le, None, None)
#                R = rg.Region((r,))
#            g.L[e, f, i, j] = R
#    else:
#        for (i, j) in g.two2:
#            a = g.a[e, f, i, j]
#            b = g.b[e, f, i, j]
#            xp = theta(a, b, x)
#            
#            if (i, j) == (0, 0):
#                r = rg.RegionBase(
#                        mmax(0,(2 * xp - d[0, 0] - d[0, 1] - lf)/(2 * le)), 
#                        mmin(p1, (xp - d[i, j])/le), 
#                        #Max(0,(2 * xp - d[0, 0] - d[0, 1] - lf)/(2 * le)), 
#                        #Min(p1, (xp - d[i, j])/le), 
#                        None, None, bd=(True, False))
#                R = rg.Region((r,))
#            elif (i, j) == (0, 1):
#                r = rg.RegionBase(
#                        mmax(0, (2 * xp - d[0, 0] - d[0, 1] - lf)/(2 * le)), 
#                        mmin(p2, (xp - d[i, j])/le), 
#                        #Max(0, (2 * xp - d[0, 0] - d[0, 1] - lf)/(2 * le)), 
#                        #Min(p2, (xp - d[i, j])/le), 
#                        None, None, bd=(True, False))
#                R = rg.Region((r,))
#            elif (i, j) == (1, 0):
#                r = rg.RegionBase(
#                        mmax(p1, (- xp + d[i, j] + le)/le), 
#                        mmin(1, (- 2 * xp + d[1, 0] + d[1, 1] + 2* le + lf)/(2 * le)), 
#                        #Max(p1, (- xp + d[i, j] + le)/le), 
#                        #Min(1, (- 2 * xp + d[1, 0] + d[1, 1] + 2* le + lf)/(2 * le)), 
#                        None, None)
#                R = rg.Region((r,))
#            elif (i, j) == (1, 1):
#                r = rg.RegionBase(
#                        mmax(p2, (- xp + d[i, j] + le)/le), 
#                        mmin(1, (- 2 * xp + d[1, 0] + d[1, 1] + 2* le + lf)/(2 * le)),
#                        #Max(p2, (- xp + d[i, j] + le)/le), 
#                        #Min(1, (- 2 * xp + d[1, 0] + d[1, 1] + 2* le + lf)/(2 * le)),
#                        None, None)
#                R = rg.Region((r,))
#
#            g.L[e, f, i, j] = R

