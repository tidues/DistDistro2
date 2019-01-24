from scipy.integrate import quad, dblquad
import numericFuncs as nf
import math

# basic components of a region
class RegionBase:
    def __init__(self, pl, pu, ql, qu, xval, sgn=1, eta=False, a=0, b=math.inf, bd=(True, True)):
        self.pl = pl
        self.pu = pu
        self.ql = ql
        self.qu = qu
        self.sgn = sgn
        self.eta = eta
        self.a = a
        self.b = b
        self.xval = xval # is x if it is a region function; or upper bound if is a region
        self.bd = bd
        # self.modules = ['numpy']

    # the mass function
    def m(self, func):
        if self.ql is None:
            region = (self.pl, self.pu)
            int_func = quad
        else:
            int_func = dblquad
            region = (self.pl, self.pu, self.ql, self.qu)

        return int_func(func, *region)[0]

    # adjusted mass function
    def m_adj(self, f):
        res = self.m(f)
        if self.eta is True:
            return self.sgn * nf.eta(self.a, self.b, self.xval) * res
        else:
            return self.sgn * res

    # conditional mass function
    def m_p(self, f, p_val):
        myeta = None
        if self.bd[0] is False:
            myeta = nf.etal
        elif self.bd[1] is False:
            myeta = nf.etar
        else:
            myeta = nf.eta

        etaval = myeta(self.pl, self.pu, p_val)
        if etaval != 0:
            return quad(f, self.ql(p_val), self.qu(p_val))[0]
        else:
            return 0

    # conditional adjusted mass function
    def m_p_adj(self, f, p_val):
        res = self.m_p(f, p_val)
        if self.eta is True:
            return self.sgn * nf.eta(self.a, self.b, self.xval) * res
        else:
            return self.sgn * res

    # print region description
    def print(self):
        print('pl:\t', self.pl)
        print('pu:\t', self.pu)
        print('ql:\t', self.ql)
        print('qu:\t', self.qu)
        print('a:\t', self.a)
        print('b:\t', self.b)
        print('xval:\t', self.xval)


# a list of region bases
class Region:
    def __init__(self, bases):
        self.bases = bases
   
    # the measure of this region
    def m(self, f):
        res = 0
        for b in self.bases:
            res = res + b.m_adj(f)
        return res

    # the conditional measure of this region
    def m_p(self, f, p_val):
        res = 0
        for b in self.bases:
            res = res + b.m_p_adj(f, p_val)
        return res

    # print region description
    def print(self):
        for i, b in enumerate(self.bases):
            print('Piece ',i+1, ':\n')
            b.print()
            print('\n')

