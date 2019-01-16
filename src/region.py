from sympy import *
import commonFuncs as cf
from sympy.abc import p, q, x


# basic components of a region
class RegionBase:
    # (pl, pu, ql, qu) tuple representation,
    # sgn: remove or add;
    # eta: a gate function to control the support
    # a, b: the eta function bounds
    # bd: indicate whether the corresponding boundary 
    #     (just pl, pu for now) is included
    #       in the region
    def __init__(self, pl, pu, ql, qu, sgn=1, eta=False, a=0, b=oo, xval=x, bd=(True, True)):
        zero = numbers.Zero()
        self.pl = cf.toSym(pl)
        self.pu = cf.toSym(pu)
        if ql is not None:
            ql = cf.toSym(ql)
        if qu is not None:
            qu = cf.toSym(qu)
        self.ql = ql
        self.qu = qu
        self.sgn = sgn
        self.eta = eta
        self.a = a
        self.b = b
        self.xval = xval
        self.bd = bd
        # self.modules = ['numpy']

    # the mass function
    def m(self, f):
        if self.ql is None:
            res = integrate(f, (p, self.pl, self.pu))
        else:
            res = integrate(f, (q, self.ql, self.qu), (p, self.pl, self.pu))
        return res

    # adjusted mass function
    def m_adj(self, f):
        res = self.m(f)
        if self.eta is True:
            return self.sgn * cf.eta(self.a, self.b, self.xval) * res
        else:
            return self.sgn * res

    # conditional mass function
    def m_p(self, f):
        myeta = None
        if self.bd[0] is False:
            myeta = cf.etal
        elif self.bd[1] is False:
            myeta = cf.etar
        else:
            myeta = cf.eta

        return myeta(self.pl, self.pu, p) * integrate(f, (q, self.ql, self.qu))

    # conditional adjusted mass function
    def m_p_adj(self, f):
        res = self.m_p(f)
        if self.eta is True:
            return self.sgn * cf.eta(self.a, self.b, self.xval) * res
        else:
            return self.sgn * res





        

    ## eval m_adj
    #def m_adj_eval(self, f):
    #    print(f)
    #    res = self.m_adj(f)
    #    print(res)
    #    print(self.xval)
    #    return self.eval(res, self.xval)

    # produce fixed region bases
    def RB_const(self, val):
        return RegionBase(
                self.eval(self.pl, val),
                self.eval(self.pu, val),
                self.eval(self.ql, val),
                self.eval(self.qu, val),
                sgn = self.sgn,
                eta = self.eta,
                a = self.a,
                b = self.b,
                xval = val
                )

    # def lambdify
    def eval(self, expr, x_val):
        return expr.subs(x, x_val)
        # res = lambdify(x, expr)(x_val)
        # print(res)
        # if isnumber(res) and self.rat:
        #     res =  Rational(res)
        # return res


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
    def m_p(self, f):
        res = 0
        for b in self.bases:
            res = res + b.m_p_adj(f)
        return res

    # print region description
    def print(self):
        for i, b in enumerate(self.bases):
            print('Piece ',i+1, ':\n')
            b.print()
            print('\n')

    # generate constant Regions 
    def R_const(self, val):
        rs = []
        for b in self.bases:
            rs.append(b.RB_const(val))
        return Region(rs)

# check whether the input is a number or not
def isnumber(x):
    tstr = str(type(x))
    return 'numbers' in tstr or 'int' in tstr or 'float' in tstr

if __name__ == '__main__':
    le = 3
    mi = Function('min')
    r1 = RegionBase(0, x/le, 0, p)
    r2 = RegionBase(x/le, mi(1,3), p-x/le, p, eta=True, a=1, b=3)

    f = 2 * p

#    print(r1.m(f))
#    print(r2.m(f))
    
    R = Region((r1,r2))
    R.assign_f(f)
    print(R.m)
    R.m_lambdify()
    print(R.m_func(0))
    print(R.m_func(0.9))
    print(R.m_func(1))
    print(R.m_func(2))
    print(R.m_func(3))
    print(R.m_func(3.1))
    print(R.m_func(4))

