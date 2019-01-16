from enums import Stats
import commonFuncs as cf
from sympy import *
from sympy.abc import p, q, x
from sympy.plotting import plot
from plot import plot1d

class Formula:
    def __init__(self, g):
        self.stat = None  # the formula type
        self.g = g
        self.fs = {}  # all formulas
        self.mods = ['numpy', {'theta': cf.theta, 'eta': cf.eta}]  ## modules for lambdify
        self.idx_num = 0  # dimension of indexes for self.fs
        self.val_keys = None  # the symbols for values
        self.plot_info = (None, x, -1, self.g.d_max + 1)  # plotting parameters
        self.unikey = 'zero'

    # evaluation
    def eval(self, *params):
        # split params into keys and vals
        keys, subval = self.__split_input(*params)

        # generate formula if nece
        self.__gen_formula_wrapper(keys)

        # select formula then return the evaluation
        return self.__subs(self.fs[self.mkkey(keys)], subval)

    # output the formula
    def formula(self, *params):
        # split params into keys and vals
        keys, subval = self.__split_input(*params)
        # generate formula if nece
        self.__gen_formula_wrapper(keys)
        return self.__subs(self.fs[self.mkkey(keys)], subval)

    # input f_info = (f, var, lb, ub)
    def plot(self, *params, method=0, step=0.1):
        if self.stat != Stats.MOMENT:
            if len(params) == 0:
                self.gen_formula()
                f = self.fs[self.unikey]
            else:
                # split params into keys and vals
                keys, subval = self.__split_input(*params)
                self.__gen_formula_wrapper(keys)
                f = self.__subs(self.fs[self.mkkey(keys)], subval)

            var = self.plot_info[1]
            lb = self.plot_info[2]
            ub = self.plot_info[3]

            if method == 0:
                try:
                    plot(f, (var, lb, ub))
                except:
                    self.plot(method=1, step=step)
            else:
                f_lambda = lambdify(var, f, modules=self.mods)
                plot1d(f_lambda, lb, ub, step)
        else:
            print('moment function cannot be plotted.')

    def __gen_formula_wrapper(self, keys):
        if keys == self.unikey:
            self.gen_formula()
        else:
            self.gen_formula(*keys)

    def mkkey(self, keys):
        if len(keys) == 1:
            return keys[0]
        else:
            return keys

    def __split_input(self, *params):
        if self.idx_num == 0:
            keys = self.unikey
            vals = params
        else:
            keys = []
            vals = []
            for i, v in enumerate(params):
                if i < self.idx_num:
                    keys.append(v)
                else:
                    vals.append(v)
            keys = tuple(keys)
        
        # make subs info
        if len(vals) == 0:
            subval = None
        elif len(vals) == 1:
            subval = [(self.val_keys[0], vals[0])]
        else:
            subval = list(zip(self.val_keys, vals))

        return (keys, subval)

    # return formula with subs
    def __subs(self, formula, subval):
        if subval is None:
            return formula
        else:
            return formula.subs(subval)

    # abstract method
    def gen_formula(self):
        # realized by the inherited classes
        pass
    
class Moment(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.MOMENT
        self.idx_num = 1

    def gen_formula(self, k):
        if k in self.fs:
            return False

        print('calculating the ', k, 'th  moment...')
        res = 0
        alphas = cf.A_curl(2, k)
        g = self.g
        for e in g.edges():
            for f in g.edges():
                for (i, j) in g.two2:
                    ip = i + 1
                    jp = j + 1
                    for alpha in alphas:
                        # produce entry info
                        if (e, f, i, j, alpha) not in g.moment_info:
                            c0 = (g.gx[e] * g.gy[f] * g.l[f] ** alpha[0] * g.l[e] ** alpha[1])
                            if e == f:
                                w = (ip + jp -1) * alpha[0] + (ip + jp) * alpha[1]
                                c1 = (-1) ** w
                            else:
                                w = j *  alpha[0] + i * alpha[1]
                                c1 = (-1) ** w

                            c = c0 * c1
                            # q = Symbol('q')
                            # p = Symbol('p')
                            expr = q ** alpha[0] * p ** alpha[1] * g.phi_pq
                            m = g.R[e, f, i, j].m(expr)
                            g.moment_info[e, f, i, j, alpha] = c * m

                        # sum up
                        if e == f:
                            c2 = cf.ncrs(k, alpha) * (i * (g.d[e[0]][e[1]] + g.l[e])) ** (k - alpha[0] - alpha[1])
                        else:
                            c2 = cf.ncrs(k, alpha) * (g.d[e[i]][f[j]] + i * g.l[e] + j * g.l[f]) ** (k - alpha[0] - alpha[1])
                        
                        res += c2 * self.g.moment_info[e, f, i, j, alpha]
        self.fs[k] = res
        return True

#    def eval(self, k):
#        self.__gen_formula(k)
#        return self.fs[k]

    
class CDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CDF
        self.val_keys = (x,)

    def gen_formula(self):
        if self.unikey in self.fs:
            return False

        print('generating the cdf...')

        g = self.g
        # x = Symbol('x')
        expr = 0

        for e in g.edges():
            for f in g.edges():
                gxy = g.gx[e] * g.gy[f]
                tmp = 0
                for (i, j) in g.two2:
                    tmp += g.Rx[e, f, i, j].m(g.phi_pq)

                expr += gxy * tmp

        self.fs[self.unikey] = expand(expr)
        return True

#    def eval(self, x_val):
#        # x = Symbol('x')
#        return self.f.subs(x, x_val)


class PDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.PDF
        self.val_keys = (x,)

    def gen_formula(self):
        if self.unikey in self.fs:
            return False

        print('generating the pdf...')
        g = self.g
        expr = 0

        # mat_D = {}
        for e in g.edges():
            for f in g.edges():
                const = g.gx[e] * g.gy[f] / g.l[f]
                tmp = 0
                # mxs = {}
                for (i, j) in g.two2:
                    c = cf.eta(g.a[e, f, i, j], g.b[e, f, i, j], x)
                    phis = g.phi_pq.subs(q, g.q[e, f, i, j])
                    m = g.L[e, f, i, j].m(phis)
                    tmp += c * m

                    # print((e, f, i, j))
                    # print('c:\t', c)
                    # print('phis:\t', phis)
                    # print('m:\t', m)
                    # print('L:\t')
                    # g.L[e, f, i, j].print()
                    # mxs[i, j] = (c * m / g.l[f]).subs(x, 3)
                    
                # if e == ('1','2') and f == ('2', '3'):
                #     print('mxs2:', mxs)

                    

                # mat_D[e, f] = (tmp / g.l[f])

                expr += const * tmp

        # self.mat_D = mat_D
        self.fs[self.unikey] = expand(expr)
        #self.f = lambdify(x, self.expr, modules=self.modules)
        return True

#    def eval(self, x_val):
#        # x = Symbol('x')
#        return self.f.subs(x, x_val)


###############################################################
# Conditional Formulas
###############################################################
class CMoment(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CMOMENT
        self.idx_num = 2
        self.val_keys = (p,)
        self.plot_info = (None, p, 0, 1)

    def gen_formula(self, k, e):
        if (k, e) in self.fs:
            return False

        print('generating the conditional ', k, 'th moment function...')
        res = 0
        alphas = cf.A_curl(2, k)
        g = self.g
        for f in g.edges():
            for (i, j) in g.two2:
                ip = i + 1
                jp = j + 1
                for alpha in alphas:
                    # produce entry info
                    c0 = (g.gy[f] * g.l[f] ** alpha[0] * g.l[e] ** alpha[1])
                    if e == f:
                        w = (ip + jp -1) * alpha[0] + (ip + jp) * alpha[1]
                        c1 = (-1) ** w
                    else:
                        w = j *  alpha[0] + i * alpha[1]
                        c1 = (-1) ** w

                    c = c0 * c1
                    expr = q ** alpha[0] * (g.phi_pq / g.phi_p)

                    m = g.R[e, f, i, j].m_p(expr)
                    tmpres = c * m * p ** alpha[1]

                    # sum up
                    if e == f:
                        c2 = cf.ncrs(k, alpha) * (i * (g.d[e[0]][e[1]] + g.l[e])) ** (k - alpha[0] - alpha[1])
                    else:
                        c2 = cf.ncrs(k, alpha) * (g.d[e[i]][f[j]] + i * g.l[e] + j * g.l[f]) ** (k - alpha[0] - alpha[1])
                    
                    res += c2 * tmpres

        self.fs[k, e] = res
        return True

#    def eval(self, k, e, p_val):
#        if (k, e) not in self.fs:
#            self.fs[k, e] = self.__gen_formula(k, e)
#            
#        return self.fs[k, e].subs(p, p_val)
        
#    def formula(self, k, e):
#        print(self.fs[k, e])

#    def plot(self, k, e):
#        info = (self.fs[k, e], p, 0, 1)
#        super().plot(f_info=info)


class CCDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CCDF
        self.idx_num = 1
        self.val_keys = (p, x)

    def gen_formula(self, e):
        if e in self.fs:
            return False

        print('generating the conditional cdf...')

        g = self.g
        # x = Symbol('x')
        expr = 0

        for f in g.edges():
            tmp = 0
            for (i, j) in g.two2:
                tmp += g.Rx[e, f, i, j].m_p(g.phi_pq / g.phi_p)
                mtmp = g.Rx[e, f, i, j].m_p(g.phi_pq / g.phi_p)
                print((f, i, j), '\t', mtmp.subs([(p, 0), (x, 3)]))


            expr += g.gy[f]* tmp
        print('total', expr.subs([(p, 0), (x, 3)]))



        self.fs[e] = expand(expr)
        return True

#    def eval(self, e, p_val, x_val):
#        self.__gen_formula(e)
#        return self.fs[e].subs([(p, p_val), (x, x_val)])

#    # formula in x and p for given e
#    def formula(self, e):
#        self.__gen_formula(e)
#        return self.fs[e]

#    # pdf for given e, p
#    def formula_p(self, e, p_val):
#        self.__gen_formula(e)
#        return self.fs[e].subs(p, p_val)

#    # plot the pdf for given p
#    def plot_p(self, e, p_val):
#        self.__gen_formula(e)
#        info = (self.fs[e].subs(p, p_val), None, None, None)
#        super().plot(f_info=info)


class CPDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CPDF
        self.idx_num = 1
        self.val_keys = (p, x)

    def gen_formula(self, e):
        if e in self.fs:
            return False

        print('generating the conditional pdf...')

        g = self.g
        expr = 0

        for f in g.edges():
            const = g.gy[f] / g.l[f]
            tmp = 0
            # mxs = {}
            for (i, j) in g.two2:
                c0 = cf.eta(g.a[e, f, i, j], g.b[e, f, i, j], x)
                c1 = cf.eta(g.L[e, f, i, j].bases[0].pl, g.L[e, f, i, j].bases[0].pu, p)
                phis = (g.phi_pq / g.phi_p).subs(q, g.q[e, f, i, j])

            expr += const * tmp

        self.fs[e] = expand(expr)
        return True

#    def eval(self, e, p_val, x_val):
#        if e not in self.fs:
#            self.fs[e] = self.__gen_formula(e)
#
#        return self.fs[e].subs([(p, p_val), (x, x_val)])

#    # formula in x and p for given e
#    def formula(self, e):
#        return self.fs[e]

#    # pdf for given e, p
#    def formula_p(self, e, p_val):
#        return self.fs[e].subs(p, p_val)

#    # plot the pdf for given p
#    def plot_p(self, e, p_val):
#        info = (self.fs[e].subs(p, p_val), None, None, None)
#        super().plot(f_info=info)
