from .enums import Stats
from . import commonFuncs as cf
from . import numericFuncs as nf
from sympy import *
from sympy.abc import p, q, x
from .plot import plot1d
from .serializedFormula import *
import os
import dill
from .pyprelude import EasyWriter as ew


class Formula:
    def __init__(self, g):
        self.stat = None  # the formula type
        self.g = g
        self.fs = {}  # all formulas
        self.N_fs = {}
        self.mods = ['numpy', {'theta': nf.theta, 'eta': nf.eta, 'etal': nf.etal, 'etar': nf.etar, 'mmin': min, 'mmax': max}]  ## modules for lambdify
        self.idx_num = 0  # dimension of indexes for self.fs
        self.val_keys = None  # the symbols for values
        self.plot_info = (None, x, -1, self.g.d_max + 1)  # plotting parameters
        self.unikey = 'zero'
        self.folder = './.formulas/'
        self.resfolder = './results/'
        self.fullsave = None

    # evaluation
    def eval(self,  *params, save=True):
        # split params into keys and vals
        keys, subval, vals = self.__split_input(*params)

        # generate formula if nece
        self.__gen_formula_wrapper(keys)

        myf = self.get_N_fs(self.mkkey(keys))

        # select formula then return the evaluation
        res = self.__apply(myf, vals)
        if save:
            g = self.g
            fname = self.resfolder + g.name + '_' + g.phi.name + '_' + str(self.stat) + '.dat'
            val = str(params) + '\t' + str(res) + '\n'
            ew.wFile(fname, val)
            
        return res
        # return self.__subs(self.fs[self.mkkey(keys)], subval)

    # lambdify
    def get_N_fs(self, keys):

        if keys not in self.N_fs:
            expr = self.fs[keys]
            if self.val_keys is None:
                self.N_fs[keys] = expr
            elif len(self.val_keys) == 1:
                self.N_fs[keys] = lambdify(self.val_keys[0], expr, modules=self.mods)
            else:
                self.N_fs[keys] = lambdify(self.val_keys, expr, modules=self.mods)
        return self.N_fs[keys]

    # output the formula
    def formula(self, *params):
        # split params into keys and vals
        keys, subval, vals = self.__split_input(*params)
        # generate formula if nece
        self.__gen_formula_wrapper(keys)
        return self.__subs(self.fs[self.mkkey(keys)], subval)

    # input f_info = (f, var, lb, ub)
    def plot(self, *params, step=0.01, save=True, show=False):
        if self.stat != Stats.MOMENT:
            if len(params) == 0:
                self.gen_formula()
                f = self.fs[self.unikey]
            else:
                # split params into keys and vals
                keys, subval, vals = self.__split_input(*params)
                self.__gen_formula_wrapper(keys)
                f = self.__subs(self.fs[self.mkkey(keys)], subval)

            var = self.plot_info[1]
            lb = self.plot_info[2]
            ub = self.plot_info[3]

            f_lambda = lambdify(var, f, modules=self.mods)

            if save:
                g = self.g
                figname = self.resfolder + g.name + '_' + g.phi.name + '_' + str(self.stat) + '_' + str(params) + '.png'
                plot1d(f_lambda, lb, ub, step, svname=figname, show=show)
            else:
                plot1d(f_lambda, lb, ub, step, show=show)

        else:
            print('moment function cannot be plotted.')

    # save formula
    def save_formulas(self):
        sf = SFormula()
        sf.stat = self.stat
        sf.fs = self.fs
        sf.N_fs = {}
        for key in self.fs:
            if key in self.N_fs:
                sf.N_fs[key] = self.N_fs[key]
            else:
                sf.N_fs[key] = self.get_N_fs(key)
        sf.mods = self.mods
        sf.idx_num = self.idx_num
        sf.val_keys = self.val_keys
        sf.plot_info = self.plot_info
        sf.unikey = self.unikey

        # create folder
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)

        # pickle sf
        dill.settings['recurse'] = True
        dill.dump(sf, open(self.fullsave, 'wb'))

        return sf


    def __gen_formula_wrapper(self, keys):
        if keys == self.unikey:
            self.gen_formula()
        else:
            self.gen_formula(*keys)

    def __apply(self, f, vals):
        if len(vals) == 0:
            return f
        else:
            return f(*vals)

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

        vals = tuple(vals)

        return (keys, subval, vals)

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

    # make save name
    def mksvname(self):
        self.fullsave = self.folder + self.g.name + '_' + self.g.phi.name + '_' + str(self.stat) + '.sav'

    
class Moment(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.MOMENT
        self.idx_num = 1
        self.mksvname()

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

                            m = g.R[e, f, i, j].m_num(expr)
                            
                            g.moment_info[e, f, i, j, alpha] = c * m

                        # sum up
                        if e == f:
                            c2 = cf.ncrs(k, alpha) * (i * (g.d[e[0]][e[1]] + g.l[e])) ** (k - alpha[0] - alpha[1])
                        else:
                            c2 = cf.ncrs(k, alpha) * (g.d[e[i]][f[j]] + i * g.l[e] + j * g.l[f]) ** (k - alpha[0] - alpha[1])
                        
                        res += c2 * self.g.moment_info[e, f, i, j, alpha]
                        
        self.fs[k] = res
        return True
    
class CDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CDF
        self.val_keys = (x,)
        self.mksvname()

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

class PDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.PDF
        self.val_keys = (x,)
        self.mksvname()

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

                expr += const * tmp

        # self.mat_D = mat_D
        self.fs[self.unikey] = expand(expr)
        return True


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
        self.mksvname()

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

        self.fs[k, e] = expand(res)
        return True

class CCDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CCDF
        self.idx_num = 1
        self.val_keys = (p, x)
        self.mksvname()

    def gen_formula(self, e):
        if e in self.fs:
            return False

        print('generating the conditional cdf...')

        g = self.g
        expr = 0

        for f in g.edges():
            tmp = 0
            for (i, j) in g.two2:
                tmp += g.Rx[e, f, i, j].m_p(g.phi_pq / g.phi_p)

            expr += g.gy[f]* tmp

        self.fs[e] = expand(expr)
        return True



class CPDF(Formula):
    def __init__(self, g):
        super().__init__(g)
        self.stat = Stats.CPDF
        self.idx_num = 1
        self.val_keys = (p, x)
        self.mksvname()

    def gen_formula(self, e):
        if e in self.fs:
            return False

        print('generating the conditional pdf...')

        g = self.g
        expr = 0

        for f in g.edges():
            const = g.gy[f] / g.l[f]
            tmp = 0

            for (i, j) in g.two2:
                c0 = cf.eta(g.a[e, f, i, j], g.b[e, f, i, j], x)
                rg = g.L[e, f, i, j].bases[0]
                myeta = None
                if rg.bd[0] is False:
                    myeta = cf.etal
                elif rg.bd[1] is False:
                    myeta = cf.etar
                else:
                    myeta = cf.eta
                c1 = myeta(rg.pl, rg.pu, p)
                phis = (g.phi_pq / g.phi_p).subs(q, g.q[e, f, i, j])
                tmp += c0 * c1 * phis

            expr += const * tmp

        self.fs[e] = expand(expr)
        return True

