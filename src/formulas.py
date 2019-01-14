from enums import Stats
import commonFuncs as cf
from sympy import *
from sympy.abc import p, q, x
from sympy.plotting import plot
from plot import plot1d

class Formula:
    def __init__(self, g):
        self.g = g
        self.f = None
        self.f_lambda = None
        self.mods = ['numpy', {'theta': cf.theta, 'eta': cf.eta}]

    def __gen_formula(self):
        # realized by the inherited classes
        pass

    def eval(self, *params):
        # realized by the inherited classes
        pass

    def formula(self):
        if self.stat != Stats.MOMENT:
            print(self.f)
        else:
            print('moment function has no symbolic formula.')

    def plot(self, method=0, step=0.1):
        if self.stat != Stats.MOMENT:
            if method == 0:
                try:
                    plot(self.f, (x, -1, self.g.d_max + 1))
                except:
                    self.plot(method=1, step=step)
            else:
                self.f_lambda = lambdify(x, self.f, modules=self.mods)
                plot1d(self.f_lambda, -1, self.g.d_max + 1, step)
        else:
            print('moment function cannot be plotted.')



class Moment(Formula):
    def __init__(self, g):
        self.stat = Stats.MOMENT
        super().__init__(g)
        self.__gen_formula()

    def __gen_formula(self):
        print('generating the moment function...')
        def moment(k):
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
                                expr = q ** alpha[0] * p ** alpha[1] * g.phi
                                m = g.R[e, f, i, j].m(expr)
                                g.moment_info[e, f, i, j, alpha] = c * m

                            # sum up
                            if e == f:
                                c2 = cf.ncrs(k, alpha) * (i * (g.d[e[0]][e[1]] + g.l[e])) ** (k - alpha[0] - alpha[1])
                            else:
                                c2 = cf.ncrs(k, alpha) * (g.d[e[i]][f[j]] + i * g.l[e] + j * g.l[f]) ** (k - alpha[0] - alpha[1])
                            
                            res += c2 * self.g.moment_info[e, f, i, j, alpha]
            return res
        self.f = moment

    def eval(self, k):
        return self.f(k)

    


class CDF(Formula):
    def __init__(self, g):
        self.stat = Stats.CDF
        super().__init__(g)
        self.__gen_formula()

    def __gen_formula(self):
        print('generating the cdf...')

        g = self.g
        # x = Symbol('x')
        expr = 0

        for e in g.edges():
            for f in g.edges():
                gxy = g.gx[e] * g.gy[f]
                tmp = 0
                for (i, j) in g.two2:
                    tmp += g.Rx[e, f, i, j].m(g.phi)

                expr += gxy * tmp

        self.f = expand(expr)

    def eval(self, x_val):
        # x = Symbol('x')
        return self.f.subs(x, x_val)


class PDF(Formula):
    def __init__(self, g):
        self.stat = Stats.PDF
        super().__init__(g)
        self.__gen_formula()

    def __gen_formula(self):
        print('generating the pdf...')
        # x = Symbol('x')
        # p = Symbol('p')
        # q = Symbol('q')
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
                    phis = g.phi.subs(q, g.q[e, f, i, j])
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
        self.f = expand(expr)
        #self.f = lambdify(x, self.expr, modules=self.modules)

    def eval(self, x_val):
        # x = Symbol('x')
        return self.f.subs(x, x_val)




        








