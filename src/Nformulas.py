import commonFuncs as cf
import basicInfo as bi
from Nregions import get_R
from pyprelude.Timer import Timer
from pyprelude.Progress import Progress


# template for all region-wise operations
class NFormulas:
    def __init__(self):
        self.g = None

    def cond_ep(self, g, e, p_val, x_val, prog, **params):
        res = 0
        # get e related info
        le = float(g.edges[e]['l'])
        px = float(g.edges[e]['x'])
        for f in g.edges():

            prog.count()

            # get f related info
            lf = float(g.edges[f]['l'])
            py = float(g.edges[f]['y'])
            # get entry ef related info
            d, p1, p2, q1, q2 = bi.entry_info(self.g, e, f, le, lf)

            for (i, j) in g.two2:
                # get region related info
                R = get_R(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val)
                res += self.region_op(g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, **params)

        return res

    def cond_region(self, e=None, p_val=None, x_val=None, **params):
        g = self.g
        res = 0

        # see progress
        prog = Progress(g.number_of_edges() ** 2)

        # treat conditional
        if e is not None:
            return self.cond_ep(g, e, p_val, x_val, prog, **params)

        for e in g.edges():
            res += self.cond_ep(g, e, p_val, x_val, prog, **params)
        return res
    
    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, **params):
        pass


class NMoment(NFormulas):
    def __init__(self, g):
        self.g = g

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, **params):
        alphas = params['alphas']
        k = params['k']
        res = 0
        # get alpha related info
        for alpha in alphas:
            c0 = cf.ncrs(k, alpha) * px * py * lf ** alpha[0] * le ** alpha[1]
            if e == f:
                w = (i + j + 1) * alpha[0] + (i + j + 2) * alpha[1]
                c1 = (-1) ** w * (i * (d + le)) ** (k - alpha[0] - alpha[1])
            else:
                w = j * alpha[0] + i * alpha[1]
                c1 = (-1) ** w * (d[i, j] + i * le + j * lf) ** (k - alpha[0] - alpha[1])
            c = c0 * c1
            func = lambda q, p: q ** alpha[0] * p ** alpha[1] * g.phi_pq(q, p)

            m = R.m(func)

            res += c * m
        return res

    def eval(self, k):
        alphas = cf.A_curl(2, k)
        return self.cond_region(alphas=alphas, k=k)


class NCDF(NFormulas):
    def __init__(self, g):
        self.g = g

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, **params):
        return px * py * R.m(g.phi_pq)

    def eval(self, x_val):
        return self.cond_region(x_val=x_val)


class NCMoment(NFormulas):
    def __init__(self, g):
        self.g = g

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, **params):
        alphas = params['alphas']
        k = params['k']
        res = 0
        for alpha in alphas:
            c0 = cf.ncrs(k, alpha) * py * lf ** alpha[0] * le ** alpha[1]
            if e == f:
                w = (i + j + 1) * alpha[0] + (i + j + 2) * alpha[1]
                c1 = (-1) ** w * (i * (d + le)) ** (k - alpha[0] - alpha[1])
            else:
                w = j * alpha[0] + i * alpha[1]
                c1 = (-1) ** w * (d[i, j] + i * le + j * lf) ** (k - alpha[0] - alpha[1])
            c = c0 * c1 * p_val ** alpha[1]

            func = lambda q: q ** alpha[0] * g.phi_qcp(q, p_val)

            m = R.m_p(func, p_val)
            res += c * m

        return res


    def eval(self, k, e, p_val):
        alphas = cf.A_curl(2, k)
        return self.cond_region(e=e, p_val=p_val, k=k, alphas=alphas)


class NCCDF(NFormulas):
    def __init__(self, g):
        self.g = g

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, **params):
        func = lambda q: g.phi_qcp(q, p_val)
        return py * R.m_p(func, p_val)

    def eval(self, e, p_val, x_val):
        return self.cond_region(e=e, p_val=p_val, x_val=x_val)
