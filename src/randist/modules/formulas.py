from . import commonFuncs as cf
from . import basicInfo as bi
from .pyprelude.Timer import Timer
from .pyprelude.Progress import Progress
from .pyprelude import EasyWriter as ew
from .enums import Stats
from sympy import symbols
from sympy.abc import x, p, q


# template for all region-wise operations
class Formula:
    def __init__(self, g):
        self.g = None
        self.stat = None
        self.keys = None
        self.idx_num = None
        self.resfolder = './results/'

    # set up the symbolic or numeric environment
    def sym_num_env(self, symbolic):
        self.symbolic = symbolic
        g = self.g
        if symbolic:
            # assgin symbolic functions
            from .Sregions import get_R, get_L, get_q
            from .commonFuncs import eta, etal, etar
            self.get_R = get_R
            self.get_L = get_L
            self.get_q = get_q
            self.eta = eta
            self.etal = etal
            self.etar = etar
            
            # assign symbolic phis
            g.phi_p = g.phi.phi_p_S
            g.phi_q = g.phi.phi_q_S
            g.phi_pq = g.phi.phi_pq_S
            g.phi_qcp = g.phi.phi_qcp_S
        else:
            # assgin symbolic functions
            from .Nregions import get_R, get_L, get_q
            from .numericFuncs import eta, etal, etar
            self.get_R = get_R
            self.get_L = get_L
            self.get_q = get_q
            self.eta = eta
            self.etal = etal
            self.etar = etar
            
            # assign symbolic phis
            g.phi_p = g.phi.phi_p_N
            g.phi_q = g.phi.phi_q_N
            g.phi_pq = g.phi.phi_pq_N
            g.phi_qcp = g.phi.phi_qcp_N

    def cond_ep(self, g, k, e, p_val, x_val, alphas, prog):
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

            # get entry coeff
            coeff = self.ef_coeff(px, py, lf)
            ef_res = 0

            for (i, j) in g.two2:
                # get region related info
                R = self.get_R(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val)
                ef_res += self.region_op(g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas)

            res += coeff * ef_res
        return res

    def cond_region(self, k=None, e=None, p_val=None, x_val=None):
        # generate extra params
        alphas = None
        if self.stat == Stats.MOMENT or self.stat == Stats.CMOMENT:
            alphas = cf.A_curl(2, k)

        g = self.g
        res = 0

        # see progress
        prog = Progress(g.number_of_edges() ** 2, response_time=10)

        # treat conditional
        if e is not None:
            return self.cond_ep(g, k, e, p_val, x_val, alphas, prog)

        for e in g.edges():
            res += self.cond_ep(g, k, e, p_val, x_val, alphas, prog)
        return res
    
    # combine results in e,f level
    def ef_coeff(self, px, py, lf):
        if self.stat == Stats.MOMENT:
            coeff = 1
        elif self.stat == Stats.CDF:
            coeff = px * py
        elif self.stat == Stats.PDF:
            coeff = px * py / lf
        elif self.stat == Stats.CMOMENT:
            coeff = 1
        elif self.stat == Stats.CCDF:
            coeff = py
        elif self.stat == Stats.CPDF:
            coeff = py / lf
        return coeff

    # save result when evaluate
    def save_res(self, keys, res):
        g = self.g
        fname = self.resfolder + g.name + '_' + g.phi.name + '_' + str(self.stat) + '.dat'
        val = str(keys) + '\t' + str(res) + '\n'
        ew.wFile(fname, val)

    # make dict from input
    def make_params(self, *params):
        p_dict = {}
        for idx, v in enumerate(params):
            if idx < len(self.keys):
                # use limits instead of undefined p value
                if self.keys[idx] == 'p_val':
                    v = self.adj_p_val(v)
                p_dict[self.keys[idx]] = v
            else:
                break
        return p_dict

    # use limit to replace undefined p
    def adj_p_val(self, p_val, eps=1e-7):
        if p_val > 0 and p_val < 1:
            return p_val
        elif p_val == 0:
            return eps
        elif p_val == 1:
            return 1 - eps
        else:
            raise Exception('p_val has to be in [0, 1]')

    ### abstract methods
    # operations in each region
    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
        pass


class Moment(Formula):
    def __init__(self, g, symbolic):
        super().__init__(g)
        self.g = g
        self.symbolic = symbolic
        self.stat = Stats.MOMENT
        self.keys = ['k']
        self.idx_num = -1
        self.sym_num_env(symbolic)

    def get_func(self, g, alpha):
        if self.symbolic:
            p, q = symbols('p,q')
            func = q ** alpha[0] * p ** alpha[1] * g.phi_pq
        else:
            func = lambda q, p: q ** alpha[0] * p ** alpha[1] * g.phi_pq(q, p)
        return func

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
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
            func = self.get_func(g, alpha)

            m = R.m(func)

            res += c * m
        return res


class CDF(Formula):
    def __init__(self, g, symbolic):
        super().__init__(g)
        self.g = g
        self.stat = Stats.CDF
        self.symbolic = symbolic
        self.keys = ['x_val']
        self.idx_num = 0
        self.sym_num_env(symbolic)

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
        return R.m(g.phi_pq)


class PDF(Formula):
    def __init__(self, g, symbolic):
        super().__init__(g)
        self.g = g
        self.stat = Stats.PDF
        self.symbolic = symbolic
        self.keys = ['x_val']
        self.idx_num = 0
        self.sym_num_env(symbolic)

    def get_func(self, q_func, g):
        if self.symbolic:
            p, q = symbols('p,q')
            func = g.phi_pq.subs(q, q_func)
        else:
            func = lambda p: g.phi_pq(q_func(p), p)
        return func

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
        c = self.eta(R.a_val, R.b_val, x_val)
        L = self.get_L(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val)
        q_func = self.get_q(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val)
        phis = self.get_func(q_func, g)
        m = L.m(phis)
        return c * m


class CMoment(Formula):
    def __init__(self, g, symbolic):
        super().__init__(g)
        self.g = g
        self.stat = Stats.CMOMENT
        self.symbolic = symbolic
        self.keys = ['k', 'e', 'p_val']
        self.idx_num = 2
        self.sym_num_env(symbolic)

    def get_func(self, g, p_val, alpha):
        if self.symbolic:
            p, q = symbols('p,q')
            func = q ** alpha[0] * g.phi_qcp
        else:
            func = lambda q: q ** alpha[0] * g.phi_qcp(q, p_val)
        return func

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
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

            func = self.get_func(g, p_val, alpha)

            m = R.m_p(func, p_val)
            res += c * m

        return res


class CCDF(Formula):
    def __init__(self, g, symbolic):
        super().__init__(g)
        self.g = g
        self.stat = Stats.CCDF
        self.symbolic = symbolic
        self.keys = ['e', 'p_val', 'x_val']
        self.idx_num = 1
        self.sym_num_env(symbolic)

    def get_func(self, g, p_val):
        if self.symbolic:
            p, q = symbols('p,q')
            func = g.phi_qcp
        else:
            func = lambda q: g.phi_qcp(q, p_val)
        return func

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
        func = self.get_func(g, p_val)
        return R.m_p(func, p_val)


class CPDF(Formula):
    def __init__(self, g, symbolic):
        super().__init__(g)
        self.g = g
        self.stat = Stats.CPDF
        self.symbolic = symbolic
        self.keys = ['e', 'p_val', 'x_val']
        self.idx_num = 1
        self.sym_num_env(symbolic)

    def get_func(self, g, q_func, p_val):
        if self.symbolic:
            p, q = symbols('p,q')
            func = g.phi_qcp.subs(q, q_func)
        else:
            func = g.phi_qcp(q_func(p_val), p_val)
        return func

    def region_op(self, g, e, f, i, j, le, lf, px, py, d, p1, p2, q1, q2, R, p_val, x_val, k, alphas):
        c0 = self.eta(R.a_val, R.b_val, x_val)
        L = self.get_L(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val)
        q_func = self.get_q(g, e, f, i, j, le, lf, d, p1, p2, q1, q2, x_val)
        myeta = None
        if L.bd[0] is False:
            myeta = self.etal
        elif L.bd[1] is False:
            myeta = self.etar
        else:
            myeta = self.eta
        c1 = myeta(L.pl, L.pu, p_val)
        phis = self.get_func(g, q_func, p_val)
        return c0 * c1 * phis


# the interface for numerical computation
class Numeric:
    # fl_cls is the formula class, e.g. Moment, CDF, CMoment, etc.
    def __init__(self, g, fl_cls):
        # gen formula
        self.formula = fl_cls(g, False)

    # evaluate value
    def eval(self, *params, save=True):
        p_dict = self.formula.make_params(*params)
        res = self.formula.cond_region(**p_dict)
        if save:
            self.formula.save_res(params, res)
        return res


# the interface for symbolic computation
class Symbolic:
    # fl_cls is the formula class, e.g. Moment, CDF, CMoment, etc.
    def __init__(self, g, fl_cls):

        # gen formula
        self.formula = fl_cls(g, True)

    # generate functions
    def gen_formula(self):
        pass

    # plot formulas
    def plot(self):
        pass

    # evaluate value
    def eval(self):
        pass

    # save formulas
    def save_formulas(self):
        pass


# the wrapper for getting formulas
class Formulas:
    # symbolic: use symbolic or numerical method, None is same as auto
    def __init__(self, gname, phi, fpath='../data/', rational=False):
        # read graph
        self.g = bi.readGraph(fpath, gname)
        # generate basic info
        bi.basicInfo(self.g, phi, rational)
        # dispatch map
        self.fl = {
                Stats.MOMENT: Moment,
                Stats.CDF: CDF,
                Stats.PDF: PDF,
                Stats.CMOMENT: CMoment,
                Stats.CCDF: CCDF,
                Stats.CPDF: CPDF
                }

    # get different stats
    def get_formula(self, stats, symbolic=None):
        g = self.g
        if Stats.is_member(stats) is False:
            raise Exception('stats must be value from the enum Stats')

        if symbolic is None:
            if stats == Stats.MOMENT or stats == Stats.CMOMENT:
                symbolic = False
            else:
                symbolic = True

        if symbolic:
            return Symbolic(g, self.fl[stats])
        else:
            return Numeric(g, self.fl[stats])
