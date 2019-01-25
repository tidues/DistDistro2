import commonFuncs as cf
import basicInfo as bi
from Nregions import get_R
from pyprelude.Timer import Timer
from pyprelude.Progress import Progress

class NMoment:
    def __init__(self, g):
        self.g = g

    def eval(self, k, dis_time=10):
        g = self.g
        alphas = cf.A_curl(2, k)
        res = 0

        # see progress
        prog = Progress(g.number_of_edges() ** 2)

        for e in g.edges():
            # get e related info
            le = float(g.edges[e]['l'])
            px = float(g.edges[e]['x'])
            for f in g.edges():

                # display progress
                prog.count()

                # get f related info
                lf = float(g.edges[f]['l'])
                py = float(g.edges[f]['y'])
                # get entry ef related info
                d, p1, p2, q1, q2 = bi.entry_info(self.g, e, f, le, lf)

                for (i, j) in g.two2:
                    # get region related info
                    R = get_R(self.g, e, f, i, j, le, lf, d, p1, p2, q1, q2)

                    # get alpha related info
                    for alpha in alphas:
                        c0 = cf.ncrs(k, alpha) * px * py * lf ** alpha[0] * le ** alpha[1]
                        if e == f:
                            w = (i + j + 1) * alpha[0] + (i + j + 2) * alpha[1]
                            c1 = (-1) ** w * (i * (d + le)) ** (k - alpha[0] - alpha[1])
                        else:
                            w = j *  alpha[0] + i * alpha[1]
                            c1 = (-1) ** w * (d[i, j] + i * le + j * lf) ** (k - alpha[0] - alpha[1])
                        c = c0 * c1
                        func = lambda q, p: q ** alpha[0] * p ** alpha[1] * g.phi_pq(q, p)

                        m = R.m(func)

                        res += c * m

        return res




class NCMoment:
    def __init__(self, g):
        self.g = g

    def eval(self, k, e, p_val, dis_time=10):
        g = self.g
        alphas = cf.A_curl(2, k)
        res = 0

        # see progress
        prog = Progress(g.number_of_edges())

        # get e related info
        le = float(g.edges[e]['l'])
#        px = float(g.edges[e]['x'])
        
        for f in g.edges():

            # display
            prog.count()

            # get f related info
            lf = float(g.edges[f]['l'])
            py = float(g.edges[f]['y'])

            # get entry ef related info
            d, p1, p2, q1, q2 = bi.entry_info(self.g, e, f, le, lf)

            for (i, j) in g.two2:
                # get region related info
                R = get_R(self.g, e, f, i, j, le, lf, d, p1, p2, q1, q2)

                for alpha in alphas:
                    c0 = cf.ncrs(k, alpha) * py * lf ** alpha[0] * le ** alpha[1]
                    if e == f:
                        w = (i + j + 1) * alpha[0] + (i + j + 2) * alpha[1]
                        c1 = (-1) ** w * (i * (d + le)) ** (k - alpha[0] - alpha[1])
                    else:
                        w = j *  alpha[0] + i * alpha[1]
                        c1 = (-1) ** w * (d[i, j] + i * le + j * lf) ** (k - alpha[0] - alpha[1])
                    c = c0 * c1 * p_val ** alpha[1]

                    func = lambda q: q ** alpha[0] * g.phi_qcp(q, p_val)

                    m = R.m_p(func, p_val)
                    res += c * m
        return res

