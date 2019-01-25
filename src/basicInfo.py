import networkx as nx
import pyprelude.GraphGen as gg
from commonFuncs import *
from regions import get_R, get_L
from sympy import *
import time

# reads a graph input file
def readGraph(fpath, gname):
    print('loading graph...')
    return gg.GraphFromFile(fpath, gname).G

def NloadInfo(g, phi_p=None, phi_q=None, phi_pq=None):
    print('generate basic info...')
    start = time.time()
    # make pdfs
    g.phi_p = phi_p
    g.phi_q = phi_q
    g.phi_pq = phi_pq
    N_gen_phi(g)

    # check graph info
    check_res = gcheck(g)
    if check_res['total'] is not True:
        raise Exception('Input did not pass all test: ', check_res)

    g.moment_info = {}
    g.cmoment_info = {}

    # index set
    g.two = (0,1)
    g.two2 = [(i, j) for i in g.two for j in g.two]
    end = time.time()
    print('basic info:', end - start)

    # all pair shortest distance
    print('get shortest path info')
    start = time.time()
    N_get_d(g)
    end = time.time()
    print('get_d:', end - start)

# entry info
def entry_info(g, e, f, le, lf):
    p1 = None
    p2 = None
    q1 = None
    q2 = None
    if e == f:
        d = g.d[e[0]][e[1]]
    else:
        d = {}
        for (i, j) in g.two2:
            d[i, j] = g.d[e[i]][f[j]]
        p1 = (-d[0, 0] + d[1, 0] + le) / (2 * le)
        p2 = (-d[0, 1] + d[1, 1] + le) / (2 * le)
        q1 = (-d[0, 0] + d[0, 1] + lf) / (2 * lf)
        q2 = (-d[1, 0] + d[1, 1] + lf) / (2 * lf)

    return(d, p1, p2, q1, q2)


# generate all basic infomation
def loadInfo(g, phi_p=None, phi_q=None, phi_pq=None, rational = False):

    print('generate basic info...')
    start = time.time()
    # make pdfs
    g.phi_p = phi_p
    g.phi_q = phi_q
    g.phi_pq = phi_pq
    gen_phi(g)

    # check graph info
    check_res = gcheck(g)
    if check_res['total'] is not True:
        raise Exception('Input did not pass all test: ', check_res)

    g.moment_info = {}
    g.cmoment_info = {}

    # if using rational number to handle division
    g.rat = rational

    # index set
    g.two = (0,1)
    g.two2 = [(i, j) for i in g.two for j in g.two]
    end = time.time()
    print('basic info:', end - start)

    ## create entry-wise info
    #for e in g.edges():
    #    le = rat(g.edges[e]['l'], g.rat)
    #    for f in g.edges():
    #        lf = rat(g.edges[f]['l'], g.rat)

    # edge lengths
    print('load edges properties')
    start = time.time()
    g.l = {}
    g.gx = {}  # the pmf of X
    g.gy = {}  # the pmf of Y

    for e in g.edges():
        g.l[e] = rat(g.edges[e]['l'], g.rat)
        g.gx[e] = rat(g.edges[e]['x'], g.rat)
        g.gy[e] = rat(g.edges[e]['y'], g.rat)
    end = time.time()
    print('load eddges:', end - start)

    # all pair shortest distance
    print('get shortest path info')
    start = time.time()
    get_d(g)
    end = time.time()
    print('get_d:', end - start)

    # compute p1 p2 q1 q2 for all edge pairs
    print('calculate pq_const')
    start = time.time()
    g.p1 = {}
    g.p2 = {}
    g.q1 = {}
    g.q2 = {}

    for e in g.edges():
        for f in g.edges():
            if e != f:
                get_pq_const(e, f, g)
    end = time.time()
    print('pq_const:', end - start)

    # the domain boundary info (a and b)
    print('calculate ab_bounds')
    start = time.time()
    g.a = {}
    g.b = {}
    g.d_max = 0
    for e in g.edges():
        for f in g.edges():
            get_ab(e, f, g)
    end = time.time()
    print('ab_bounds:', end - start)

    # the region info R, L, and q
    print('make regions')
    start = time.time()
    g.Rx = {}
    g.R = {}
    g.L = {}
    g.q = {}
    cnt = 0
    tR = 0
    tL = 0
    tq = 0
    for e in g.edges():
        for f in g.edges():
            cnt += 1
            if cnt % 1000 == 0:
                print(cnt)
            t0, t1, t2 = get_RLq(e, f, g)
            tR += t0
            tL += t1
            tq += t2
    end = time.time()
    print('regions:', end - start)
    print('tR:', tR)
    print('tL:', tL)
    print('tq:', tq)


# assign pdf phi_{P,Q}
def update_phi(phi):
    g.phi = sympify(phi)
    g.moment_info = {}


#########################################################################
# help functions                                                        #
#########################################################################

# sort edges
#def sort_edges(g):
#    for e in g.edges():
#        g.oedges.append(e_repr(e))

# generate all phis
def gen_phi(g):
    if g.phi_p is not None and g.phi_q is not None:
        g.phi_p = sympify(g.phi_p)
        g.phi_q = sympify(g.phi_q)
        g.phi_pq = g.phi_p * g.phi_q
    elif g.phi_pq is not None:
        g.phi_pq = sympify(g.phi_pq)
        p, q = symbols('p,q')
        g.phi_p = integrate(g.phi_pq, (q, 0 , 1))
        g.phi_q = integrate(g.phi_pq, (p, 0 , 1))
    else:
        raise Exception('please provide both phi_p and phi_q, or the joint pdf phi_pq.')
    g.phi_pq_N = lambdify((q, p), g.phi_pq)

# generate numerica phi
def N_gen_phi(g):
    p, q = symbols('p,q')
    if g.phi_p is not None and g.phi_q is not None:
        g.phi_p = sympify(g.phi_p)
        g.phi_q = sympify(g.phi_q)
        g.phi_pq = g.phi_p * g.phi_q
        g.phi_qcp = g.phi_q
    elif g.phi_pq is not None:
        g.phi_pq = sympify(g.phi_pq)
        g.phi_p = integrate(g.phi_pq, (q, 0 , 1))
        g.phi_q = integrate(g.phi_pq, (p, 0 , 1))
        g.phi_qcp = g.phi_pq / g.phi_p
    else:
        raise Exception('please provide both phi_p and phi_q, or the joint pdf phi_pq.')

    g.phi_p = lambdify(p, g.phi_p)
    g.phi_q = lambdify(q, g.phi_q)
    g.phi_pq = lambdify((q, p), g.phi_pq)
    g.phi_pq_N = g.phi_pq
    g.phi_qcp = lambdify((q, p), g.phi_qcp)

# get shortest path matrix D
def get_d(g):
    g.d = dict(nx.all_pairs_dijkstra_path_length(g, weight = 'l'))
    print("load shortest distance")
    if g.rat:
        for i in g.nodes():
            for j in g.nodes():
                g.d[i][j] = rat(g.d[i][j], g.rat)

# get shortest path matrix D
def N_get_d(g):
    g.d = dict(nx.all_pairs_dijkstra_path_length(g, weight = 'l'))

# get x bounds info
def get_ab(e, f, g):
    if e == f:
        duv = g.d[e[0]][e[1]]
        le = g.l[e]
        g.a[(e, e, 0, 0)] = rat(0, g.rat)
        g.a[(e, e, 0, 1)] = rat(0, g.rat)
        g.a[(e, e, 1, 0)] = g.d[e[0]][e[1]]
        g.a[(e, e, 1, 1)] = g.d[e[0]][e[1]]

        for (i,j) in g.two2:
            bval = (g.d[e[0]][e[1]] + g.l[e]) / 2
            g.b[(e, f, i, j)] = bval

            if bval > g.d_max:
                g.d_max = bval
    else:

        for (i,j) in g.two2:
            g.a[e, f, i, j] = g.d[e[i]][f[j]]
            bval = (g.l[e] + g.l[f] + 
                    min(g.d[e[0]][f[0]] + g.d[e[1]][f[1]], 
                        g.d[e[0]][f[1]] + g.d[e[1]][f[0]])) / 2
            g.b[e, f, i, j] = bval
            if bval > g.d_max:
                g.d_max = bval

# get p1 p2 q1 q2
def get_pq_const(e, f, g):
    g.p1[e, f] = (-g.d[e[0]][f[0]] + g.d[e[1]][f[0]] + g.l[e]) / (2 * g.l[e])
    g.p2[e, f] = (-g.d[e[0]][f[1]] + g.d[e[1]][f[1]] + g.l[e]) / (2 * g.l[e])
    g.q1[e, f] = (-g.d[e[0]][f[0]] + g.d[e[0]][f[1]] + g.l[f]) / (2 * g.l[f])
    g.q2[e, f] = (-g.d[e[1]][f[0]] + g.d[e[1]][f[1]] + g.l[f]) / (2 * g.l[f])


# get region info R L q
def get_RLq(e, f, g):
    
    # the delta
    d = {}
    for (i, j) in g.two2:
        d[i, j] = g.d[e[i]][f[j]]

    if e != f:
        start = time.time()
        get_R(e, f, g, g.l[e], g.l[f], g.d[e[0]][e[1]], d, g.p1[e,f], g.p2[e,f], g.q1[e,f], g.q2[e,f])
        end = time.time()
        tR = end - start
        start = time.time()
        get_L(e, f, g, g.l[e], g.l[f], g.d[e[0]][e[1]], d, g.p1[e,f], g.p2[e,f], g.q1[e,f], g.q2[e,f])
        end = time.time()
        tL = end - start
        start = time.time()
        get_q(e, f, g, g.l[e], g.l[f], g.d[e[0]][e[1]], d, g.p1[e,f], g.p2[e,f], g.q1[e,f], g.q2[e,f])
        end = time.time()
        tq = end - start
    else:
        start = time.time()
        get_R(e, f, g, g.l[e], g.l[f], g.d[e[0]][e[1]], d, None, None, None, None)
        end = time.time()
        tR = end - start
        start = time.time()
        get_L(e, f, g, g.l[e], g.l[f], g.d[e[0]][e[1]], d, None, None, None, None)
        end = time.time()
        tL = end - start
        start = time.time()
        get_q(e, f, g, g.l[e], g.l[f], g.d[e[0]][e[1]], d, None, None, None, None)
        end = time.time()
        tq = end - start

    return (tR, tL, tq)


# the q function for pdf
def get_q(e, f, g, le, lf, duv, d, p1, p2, q1, q2):
    x = Symbol('x')
    p = Symbol('p')

    if e == f:
        for (i, j) in g.two2:
            if (i, j) == (0, 0):
                g.q[e, e, i, j] = p - x / le
            elif (i, j) == (0, 1):
                g.q[e, e, i, j] = p + x / le
            elif (i, j) == (1, 0):
                g.q[e, e, i, j] = p - 1 + (x - duv) / le
            else:
                g.q[e, e, i, j] = p + 1 - (x - duv) / le
    else:
        for (i, j) in g.two2:
            if (i, j) == (0, 0):
                g.q[e, f, i, j] = (x - le * p - d[i, j]) / lf
            elif (i, j) == (0, 1):
                g.q[e, f, i, j] = (- x + le * p + d[i, j] + lf) / lf
            elif (i, j) == (1, 0):
                g.q[e, f, i, j] = (x + le * p - d[i, j] - le) / lf
            else:
                g.q[e, f, i, j] = (- x - le * p + d[i, j] + le + lf) / lf

#############################################################################
# verify functions
#############################################################################
def g_info_show(g):
    checkLst = {
            'phi': g.phi,
            'two2': g.two2,
            'edges()': g.edges(),
            'l': g.l,
            'gx': g.gx,
            'gy': g.gy,
            'd': g.d,
            'p1': g.p1,
            'p2': g.p2,
            'q1': g.q1,
            'q2': g.q2,
            'a': g.a,
            'b': g.b,
            'R': g.R,
            'Rx': g.Rx,
            'L': g.L,
            'q': g.q
            }

    for info in checkLst:
        print(info, '\n')
        print(checkLst[info])
        print('\n\n')


def check_R_measure_one(g, verb=False):
    for e in g.edges():
        for f in g.edges():
            if e <= f:
                total = 0
                for (i, j) in g.two2:
                    if verb is True:
                        for b in g.R[e, f, i, j].bases:
                            print(e, f, i, j, b.m_adj(1))
                    total += g.R[e, f, i, j].m(1)
                print(e, f, total)

def print_R(g):
    for e in g.edges():
        for f in g.edges():
            if e <= f:
                for (i, j) in g.two2:
                    print('R of ', str((e, f, i, j)), 'is :')
                    g.R[e, f, i, j].print()



