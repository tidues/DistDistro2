import networkx as nx
import pyprelude.GraphGen as gg
from commonFuncs import *
from regions import get_R, get_L
from sympy import *

# reads a graph input file
def readGraph(fpath, gname):
    return gg.GraphFromFile(fpath, gname).G

# generate all basic infomation
def loadInfo(g, phi, rational = False):
    # if using rational number to handle division
    g.rat = rational

    # save pdf
    g.phi = toSym(phi)
    g.moment_info = {}

    # index set
    g.two = (0,1)
    g.two2 = [(i, j) for i in g.two for j in g.two]

    # generate ordered edges
    # g.oedges = []
    # sort_edges(g)

    # edge lengths
    g.l = {}
    g.gx = {}  # the pmf of X
    g.gy = {}  # the pmf of Y

    for e in g.edges():
        g.l[e] = rat(g.edges[e]['l'], g.rat)
        g.gx[e] = rat(g.edges[e]['x'], g.rat)
        g.gy[e] = rat(g.edges[e]['y'], g.rat)

    # all pair shortest distance
    get_d(g)

    # compute p1 p2 q1 q2 for all edge pairs
    g.p1 = {}
    g.p2 = {}
    g.q1 = {}
    g.q2 = {}

    for e in g.edges():
        for f in g.edges():
            if e != f:
                get_pq_const(e, f, g)

    # the domain boundary info (a and b)
    g.a = {}
    g.b = {}
    for e in g.edges():
        for f in g.edges():
            get_ab(e, f, g)

    # the region info R, L, and q
    g.Rx = {}
    g.R = {}
    g.L = {}
    g.q = {}
    for e in g.edges():
        for f in g.edges():
            get_RLq(e, f, g)


# assign pdf phi_{P,Q}
def update_phi(phi):
    g.phi = toSym(phi)
    g.moment_info = {}


#########################################################################
# help functions                                                        #
#########################################################################

# sort edges
#def sort_edges(g):
#    for e in g.edges():
#        g.oedges.append(e_repr(e))


# get shortest path matrix D
def get_d(g):
    g.d = dict(nx.all_pairs_dijkstra_path_length(g, weight = 'l'))
    for i in g.nodes():
        for j in g.nodes():
            g.d[i][j] = rat(g.d[i][j], g.rat)

# get x bounds info
def get_ab(e, f, g):
    if e == f:
        duv = g.d[e[0]][e[1]]
        le = g.l[e]
        g.a[(e, e, 0, 0)] = rat(0, g.rat)
        g.a[(e, e, 0, 1)] = rat(0, g.rat)
        g.a[(e, e, 1, 0)] = duv
        g.a[(e, e, 1, 1)] = duv

        for (i,j) in g.two2:
            g.b[(e, f, i, j)] = (duv + le) / 2
    else:
        le = g.l[e]
        lf = g.l[f]

        for (i,j) in g.two2:
            g.a[e, f, i, j] = g.d[e[i]][f[j]]
            g.b[e, f, i, j] = (le + lf + 
                    min(g.d[e[0]][f[0]] + g.d[e[1]][f[1]], 
                        g.d[e[0]][f[1]] + g.d[e[1]][f[0]])) / 2

# get p1 p2 q1 q2
def get_pq_const(e, f, g):
    le = g.l[e]
    lf = g.l[f]
    g.p1[e, f] = (-g.d[e[0]][f[0]] + g.d[e[1]][f[0]] + le) / (2 * le)
    g.p2[e, f] = (-g.d[e[0]][f[1]] + g.d[e[1]][f[1]] + le) / (2 * le)
    g.q1[e, f] = (-g.d[e[0]][f[0]] + g.d[e[0]][f[1]] + lf) / (2 * lf)
    g.q2[e, f] = (-g.d[e[1]][f[0]] + g.d[e[1]][f[1]] + lf) / (2 * lf)


# get region info R L q
def get_RLq(e, f, g):
    
    # the delta
    d = {}
    for (i, j) in g.two2:
        d[i, j] = g.d[e[i]][f[j]]

    le = g.l[e]
    lf = g.l[f]
    duv = g.d[e[0]][e[1]]

    p1 = None
    p2 = None
    q1 = None
    q2 = None

    if e != f:
        p1 = g.p1[e, f]
        p2 = g.p2[e, f]
        q1 = g.q1[e, f]
        q2 = g.q2[e, f]

    get_R(e, f, g, le, lf, duv, d, p1, p2, q1, q2)
    get_L(e, f, g, le, lf, duv, d, p1, p2, q1, q2)
    get_q(e, f, g, le, lf, duv, d, p1, p2, q1, q2)

# the q function for pdf
def get_q(e, f, g, le, lf, duv, d, p1, p2, q1, q2):
    x = Symbol('x')
    p = Symbol('p')

    if e == f:
        for (i, j) in g.two2:
            a = g.a[e, f, i, j]
            b = g.b[e, f, i, j]
            xp = theta(a, b, x)
            if (i, j) == (0, 0):
                g.q[e, e, i, j] = p - xp / le
            elif (i, j) == (0, 1):
                g.q[e, e, i, j] = p + xp / le
            elif (i, j) == (1, 0):
                g.q[e, e, i, j] = p - 1 + (xp - duv) / le
            else:
                g.q[e, e, i, j] = p + 1 - (xp - duv) / le
    else:
        for (i, j) in g.two2:
            a = g.a[e, f, i, j]
            b = g.b[e, f, i, j]
            xp = theta(a, b, x)
            if (i, j) == (0, 0):
                g.q[e, f, i, j] = (xp - le * p - d[i, j]) / lf
            elif (i, j) == (0, 1):
                g.q[e, f, i, j] = (- xp + le * p + d[i, j] + lf) / lf
            elif (i, j) == (1, 0):
                g.q[e, f, i, j] = (xp + le * p - d[i, j] - le) / lf
            else:
                g.q[e, f, i, j] = (- xp - le * p + d[i, j] + le + lf) / lf

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



