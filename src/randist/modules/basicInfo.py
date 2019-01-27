import networkx as nx
from .pyprelude import GraphGen as gg
from .commonFuncs import *
from sympy import *
import time


# reads a graph input file
def readGraph(fpath, gname):
    print('loading graph...')
    g = gg.GraphFromFile(fpath, gname).G
    g.name = gname
    return g

def basicInfo(g, phi, rational=False):
    print('generate basic info...')
    g.rat = rational
    g.phi = phi

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

    get_d(g)
    g.d_max = 0

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


#########################################################################
# help functions                                                        #
#########################################################################

# get shortest path matrix D
def get_d(g):
    g.d = dict(nx.all_pairs_dijkstra_path_length(g, weight = 'l'))
    print("load shortest distance")
    if g.rat:
        for i in g.nodes():
            for j in g.nodes():
                g.d[i][j] = rat(g.d[i][j], g.rat)


