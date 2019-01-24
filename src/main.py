import basicInfo as bi
import formulas as fl
import simulation as sl
from scipy import integrate
from sympy.abc import x, p, q
from sympy import lambdify, expand
import numericFuncs as nf
from networkx import is_connected
from enums import Stats
import numpy as np
import time
from commonFuncs import load_formulas, get_largest_component

# module test switch
switches = {
        'cycuni': 0,
        'cycnuni': 0,
        'clinuni': 0,
        'manhatton': 1
        }

if switches['cycuni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'gcyc'
    g = bi.readGraph(fpath, gname)

    phi_pq = 1

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
    moment = fl.Moment(g)
    m1 = moment.eval(1)
    m2 = moment.eval(2)
    m3 = moment.eval(3)
    print(m1)
    print(m2)
    print(m3)
    print(m2 - m1 ** 2)

    cdf = fl.CDF(g)
    #print(cdf.formula())
    cdf.plot()

    pdf = fl.PDF(g)
    #print(pdf.formula())
    pdf.plot()


if switches['cycnuni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'gcyc'
    g = bi.readGraph(fpath, gname)

    phi_pq = 36 * p * q * (1-p) * (1-q)

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
    moment = fl.Moment(g)
    m1 = moment.eval(1)
    m2 = moment.eval(2)
    m3 = moment.eval(3)
    print(m1)
    print(m2)
    print(m3)
    print(m2 - m1 ** 2)

    cdf = fl.CDF(g)
    #print(cdf.formula())
    cdf.plot()

    pdf = fl.PDF(g)
    #print(pdf.formula())
    pdf.plot()

if switches['clinuni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'gcli'
    g = bi.readGraph(fpath, gname)

    phi_pq = 6 * q * (1-q)

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
    moment = fl.Moment(g)
    m1 = moment.eval(1)
    m2 = moment.eval(2)
    m3 = moment.eval(3)
    print(m1)
    print(m2)
    print(m3)
    print(m2 - m1 ** 2)

    pdf = fl.PDF(g)
    #print(pdf.formula())
    pdf.plot()

    e = ('1', '2')
    p = 0.2
    cpdf = fl.CPDF(g)
    cpdf.plot(e, p)

    e = ('1', '3')
    p = 0.5
    cpdf = fl.CPDF(g)
    cpdf.plot(e, p)

    e = ('3', '4')
    p = 0
    cpdf = fl.CPDF(g)
    cpdf.plot(e, p)

if switches['manhatton'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'planar_side_11'
    g = bi.readGraph(fpath, gname)

    g = get_largest_component(g)

    print(len(g.nodes()))
    print(len(g.edges()))

    print('is connected: ', is_connected(g))

    phi_pq = 1

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
    moment = fl.Moment(g)
    m1 = moment.eval(1)
    m2 = moment.eval(2)
    m3 = moment.eval(3)
    print(m1)
    print(m2)
    print(m3)
    print(m2 - m1 ** 2)
    
    pdf = fl.PDF(g)
    pdf.plot()

    es = [('1', '2'), ('6', '7'), ('61', '62'), ('62', '63')]
    ps = [0, 0, 0, 0.5]
    cpdf = fl.CPDF(g)
    for idx, e in enumerate(es):
        cpdf.plot(e, ps[idx])

    #moment = fl.Moment(g)
    #m1 = moment.eval(1)
    #m2 = moment.eval(2)
    #m3 = moment.eval(3)
    #print(m1)
    #print(m2)
    #print(m3)
    #print(m2 - m1 ** 2)

    #cdf = fl.CDF(g)
    ##print(cdf.formula())
    #cdf.plot()

    #pdf = fl.PDF(g)
    ##print(pdf.formula())
    #pdf.plot()
