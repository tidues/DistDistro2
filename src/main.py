import basicInfo as bi
import formulas as fl
import Nformulas as Nfl
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

#def get_data(g, phi_pq, ddict):
#    fpath = '../data/'
#    g = bi.readGraph(fpath, gname)

    
    



# module test switch
switches = {
        'cycuni': 0,
        'cycnuni': 0,
        'clinuni': 0,
        'griduni': 1,
        'gridnuni': 1,
        'manhatton': 0
        }

if switches['cycuni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'gcyc'
    g = bi.readGraph(fpath, gname)

    phi_pq = 1

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
#    moment = fl.Moment(g)
#    m1 = moment.eval(1)
#    m2 = moment.eval(2)
#    m3 = moment.eval(3)
#    print(m1)
#    print(m2)
#    print(m3)
#    print(m2 - m1 ** 2)

    cdf = fl.CDF(g)
    #print(cdf.formula())
    cdf.plot()
    cdf.save()

    pdf = fl.PDF(g)
    #print(pdf.formula())
    pdf.plot()
    pdf.save()

if switches['cycnuni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'gcyc'
    g = bi.readGraph(fpath, gname)

    phi_pq = 36 * p * q * (1-p) * (1-q)

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
    #moment = fl.Moment(g)
    #m1 = moment.eval(1)
    #m2 = moment.eval(2)
    #m3 = moment.eval(3)
    #print(m1)
    #print(m2)
    #print(m3)
    #print(m2 - m1 ** 2)

    #cdf = fl.CDF(g)
    #print(cdf.formula())
    #cdf.plot()
    

    pdf = fl.PDF(g)
    #print(pdf.formula())
    pdf.plot()
    pdf.save()

if switches['clinuni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'gcli'
    g = bi.readGraph(fpath, gname)

    phi_pq = 6 * q * (1-q)

    bi.loadInfo(g, phi_pq=phi_pq, rational=False)
    #moment = fl.Moment(g)
    #m1 = moment.eval(1)
    #m2 = moment.eval(2)
    #m3 = moment.eval(3)
    #print(m1)
    #print(m2)
    #print(m3)
    #print(m2 - m1 ** 2)

    pdf = fl.PDF(g)
    #print(pdf.formula())
    pdf.plot()
    pdf.save()

    e = ('1', '2')
    p = 0.2
    cpdf = fl.CPDF(g)
    cpdf.plot(e, p)
    cpdf.save()

    e = ('1', '3')
    p = 0.5
    cpdf = fl.CPDF(g)
    cpdf.plot(e, p)
    cpdf.save()

    e = ('3', '4')
    p = 0
    cpdf = fl.CPDF(g)
    cpdf.plot(e, p)
    cpdf.save()

if switches['griduni'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'planar_side_10'

    #gname = 'g0'
    g = bi.readGraph(fpath, gname)

    #phi_pq = 1
    #phi_pq = 6 * q * (1-q)
    phi_pq = 36 * p * q * (1-p) * (1-q)

    # get moments info
    #bi.NloadInfo(g, phi_pq=phi_pq)
    #moment = Nfl.NMoment(g)
    #m0 = moment.eval(0)
    #print(m0)
    #m1 = moment.eval(1)
    #print(m1)
    #m2 = moment.eval(2)
    #print(m2)
    #print(m2 - m1 ** 2)
    #m3 = moment.eval(3)
    #print(m3)
    
    # get pdf info
    bi.loadInfo(g, phi_pq=phi_pq, rational=False)

    #pdf = fl.PDF(g)
    ##print(pdf.formula())
    #pdf.plot()

    #es = [('1', '2'), ('2', '3'), ('3', '4'), ('1', '2')]
    #es = [('1', '2'), ('5', '6'), ('45', '46')]
    es = [('45', '46')]
    ps = [0.6]
    cpdf = fl.CPDF(g)
    for idx, e in enumerate(es):
        cpdf.plot(e, ps[idx])

if switches['manhatton'] == 1:
    # graph file
    fpath = '../data/'
    gname = 'manhattan'
    g = bi.readGraph(fpath, gname)

    g = get_largest_component(g)

    phi_pq = 1

    bi.NloadInfo(g, phi_pq=phi_pq)

    # conditional moments
    cmoment = Nfl.NCMoment(g)
    ks = [1]
    es = [('1', '2')]
    ps = [0.5]
    for k in ks:
        for e in es:
            for p in ps:
                print(cmoment.eval(k, e, p))

