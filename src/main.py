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
from commonFuncs import load_formulas

# module test switch
switches = {
        Stats.MOMENT: 0,
        Stats.CDF: 0,
        Stats.PDF: 0,
        Stats.CMOMENT: 0,
        Stats.CCDF: 0,
        Stats.CPDF: 0,
        Stats.SIMULATION: 0,
        Stats.SAVE: 1,
        Stats.TIMING: 0
        }

if switches[Stats.SAVE] == 1:
    moment = load_formulas(Stats.MOMENT)
    print(moment.stat)
    print(moment.eval(0))
    print(moment.eval(1))
    print(moment.eval(2))
    cdf = load_formulas(Stats.CDF)
    cdf.plot()

# graph file
fpath = '../data/'
gname = 'g3'

# read graph
g = bi.readGraph(fpath, gname)

# connected test
print('is connected: ', is_connected(g))

# set either joint pdf phi_pq, 
# or independent pdfs phi_p and phi_q

phi_pq = 1
# phi_pq = 4 * p * q

# load basic info
bi.loadInfo(g, phi_pq=phi_pq, rational=False)

# test values
valLst = [-0.1, -0.01, -0.001, 0, 1, 2, 3, -1, 9, 9.5, 9.9, 10]

# simulation 
if switches[Stats.SIMULATION] == 1:
    sim = sl.events(g)
    print(sim.E_d(f=lambda x: x**3))

# moments stats
if switches[Stats.MOMENT] == 1:
    moment = fl.Moment(g)
    # print(moment.eval(0))
    print(moment.eval(1))
    print(moment.eval(2))
    sm = moment.save()
    print(sm.eval(1))
    print(sm.eval(2))
    print(sm.eval(3))

# cdf stats
if switches[Stats.CDF] == 1:
    cdf = fl.CDF(g)
    cdf.formula()
    cdf.plot()
    cdf.save()
    mcdf = load_formulas(Stats.CDF)
    cdf.plot()
    

# pdf stats
if switches[Stats.PDF] == 1:
    pdf = fl.PDF(g)
    f = pdf.formula()
    print(nf.pdf_check(f, (x, 0, 12)))
    pdf.plot()
    for x in valLst:
        print(x, '\t', pdf.eval(x))

# conditional moments stats
if switches[Stats.CMOMENT] == 1:
    e = ('1', '2')
    f = ('1', '2')
    i, j = (0, 1)
    cmoment = fl.CMoment(g)
    print(cmoment.eval(0, ('1','2'), 0))
    print(cmoment.eval(0, ('1','2'), 0.5))
    print(cmoment.eval(0, ('1','2'), 1))
    print(cmoment.eval(1, ('1','2'), 0))
    print(cmoment.eval(1, ('1','2'), 0.5))
    print(cmoment.eval(1, ('1','2'), 1))
    print(cmoment.eval(2, ('1','2'), 0))
    print(cmoment.eval(2, ('1','2'), 0.5))
    print(cmoment.eval(2, ('1','2'), 1))
    cmoment.plot(0, ('1', '2'))
    cmoment.plot(1, ('1', '2'))
    cmoment.plot(2, ('1', '2'))


# conditional cdf
if switches[Stats.CCDF] == 1:
    e = ('1', '2')
    f = ('2', '3')
    ccdf = fl.CCDF(g)
    ccdf.plot(e, 0)
    ccdf.plot(e, 0.2)
    ccdf.plot(e, 0.4)
    ccdf.plot(e, 0.5)
    ccdf.plot(e, 0.8)
    ccdf.plot(e, 1)

# conditional pdf
if switches[Stats.CPDF] == 1:
    e = ('1', '2')
    f = ('2', '3')
    cpdf = fl.CPDF(g)
    plst = [0, 0.2, 0.4, 0.5, 0.8, 1]
    for p in plst:
        cpdf.plot(e, p)
        print(nf.pdf_check(cpdf.formula(e, p), (x, 0, g.d_max)))


# time comparison
# take an region
if switches[Stats.TIMING] == 1:
    e = ('1', '2')
    f = ('3', '4')
    mods = ['numpy', {'theta': nf.theta, 'eta': nf.eta, 'etal': nf.etal, 'etar': nf.etar}]

# symbolic integration
    # flow1: symbolic -> sum -> lambdify -> evaluation
    # symbolic + sum
    print('starting flow 1')
    start = time.time()
    expr = 0
    for t in range(1000):
        for (i, j) in g.two2:
            expr += g.Rx[e, f, i, j].m(phi_pq)
    end = time.time()
    f1sym = end - start
    
    # lambdify
    start = time.time()
    expr = expand(expr)
    myf = lambdify(x, expr, modules=mods)
    end = time.time()
    f1lbd = end - start

    # evaluation
    start = time.time()
    val1 = 0
    for x_val in np.arange(4, 5, 1):
        val1 += myf(x_val)
    end = time.time()
    f1eval = end - start
    res1 = [('sym', f1sym), ('lbd', f1lbd), ('eval', f1eval)]

    
    # flow2: symbolic -> lambdify -> evaluation -> sum
    # symbolic + lambdify
    #start = time.time()
    #myf = {}
    #for t in range(1):
    #    for (i, j) in g.two2:
    #        myf[i, j, t] = lambdify(x, g.Rx[e, f, i, j].m(phi_pq), mods)
    #end = time.time()
    #f2sym = end - start

    ## evaluation + sum
    #start = time.time()
    #val2 = 0
    #for x_val in np.arange(0, 5, 0.001):
    #    for t in range(1):
    #        for (i, j) in g.two2:
    #            val2 += myf[i, j, t](x_val)
    #end = time.time()
    #f2eval = end - start
    #
    #res2 = [('sym', f2sym), ('eval', f2eval)]

    # flow3: numerical integration
    # symbolic + sum
    print('starting flow 3')
    start = time.time()
    val3 = 0
    for x_val in np.arange(4, 5, 1):
        for t in range(1000):
            for (i, j) in g.two2:
                val3 += g.Rx[e, f, i, j].m_num(phi_pq, x_val)
    end = time.time()
    eval3 = end - start

    print('val1:\t', val1)
    # print('val2:\t', val2)
    print('val2:\t', val3)
    print('\nflow1:')
    total = 0
    for it in res1:
        total += it[1]
        print(it[0], ':\t', it[1])
    print('total:\t', total)
    #print('\nflow2:')
    #total = 0
    #for it in res2:
    #    total += it[1]
    #    print(it[0], ':\t', it[1])
    #print('total:\t', total)
    print('\nflow3:')
    print('total:\t', eval3)
