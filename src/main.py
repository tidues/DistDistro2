import basicInfo as bi
import formulas as fl
import simulation as sl
from scipy import integrate
from sympy.abc import x, p, q
from sympy import lambdify
import numericFuncs as nf
from networkx import is_connected
from enums import Stats

# module test switch
switches = {
        Stats.MOMENT: 0,
        Stats.CDF: 0,
        Stats.PDF: 0,
        Stats.CMOMENT: 0,
        Stats.CCDF: 1,
        Stats.CPDF: 0,
        Stats.SIMULATION: 0
        }

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
bi.loadInfo(g, phi_pq=phi_pq, rational=True)

# bi.g_info_show(g)
# bi.check_R_measure_one(g, verb=False)
# bi.print_R(g)

# test values
valLst = [-0.1, -0.01, -0.001, 0, 1, 2, 3, -1, 9, 9.5, 9.9, 10]

# simulation 
if switches[Stats.SIMULATION] == 1:
    sim = sl.events(g)
    print(sim.E_d(f=lambda x: x**3))

# moments stats
if switches[Stats.MOMENT] == 1:
    moment = fl.Moment(g)
    print(moment.eval(0))
    print(moment.eval(1))
    print(moment.eval(2))

# cdf stats
if switches[Stats.CDF] == 1:
    cdf = fl.CDF(g)
    cdf.formula()
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
    cmoment = fl.CMoment(g)
    print(cmoment.eval(1, ('1','2'), 0))
    print(cmoment.eval(1, ('1','2'), 0.5))
    print(cmoment.eval(1, ('1','2'), 0.9))
    print(cmoment.eval(2, ('1','2'), 0.5))
    print(cmoment.eval(2, ('1','2'), 0.8))
    print(cmoment.eval(2, ('1','2'), 0.9))
    cmoment.formula(1, ('1', '2'))
    cmoment.plot(1, ('1', '2'))
    cmoment.plot(2, ('1', '2'))

# conditional cdf
if switches[Stats.CCDF] == 1:
    e = ('1', '2')
    f = ('2', '3')
    print(g.p1[('1', '2'), ('2', '3')])
    print(g.p2[('1', '2'), ('2', '3')])
    print(g.q1[('1', '2'), ('2', '3')])
    print(g.q2[('1', '2'), ('2', '3')])
    for (i, j) in g.two2:
        print(g.d[e[i]][f[j]])
    ccdf = fl.CCDF(g)
    print(ccdf.eval(('1', '2'), 0, 0))
    print(ccdf.eval(('1', '2'), 0, 2))
    print(ccdf.eval(('1', '2'), 0, 5))
    print(ccdf.formula(('1', '2')))
    ccdf.plot(('1', '2'), 0)
    print(ccdf.eval(('1', '2'), 0.5, 0))
    print(ccdf.eval(('1', '2'), 0.5, 2))
    print(ccdf.eval(('1', '2'), 0.5, 5))
    ccdf.plot(('1', '2'), 0.5)
    print(ccdf.eval(('1', '2'), 1, 0))
    print(ccdf.eval(('1', '2'), 1, 2))
    print(ccdf.eval(('1', '2'), 1, 5))
    ccdf.plot(('1', '2'), 1)

# conditional pdf
if switches[Stats.CPDF] == 1:
    print(g.p1[('1', '2'), ('2', '3')])
    print(g.p2[('1', '2'), ('2', '3')])
    print(g.q1[('1', '2'), ('2', '3')])
    print(g.q2[('1', '2'), ('2', '3')])
    cpdf = fl.CPDF(g)
    print(cpdf.eval(('1', '2'), 0, 0))
    print(cpdf.eval(('1', '2'), 0, 2))
    print(cpdf.eval(('1', '2'), 0, 5))
    cpdf.plot(('1', '2'), 0)
    print(cpdf.eval(('1', '2'), 0.5, 0))
    print(cpdf.eval(('1', '2'), 0.5, 2))
    print(cpdf.eval(('1', '2'), 0.5, 5))
    cpdf.plot(('1', '2'), 0.5)
    print(cpdf.eval(('1', '2'), 1, 0))
    print(cpdf.eval(('1', '2'), 1, 2))
    print(cpdf.eval(('1', '2'), 1, 5))
    cpdf.plot(('1', '2'), 1)




