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
        Stats.CCDF: 0,
        Stats.CPDF: 1,
        Stats.SIMULATION: 0
        }

# graph file
fpath = '../data/'
gname = 'g0'

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



