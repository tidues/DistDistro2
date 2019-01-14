import basicInfo as bi
import formulas as fl
import simulation as sl
from scipy import integrate
from sympy.abc import x, p, q
from sympy import lambdify
import numericFuncs as nf
from networkx import is_connected

# graph file
fpath = '../data/'
gname = 'g0'

# read graph
g = bi.readGraph(fpath, gname)

# connected test
print('is connected: ', is_connected(g))

# phi = 1
phi = 4 * p * q

# load basic info
bi.loadInfo(g, phi, rational=False)
print('d_max:', g.d_max)
# bi.g_info_show(g)
# bi.check_R_measure_one(g, verb=False)
# bi.print_R(g)

# generate formulas
moment = fl.Moment(g)
cdf = fl.CDF(g)
cdf.formula()
cdf.plot()
pdf = fl.PDF(g)
pdf.formula()
pdf.plot()

#mat_Dx = []
#for e in g.edges():
#    row_Dx = []
#    for f in g.edges():
#        row_Dx.append(pdf.mat_D[e, f].subs(x, 3))
#    mat_Dx.append(row_Dx)
#print(mat_Dx)

print(nf.pdf_check(pdf.f, (x, 0, 12)))

# simulation
#sim = sl.events(g)
#print(sim.E_d(f=lambda x: x**3))

# evaluate value
print('evaluating 0')
print(moment.eval(0))
print('evaluating 1')
print(moment.eval(1))
print('evaluating 2')
print(moment.eval(2))

valLst = [-0.1, -0.01, -0.001, 0, 1, 2, 3, -1, 9, 9.5, 9.9, 10]
#for x in valLst:
#    print(x, '\t', cdf.eval(x))

for x in valLst:
    print(x, '\t', pdf.eval(x))
