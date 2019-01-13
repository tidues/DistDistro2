import basicInfo as bi
import formula as fl
import simulation as sl
from scipy import integrate
from sympy.abc import x, p, q
from sympy import lambdify
import numericFuncs as nf

# graph file
fpath = '../data/'
gname = 'g3'

# read graph
print('reading graph')
g = bi.readGraph(fpath, gname)
phi = 4 * p * q

# load basic info
print('loading info')
bi.loadInfo(g, phi, rational=True)
# bi.g_info_show(g)
# bi.check_R_measure_one(g, verb=False)
# bi.print_R(g)

# generate formulas
print('generating moment function')
moment = fl.Moment(g)
print('generating cdf function')
cdf = fl.CDF(g)
print(cdf.f)
print('generating pdf function')
pdf = fl.PDF(g)
print(pdf.f)

mat_Dx = []
for e in g.edges():
    row_Dx = []
    for f in g.edges():
        row_Dx.append(pdf.mat_D[e, f].subs(x, 3))
    mat_Dx.append(row_Dx)
print(mat_Dx)

mods = ['numpy', {'theta': nf.theta, 'eta': nf.eta}]
f = lambdify(x, pdf.f, modules=mods)
print('f(0):', f(0))
print(integrate.quad(f, 0, 12))

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
print('evaluating 3')

valLst = [-0.1, -0.01, -0.001, 0, 1, 2, 3, -1, 9, 9.5, 9.9, 10]
#for x in valLst:
#    print(x, '\t', cdf.eval(x))

for x in valLst:
    print(x, '\t', pdf.eval(x))
