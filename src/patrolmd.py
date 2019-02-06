from gurobipy import Model, GRB, quicksum
import randist as rt

# all inputs
gname = 'patrol01'
fpath = './data/'

W_dict = {
        'w1': [('1', '2'), ('2', '3'), ('3', '4'), ('4', '1')],
        'w2': [('1', '2'), ('2', '4'), ('4', '3'), ('3', '1')],
        'w3': [('2', '3'), ('3', '4'), ('4', '2')],
        'w4': [('1', '3'), ('3', '4'), ('4', '1')],
        'w5': [('1', '2'), ('2', '1')]
        }

coeff_type = rt.Stats.CDF
coeff_param = {'x_val': 3}
sense = GRB.MAXIMIZE

#coeff_type = rt.Stats.MOMENT
#coeff_param = {'k_val': 1}
#sense = GRB.MINIMIZE

# assume uniform
phi = rt.Phi('uniform', phi_pq=1)

# get coeff
fls = rt.Formulas(gname, phi, memorize=True)
fl = fls.get_formula(coeff_type)
c = fl.X_coeff(**coeff_param)

# some functions
def same_e(e, f):
    if e == f or (e[0] == f[1] and e[1] == f[0]):
        return True
    else:
        return False

# calculate a_w^e: the percentage of length of e
# in the closed walk w
def get_a(w, mye):
    tot_w = 0.0
    tot_e = 0.0
    for e in w:
        tmpl = float(g.edges[e]['l'])
        tot_w += tmpl
        if same_e(e, mye):
            tot_e += tmpl
    return tot_e / tot_w

# read graph
g = rt.readGraph(fpath, gname)

# basic sets
W = W_dict.keys()
E = g.edges()

# get the parameter a
a = {}

for w in W:
    for e in E:
        a[w, e] = get_a(W_dict[w], e)


# building the model
md = Model('patrol')

x = {}
for e in E:
    x[e] = md.addVar(vtype=GRB.CONTINUOUS, obj=c[e])

y = {}
for w in W:
    y[w] = md.addVar(vtype=GRB.CONTINUOUS)

md.modelSense = sense

md.update()

# constraints
md.addConstr(quicksum(y[w] for w in W) == 1)

for e in E:
    md.addConstr(quicksum(y[w] * a[w, e] for w in W) == x[e])

# solve
md.optimize()

# print solutions
if md.status == GRB.status.OPTIMAL:
    print('Obj: %g' % md.objVal)
    for w in W:
        print('y[%s]=%g' % (w, y[w].x))
    for e in E:
        print('x[%s]=%g' % (e, x[e].x))

