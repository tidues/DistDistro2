import basicInfo as bi
import Nformulas as fl
from enums import Stats
# from commonFuncs import gcheck

# module test switch
switches = {
        Stats.MOMENT: 0,
        Stats.CDF: 0,
        Stats.CMOMENT: 0,
        Stats.CCDF: 1
        }


# graph file
fpath = '../data/'
gname = 'g0'

# read graph
g = bi.readGraph(fpath, gname)

phi_pq = 1
# phi_pq = 4 * p * q

# load basic info
bi.NloadInfo(g, phi_pq=phi_pq)

# moments stats
if switches[Stats.MOMENT] == 1:
    moment = fl.NMoment(g)
    vals = [0, 1, 2, 3]
    for v in vals:
        print(moment.eval(v))

# cdf stats
if switches[Stats.CDF] == 1:
    cdf = fl.NCDF(g)
    vals = [0, 0.5, 1, 3, 5]
    for v in vals:
        print(cdf.eval(v))

# conditional moments stats
if switches[Stats.CMOMENT] == 1:
    cmoment = fl.NCMoment(g)
    ks = [0, 1, 2]
    es = [('1', '2')]
    ps = [0, 0.5, 1]

    for k in ks:
        for e in es:
            for p in ps:
                print(cmoment.eval(k, e, p))

# conditional cdf stats
if switches[Stats.CCDF] == 1:
    ccdf = fl.NCCDF(g)
    es = [('1', '2')]
    ps = [0, 0.5, 1]
    xs = [0, 0.5, 1, 3, 5]

    for e in es:
        for p in ps:
            for x in xs:
                print(ccdf.eval(e, p, x))



