import basicInfo as bi
import Nformulas as fl
from networkx import is_connected
from enums import Stats
import time

# module test switch
switches = {
        Stats.MOMENT: 0,
        Stats.CMOMENT: 1
        }


# graph file
fpath = '../data/'
gname = 'g1'

# read graph
g = bi.readGraph(fpath, gname)

# connected test
print('is connected: ', is_connected(g))

# or independent pdfs phi_p and phi_q

phi_pq = 1
# phi_pq = 4 * p * q

# load basic info
bi.NloadInfo(g, phi_pq=phi_pq)

# moments stats
if switches[Stats.MOMENT] == 1:
    moment = fl.NMoment(g)
    # print(moment.eval(0))
    print(moment.eval(1))
    print(moment.eval(2))
    print(moment.eval(3))
    #print(moment.eval(1))
    #print(moment.eval(2))


# conditional moments stats
if switches[Stats.CMOMENT] == 1:
    ks = [0, 1, 2]
    es = [('1', '2')]
    ps = [0, 0.5, 1]
    cmoment = fl.NCMoment(g)

    for k in ks:
        for e in es:
            for p in ps:
                print(cmoment.eval(k, e, p))



