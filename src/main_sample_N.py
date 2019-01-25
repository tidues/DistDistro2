import basicInfo as bi
import Nformulas as fl
from enums import Stats
# from commonFuncs import gcheck

# module test switch
switches = {
        Stats.MOMENT: 1,
        Stats.CMOMENT: 1
        }


# graph file
fpath = '../data/'
gname = 'planar_side_5'

# read graph
g = bi.readGraph(fpath, gname)

phi_pq = 1
# phi_pq = 4 * p * q

# load basic info
bi.NloadInfo(g, phi_pq=phi_pq)

## check all status
#check_res = gcheck(g)
#if check_res['total'] is not True:
#    print(check_res)

# moments stats
if switches[Stats.MOMENT] == 1:
    moment = fl.NMoment(g)
    m0 = moment.eval(0)
    print(m0)
    m1 = moment.eval(1)
    print(m1)
    m2 = moment.eval(2)
    print(m2)
    m3 = moment.eval(3)
    print(m3)
    print(m2 - m1 ** 2)
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



