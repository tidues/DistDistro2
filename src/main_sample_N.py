import randist as rt
# from commonFuncs import gcheck
Stats = rt.Stats

# module test switch
switches = {
        Stats.MOMENT: 1,
        Stats.CDF: 1,
        Stats.PDF: 0,
        Stats.CMOMENT: 1,
        Stats.CCDF: 1,
        Stats.CPDF: 0
        }


# graph file
gname = 'g3'

phi = rt.Phi('uniform', phi_pq=1)
# phi_pq = 4 * p * q

fls = rt.Formulas(gname, phi)

# moments stats
if switches[Stats.MOMENT] == 1:
    moment = fls.get_formula(Stats.MOMENT)
    vals = [0, 1, 2, 3]
    for v in vals:
        print(moment.eval(v))

# cdf stats
if switches[Stats.CDF] == 1:
    cdf = fls.get_formula(Stats.CDF, symbolic=False)
    vals = [0, 0.5, 1, 3, 5]
    for v in vals:
        print(cdf.eval(v))

# conditional moments stats
if switches[Stats.CMOMENT] == 1:
    cmoment = fls.get_formula(Stats.CMOMENT)
    ks = [0, 1, 2]
    es = [('1', '2')]
    ps = [0, 0.5, 1]

    for k in ks:
        for e in es:
            for p in ps:
                print(cmoment.eval(k, e, p))

# conditional cdf stats
if switches[Stats.CCDF] == 1:
    ccdf = fls.get_formula(Stats.CCDF, symbolic=False)
    es = [('1', '2')]
    ps = [0, 0.5, 1]
    xs = [0, 0.5, 1, 3, 5]

    for e in es:
        for p in ps:
            for x in xs:
                print(ccdf.eval(e, p, x))



