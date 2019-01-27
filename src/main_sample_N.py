import randist as rt
# from commonFuncs import gcheck
Stats = rt.Stats

# module test switch
switches = {
        Stats.MOMENT: 0,
        Stats.CDF: 0,
        Stats.PDF: 1,
        Stats.CMOMENT: 0,
        Stats.CCDF: 0,
        Stats.CPDF: 1
        }


# graph file
gname = 'g2'

phi = rt.Phi('uniform', phi_pq=1)
# phi_pq = 4 * p * q

fls = rt.Formulas(gname, phi)

# test values
valLst = [-0.1, -0.01, -0.001, 0, 1, 2, 3, -1, 9, 9.5, 9.9, 10]

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

# pdf stats
if switches[Stats.PDF] == 1:
    pdf = fls.get_formula(Stats.PDF, symbolic=False)
    for x_val in valLst:
        print(x_val, '\t', pdf.eval(x_val))

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

# conditional pdf stats
if switches[Stats.CPDF] == 1:
    cpdf = fls.get_formula(Stats.CPDF, symbolic=False)
    e = ('1', '2')
    #ps = [0]
    ps = [0, 0.2, 0.4, 0.5, 0.8, 1]
    x = 3

    for p in ps:
        print(cpdf.eval(e, p, x))
