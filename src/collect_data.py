import randist as rt
from sympy.abc import p, q
from Timer import Timer
from sympy import Piecewise, Max

# module test switch
switches = {
        'test': 0,
        'cycuni': 0,
        'cycnuni': 0,
        'clinuni': 0,
        'griduni': 0,
        'gridbeta': 0,
        'griddiri': 0,
        'manhattonCondi': 1,
        'manhattonUncondi': 0
        }

Stats = rt.Stats

if switches['test'] == 1:
    # collection params
    gname = 'g0'
    phi = rt.Phi('uniform', 1)
    ks = [1, 2, 3]
    loc1 = (('1', '2'), 0.2)
    loc2 = (('1', '3'), 0.5)
    loc3 = (('3', '4'), 0)
    locs = [loc1, loc2, loc3]

    xs = [3, 5, 8]
    locx1 = list(loc1) + [xs[1]]
    locx2 = list(loc2) + [xs[1]]
    locxs = [locx1, locx2]

    mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    cdfp = {'collect': True, 'symbolic': None}
    cdfe = {'collect': True, 'symbolic': None, 'valst': xs}
    pdfp = {'collect': True, 'symbolic': None}
    pdfe = {'collect': True, 'symbolic': None, 'valst': xs}
    cmmtp ={'collect': True, 'symbolic': None, 'valst': (ks, locs)}
    ccdfp = {'collect': True, 'symbolic': None, 'valst': locs}
    ccdfe = {'collect': True, 'symbolic': None, 'valst': locxs}
    cpdfp = {'collect': True, 'symbolic': False, 'valst': locs}
    cpdfe = {'collect': True, 'symbolic': False, 'valst': locxs}

    d_jit = False
    memorize = True

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe, d_jit=d_jit, memorize=memorize)

if switches['cycuni'] == 1:
    # collection params
    gname = 'gcyc'
    phi = rt.Phi('uniform', 1)

    ks = [1, 2, 3]
    xs = [3, 4, 5]
    loc1 = (('1', '2'), 0.5)

    mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    cdfp = {'collect': True, 'symbolic': None}
    cdfe = False
    pdfp = {'collect': True, 'symbolic': None}
    pdfe = False
    cmmtp = False
    ccdfp = False
    ccdfe = False
    cpdfp = False
    cpdfe = False

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)

if switches['cycnuni'] == 1:
    # collection params
    gname = 'gcyc'
    phi_pq = 36 * p * q * (1-p) * (1-q)
    phi = rt.Phi('pqbeta', phi_pq)

    ks = [1, 2, 3]

    mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    cdfp = False
    cdfe = False
    pdfp = {'collect': True, 'symbolic': None}
    pdfe = False
    cmmtp = False
    ccdfp = False
    ccdfe = False
    cpdfp = False
    cpdfe = False

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)

if switches['clinuni'] == 1:
    # collection params
    gname = 'gcli'
    phi_pq = 6 * q * (1-q)
    phi = rt.Phi('qbeta', phi_pq)

    ks = [1, 2, 3]
    loc1 = (('1', '2'), 0.2)
    loc2 = (('1', '3'), 0.5)
    loc3 = (('3', '4'), 0)
    locs = [loc1, loc2, loc3]

    mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    cdfp = False
    cdfe = False
    pdfp = {'collect': True, 'symbolic': None}
    pdfe = False
    cmmtp = False
    ccdfp = False
    ccdfe = False
    cpdfp = {'collect': True, 'symbolic': None, 'valst': locs}
    cpdfe = False

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)

if switches['griduni'] == 1:
    # collection params
    gname = 'planar_side_10'
    # phi = rt.Phi('uniform', 1)

    phi_pq = lambda p, q: 1
    phi_p = lambda p: 1

    phi = rt.Phi('uniform', phi_pq, phi_p, phi_p, symbolic=False)

    ks = [1, 2, 3]
    xs = [3, 6, 9]
    # loc1 = (('1', '2'), 0)
    # loc2 = (('5', '6'), 0.5)
    # loc3 = (('45', '46'), 0.5)
    # locs = [loc1, loc2, loc3]
    loc2 = (('5', '6'), 0.5)
    locs = [loc2]
    locxs = []
    for x in xs:
        locx = tuple(list(loc2) + [x])
        locxs.append(locx)

    mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    cdfp = False
    cdfe = {'collect': True, 'symbolic': False, 'valst': xs}
    #pdfp = {'collect': True, 'symbolic': None}
    pdfp = False
    pdfe = {'collect': True, 'symbolic': False, 'valst': xs}
    cmmtp = {'collect': True, 'symbolic': None, 'valst': (ks, locs)}
    ccdfp = False
    ccdfe = {'collect': True, 'symbolic': False, 'valst': locxs}
    #cpdfp = {'collect': True, 'symbolic': None, 'valst': locs}
    cpdfp = False
    cpdfe = {'collect': True, 'symbolic': False, 'valst': locxs}

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)

if switches['gridbeta'] == 1:
    # collection params
    gname = 'planar_side_10'
    #phi_pq = 6 * q * (1-q)
    #phi = rt.Phi('qbeta', phi_pq)

    def phi_pq(p, q):
        return 6 * q * (1 - q)

    def phi_p(p):
        return 1

    def phi_q(q):
        return 6 * q * (1 - q)

    phi = rt.Phi('qbeta', phi_pq, phi_p, phi_q, symbolic=False)

    ks = [1, 2, 3]
    xs = [3, 6, 9]
    # loc1 = (('1', '2'), 0)
    # loc2 = (('5', '6'), 0.5)
    # loc3 = (('45', '46'), 0.5)
    # locs = [loc1, loc2, loc3]
    loc2 = (('5', '6'), 0.5)
    locs = [loc2]
    locxs = []
    for x in xs:
        locx = tuple(list(loc2) + [x])
        locxs.append(locx)

    mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    cdfp = False
    cdfe = {'collect': True, 'symbolic': False, 'valst': xs}
    #pdfp = {'collect': True, 'symbolic': None}
    pdfp = False
    pdfe = {'collect': True, 'symbolic': False, 'valst': xs}
    cmmtp = {'collect': True, 'symbolic': None, 'valst': (ks, locs)}
    ccdfp = False
    ccdfe = {'collect': True, 'symbolic': False, 'valst': locxs}
    #cpdfp = {'collect': True, 'symbolic': None, 'valst': locs}
    cpdfp = False
    cpdfe = {'collect': True, 'symbolic': False, 'valst': locxs}

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)

if switches['griddiri'] == 1:
    # collection params
    gname = 'planar_side_10'

    #phi_pq = Piecewise(
    #        (180 * p * q * (p + q - 1) ** 2, (p + q <= 1)),
    #        (180 * (p - 1) * (q - 1) * (p + q - 1) ** 2, (p + q > 1))
    #        )
    #phi_p = -15 * (-p + 4 * p ** 2 - 6 * p ** 3 + 3 * p ** 4)
    #phi_q = -15 * (-q + 4 * q ** 2 - 6 * q ** 3 + 3 * q ** 4)
    #phi = rt.Phi('dirichlet', phi_pq, phi_p, phi_q)

    def phi_pq(p, q):
        if p + q <= 1:
            return 180 * p * q * (p + q - 1) ** 2
        else:
            return 180 * (p - 1) * (q - 1) * (p + q - 1) ** 2

    def phi_p(p):
        return -15 * (-p + 4 * p ** 2 - 6 * p ** 3 + 3 * p ** 4)

    phi = rt.Phi('dirichlet', phi_pq, phi_p, phi_p, symbolic=False)

    ks = [1, 2, 3]
    xs = [3, 6, 9]
    # loc1 = (('1', '2'), 0)
    # loc2 = (('5', '6'), 0.5)
    # loc3 = (('45', '46'), 0.5)
    # locs = [loc1, loc2, loc3]
    loc2 = (('5', '6'), 0.5)
    locs = [loc2]
    locxs = []
    for x in xs:
        locx = tuple(list(loc2) + [x])
        locxs.append(locx)

    #mmtp = {'collect': True, 'symbolic': None, 'valst': ks}
    #cdfp = False
    #cdfe = {'collect': True, 'symbolic': False, 'valst': xs}
    ##pdfp = {'collect': True, 'symbolic': None}
    #pdfp = False
    #pdfe = {'collect': True, 'symbolic': False, 'valst': xs}
    #cmmtp = {'collect': True, 'symbolic': None, 'valst': (ks, locs)}
    #ccdfp = False
    #ccdfe = {'collect': True, 'symbolic': False, 'valst': locxs}
    ##cpdfp = {'collect': True, 'symbolic': None, 'valst': locs}
    #cpdfp = False
    #cpdfe = {'collect': True, 'symbolic': False, 'valst': locxs}

    mmtp = False
    cdfp = False
    cdfe = False
    #pdfp = {'collect': True, 'symbolic': False}
    pdfp = False
    pdfe = False
    cmmtp = False
    ccdfp = False
    ccdfe = False
    cpdfp = {'collect': True, 'symbolic': False, 'valst': locs}
    cpdfe = {'collect': True, 'symbolic': False, 'valst': locxs}

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)

if switches['manhattonCondi'] == 1:
    # collection params
    gname = 'manhattan_cc'
    phi = rt.Phi('uniform', 1)
#    def phi_pq(p, q):
#        if p + q <= 1:
#            return 180 * p * q * (p + q - 1) ** 2
#        else:
#            return 180 * (p - 1) * (q - 1) * (p + q - 1) ** 2
#
#    def phi_p(p):
#        return -15 * (-p + 4 * p ** 2 - 6 * p ** 3 + 3 * p ** 4)
#
#    phi = rt.Phi('dirichlet', phi_pq, phi_p, phi_p, symbolic=False)

    ks = [1, 2]
    x = 2000

    fire = (('346','370'), 0.353)
    bull = (('433', '468'), 0.675)
    police = (('807','866'), 0.568)
    empire = (('3369', '3419'), 0.384)
    tmsq = (('4116', '4168'), 0.773)
    bdy = (('4474', '4511'), 0.683)
    hosp = (('4564','4586'), 0.438)
    met = (('6019', '6065'), 0.328)

    locs = [fire, police, hosp]

    locxs = []

    for loc in locs:
        locx = tuple(list(loc) + [x])
        locxs.append(locx)

    mmtp = False
    cdfp = False
    cdfe = False
    pdfp = False
    pdfe = False
    #cmmtp = {'collect': True, 'symbolic': False, 'valst': (ks, locs)} 
    cmmtp = False
    ccdfp = False
    #ccdfe = {'collect': True, 'symbolic': False, 'valst': locxs}
    ccdfe = False
    cpdfp = {'collect': True, 'symbolic': False, 'valst': locs}
    #cpdfe = {'collect': True, 'symbolic': False, 'valst': locxs}
    cpdfe = False

    d_jit = True

    # collect data
    #rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp, d_jit=d_jit)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe, d_jit=d_jit)

if switches['manhattonUncondi'] == 1:
    # collection params
    gname = 'manhattan_cc'
    phi = rt.Phi('uniform', 1)

    ks = [1, 2]
    mmtp = {'collect': True, 'symbolic': False, 'valst': ks}
    cdfp = False
    cdfe = False
    pdfp = False
    pdfe = False
    cmmtp = False
    ccdfp = False
    ccdfe = False
    cpdfp = False
    cpdfe = False

    d_jit = False
    memorize = False

    # collect data
    # rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp, d_jit=d_jit, memorize=memorize)
    rt.data_collector(gname, phi, mmtp, cdfp, cdfe, pdfp, pdfe, cmmtp, ccdfp, ccdfe, cpdfp, cpdfe)
