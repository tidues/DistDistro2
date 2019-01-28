import randist as rt
from sympy.abc import p, q

# module test switch
switches = {
        'test': 1,
        'cycuni': 0,
        'cycnuni': 0,
        'clinuni': 0,
        'griduni': 0,
        'gridnuni': 0,
        'manhattonCondi': 0,
        'manhattonUncondi': 0
        }

Stats = rt.Stats

if switches['test'] == 1:
    # collection params
    #gname = 'g0'
    #phi = rt.Phi('uniform', 1)
    #mmtp = [1, 2, 3]
    #cdfp = True
    #pdfp = True
    #loc1 = (('1', '2'), 0.2)
    #loc2 = (('1', '3'), 0.5)
    #loc3 = (('3', '4'), 0)
    #locs = [loc1, loc2, loc3]
    #cmmtp = (mmtp, locs)
    #ccdfp = locs
    #cpdfp = locs
    #condi_e = ('1', '3')
    gname = 'g3'
    phi = rt.Phi('uniform', 1)
    ks = [1, 2, 3]
    mmtp = False
    cdfp = False
    pdfp = False
    loc1 = (('1', '2'), 0.2)
    loc2 = (('1', '3'), 0.5)
    loc3 = (('2', '3'), 0)
    locs = [loc1, loc2, loc3]
    cmmtp = (ks, locs)
    ccdfp = locs
    cpdfp = locs
    d_jit = False
    #condi_e = ('1', '3')

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp, d_jit=d_jit)

if switches['cycuni'] == 1:
    # collection params
    gname = 'gcyc'
    phi = rt.Phi('uniform', 1)
    mmtp = [1, 2, 3]
    cdfp = True
    pdfp = True
    cmmtp = False
    ccdfp = False
    cpdfp = False

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)

if switches['cycnuni'] == 1:
    # collection params
    gname = 'gcyc'
    phi_pq = 36 * p * q * (1-p) * (1-q)
    phi = rt.Phi('pqbeta', phi_pq)
    mmtp = [1, 2, 3]
    cdfp = False
    pdfp = True
    cmmtp = False
    ccdfp = False
    cpdfp = False

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)

if switches['clinuni'] == 1:
    # collection params
    gname = 'gcli'
    phi_pq = 6 * q * (1-q)
    phi = rt.Phi('qbeta', phi_pq)
    mmtp = [1, 2, 3]
    cdfp = False
    pdfp = True
    cmmtp = False
    ccdfp = False
    loc1 = (('1', '2'), 0.2)
    loc2 = (('1', '3'), 0.5)
    loc3 = (('3', '4'), 0)
    cpdfp = [loc1, loc2, loc3]

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)

if switches['griduni'] == 1:
    # collection params
    gname = 'planar_side_10'
    phi = rt.Phi('uniform', 1)
    mmtp = [1, 2, 3]
    cdfp = False
    pdfp = True
    cmmtp = False
    ccdfp = False
    loc1 = (('1', '2'), 0)
    loc2 = (('5', '6'), 0.5)
    loc3 = (('45', '46'), 0.5)
    cpdfp = [loc1, loc2, loc3]

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)

if switches['gridnuni'] == 1:
    # collection params
    gname = 'planar_side_10'
    phi_pq = 6 * q * (1-q)
    phi = rt.Phi('qbeta', phi_pq)
    mmtp = [1, 2, 3]
    cdfp = False
    pdfp = True
    cmmtp = False
    ccdfp = False
    loc1 = (('1', '2'), 0)
    loc2 = (('5', '6'), 0.5)
    loc3 = (('45', '46'), 0.5)
    loc4 = (('45', '46'), 0.6)
    cpdfp = [loc1, loc2, loc3, loc4]

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)

if switches['manhattonCondi'] == 1:
    # collection params
    gname = 'manhattan_cc'
    phi = rt.Phi('uniform', 1)
    mmtp = False
    cdfp = False
    pdfp = False

    ks = [1, 2]

    bull = (('433', '468'), 0.675)
    empire = (('3369', '3419'), 0.384)
    tmsq = (('4116', '4168'), 0.773)
    bdy = (('4474', '4511'), 0.683)
    met = (('6019', '6065'), 0.328)
    locs = [bull, empire, tmsq, bdy, met]

    cmmtp = (ks, locs)
    ccdfp = False
    cpdfp = locs

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)

if switches['manhattonUncondi'] == 1:
    # collection params
    gname = 'manhattan_cc'
    phi = rt.Phi('uniform', 1)
    mmtp = [1, 2]
    cdfp = False
    pdfp = False
    cmmtp = False
    ccdfp = False
    cpdfp = False

    # collect data
    rt.data_collector(gname, phi, mmtp, cdfp, pdfp, cmmtp, ccdfp, cpdfp)
