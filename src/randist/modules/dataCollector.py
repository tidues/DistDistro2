from .formulas import Formulas
from .enums import Stats

# info_dict keys and values
# moment: ks; cdf: nothing; pdf: nothing;
# cmoment: (ks, locs); ccdf: locs; cpdf: locs
def data_collector(gname, phi, mmtp=False, cdfp=False, pdfp=False, cmmtp=False, ccdfp=False, cpdfp=False, d_jit=False, memorize=True):
    fls = Formulas(gname, phi, d_jit=d_jit, memorize=memorize)

    # moments
    if mmtp is not False:
        moment = fls.get_formula(Stats.MOMENT)
        for k in mmtp:
            moment.eval(k)

    # cdf
    if cdfp is not False:
        cdf = fls.get_formula(Stats.CDF)
        cdf.plot()
        cdf.save_formulas()

    # pdf
    if pdfp is not False:
        pdf = fls.get_formula(Stats.PDF)
        pdf.plot()
        pdf.save_formulas()

    # cmoments
    if cmmtp is not False:
        cmoment = fls.get_formula(Stats.CMOMENT)
        for k in cmmtp[0]:
            for loc in cmmtp[1]:
                cmoment.eval(k, *loc)

    # ccdf
    if ccdfp is not False:
        ccdf = fls.get_formula(Stats.CCDF)
        for loc in ccdfp:
            ccdf.plot(*loc)
        ccdf.save_formulas()

    # cpdf
    if cpdfp is not False:
        cpdf = fls.get_formula(Stats.CPDF, symbolic=False)
        for loc in cpdfp:
            cpdf.plot(*loc)
        cpdf.save_formulas()

