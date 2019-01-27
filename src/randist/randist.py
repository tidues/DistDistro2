from .modules import formulas as fl
from .modules import Nformulas as Nfl
from .modules.enums import Stats
from .modules import basicInfo as bi
from .modules.basicInfo import readGraph
from sympy import lambdify
from sympy.abc import p, q
from .modules.numericFuncs import pdf_check
from .modules.commonFuncs import load_formulas, get_largest_component
from .modules.phis import Phi

class Formulas:
    # symbolic: use symbolic or numerical method, None is same as auto
    def __init__(self, gname, phi, fpath='../data/', rational=False):
        # read graph
        self.g = bi.readGraph(fpath, gname)
        self.g.symInfo = False
        # generate basic info
        bi.basicInfo(self.g, phi, rational)
        # dispatch map
        self.flS = {
                Stats.MOMENT: fl.Moment,
                Stats.CDF: fl.CDF,
                Stats.PDF: fl.PDF,
                Stats.CMOMENT: fl.CMoment,
                Stats.CCDF: fl.CCDF,
                Stats.CPDF: fl.CPDF
                }
        self.flN = {
                Stats.MOMENT: Nfl.NMoment,
                Stats.CDF: Nfl.NCDF,
                Stats.PDF: None,
                Stats.CMOMENT: Nfl.NCMoment,
                Stats.CCDF: Nfl.NCCDF,
                Stats.CPDF: None
                }

    # get different stats
    def get_formula(self, stats, symbolic=None):
        g = self.g
        if Stats.is_member(stats) is False:
            raise Exception('stats must be value from the enum Stats')

        # auto treat symbolic
        if symbolic is None:
            if stats == Stats.MOMENT or stats == Stats.CMOMENT:
                symbolic = False
            else:
                symbolic = True

        # dispatch to generate formulas
        if symbolic:
            # update phis info
            g.phi_p = g.phi.phi_p_S
            g.phi_q = g.phi.phi_q_S
            g.phi_pq = g.phi.phi_pq_S
            g.phi_qcp = g.phi.phi_qcp_S

            # add extra symbolic infos if necessary
            if g.symInfo is False:
                bi.symInfo(g)

            # dispatch
            return self.flS[stats](g)
        else:
            # update phis info
            g.phi_p = g.phi.phi_p_N
            g.phi_q = g.phi.phi_q_N
            g.phi_pq = g.phi.phi_pq_N
            g.phi_qcp = g.phi.phi_qcp_N

            # dispatch
            return self.flN[stats](g)


# info_dict keys and values
# moment: ks; cdf: nothing; pdf: nothing;
# cmoment: (ks, locs); ccdf: locs; cpdf: locs
def data_collector(gname, phi, mmtp=False, cdfp=False, pdfp=False, cmmtp=False, ccdfp=False, cpdfp=False):
    fls = Formulas(gname, phi)

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
        cpdf = fls.get_formula(Stats.CPDF)
        for loc in cpdfp:
            cpdf.plot(*loc)
        cpdf.save_formulas()

