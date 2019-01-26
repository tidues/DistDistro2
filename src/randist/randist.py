from .modules import formulas as fl
from .modules import Nformulas as Nfl
from .modules.enums import Stats
from .modules import basicInfo as bi
from sympy import lambdify
from sympy.abc import p, q
from .modules.numericFuncs import pdf_check
from .modules.commonFuncs import load_formulas

class Formulas:
    # symbolic: use symbolic or numerical method, None is same as auto
    def __init__(self, gname, fpath='../data/', phi_p=None, phi_q=None, phi_pq=None, rational=False):
        # read graph
        self.g = bi.readGraph(fpath, gname)
        self.g.symInfo = False
        # generate basic info
        bi.basicInfo(self.g, phi_p, phi_q, phi_pq, rational)
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
            g.phi_p = g.phi_p_S
            g.phi_q = g.phi_q_S
            g.phi_pq = g.phi_pq_S
            g.phi_qcp = g.phi_qcp_S

            # add extra symbolic infos if necessary
            if g.symInfo is False:
                bi.symInfo(g)

            # dispatch
            return self.flS[stats](g)
        else:
            # update phis info
            g.phi_p = lambdify(p, g.phi_p_S)
            g.phi_q = lambdify(q, g.phi_q_S)
            g.phi_pq = lambdify((q, p), g.phi_pq_S)
            g.phi_qcp = lambdify((q, p), g.phi_qcp_S)

            # dispatch
            return self.flN[stats](g)




