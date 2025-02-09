import sys
import numpy as np

from bin.conv.ray_line import Ray_py, RayConv_py, RayXill_py
from bin.conv.rebin_spectrum_wrapper import rebin_spectrum_py as RebinSpec

from zijiray import TransitRay
from xillver import Xillver


def FluxToPhot(ener, flux):
    el = ener[:-1]
    eh = ener[1:]
    de = eh - el
    En_bin = (eh + el) / 2

    phot = flux * de / En_bin

    return phot


def PhotToFlux(ener, phot):
    el = ener[:-1]
    eh = ener[1:]
    de = eh - el
    En_bin = (eh + el) / 2

    return En_bin, phot * En_bin / de






class RayLine:
    def __init__(self, alpha=3., Eline=6.4, spin=0.998, incl = 60.) -> None:
        self.alpha = alpha
        self.Eline = Eline
        self.spin = spin
        self.Incl = incl

        self.transit_data_dict = None

    def Run(self, Npar, energy, cache = 0):

        TransitRay(Npar= Npar, spin=self.spin, inc=self.Incl, cache=cache)
        fname = "data/data_cache/transit/transit_data_dic_a%.3f_inc_%d.npz" % (self.spin, self.Incl)
        transit = np.load(fname, allow_pickle=True)
        self.transit_data_dict = transit
        radi_2d = self.transit_data_dict["radi_emitters"]
        garr_2d = self.transit_data_dict["redshift"]
        dSsc_2d = self.transit_data_dict["area_pixel_screen"]
        rem = radi_2d.flatten()
        garr = garr_2d.flatten()
        dS = dSsc_2d.flatten()

        phot_out = Ray_py(rem, garr, dS, self.alpha, self.Eline, energy, Npar)

        return phot_out


class RayConv:
    def __init__(self, alpha=3., spin=0.998, incl = 60.) -> None:
        self.alpha = alpha
        self.spin = spin
        self.Incl = incl

        self.transit_data_dict = None


    def Run(self, Npar, energy, phot_reflec, cache = 0):
        TransitRay(Npar= Npar, spin=self.spin, inc=self.Incl, cache=cache)
        fname = "data/data_cache/transit/transit_data_dic_a%.3f_inc_%d.npz" % (self.spin, self.Incl)
        transit = np.load(fname, allow_pickle=True)
        self.transit_data_dict = transit

        phot_reflec_copy = np.copy(phot_reflec)
        radi_2d = self.transit_data_dict["radi_emitters"]
        garr_2d = self.transit_data_dict["redshift"]
        dSsc_2d = self.transit_data_dict["area_pixel_screen"]
        rem = radi_2d.flatten()
        garr = garr_2d.flatten()
        dS = dSsc_2d.flatten()

        phot_out = RayConv_py(rem, garr, dS, self.alpha, energy, phot_reflec_copy, Npar)

        return phot_out


class RayXill:
    def __init__(self, alpha=3., gamma = 2., Ecut = 300., logXi= 3., Afe = 1.,
                spin=0.998, incl = 60., xtable_path = "xillver-a-Ec5.fits") -> None:
        self.alpha = alpha
        self.gamma = gamma
        self.Ecut = Ecut
        self.logXi = logXi
        self.Afe = Afe
        self.spin = spin
        self.Incl = incl

        self.xtable_path = xtable_path

        self.transit_data_dict = None

    def Run(self, Npar, energy, cache = 0):

        TransitRay(Npar= Npar, spin=self.spin, inc=self.Incl, cache=cache)
        fname = "data/data_cache/transit/transit_data_dic_a%.3f_inc_%d.npz" % (self.spin, self.Incl)
        transit = np.load(fname, allow_pickle=True)
        self.transit_data_dict = transit

        xillver = Xillver(gamma=self.gamma, Ecut=self.Ecut, logXi=self.logXi, 
                          Afe=self.logXi, xtable_path=self.xtable_path)
        xillver.Get()

        xill_ener = xillver.xill_ener
        xill_spec10 = xillver.xill_spec10

        radi_2d = self.transit_data_dict["radi_emitters"]
        garr_2d = self.transit_data_dict["redshift"]
        cosem_2d = self.transit_data_dict["cosem"]
        dSsc_2d = self.transit_data_dict["area_pixel_screen"]
        rem = radi_2d.flatten()
        garr = garr_2d.flatten()
        cosem = cosem_2d.flatten()
        dS = dSsc_2d.flatten()
        phot_reflec = xill_spec10.flatten()

        phot_out = RayXill_py(
            rem, garr, cosem, dS, self.alpha, xill_ener, phot_reflec, Npar
        )

        return RebinSpec(energy, xill_ener, phot_out)
