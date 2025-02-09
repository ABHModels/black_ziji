
import numpy as np
import matplotlib.pyplot as plt

from conv_core import RayXill, RayLine
from zijiray import LampostGeom


def LogSpace(start, stop, num_points):
    start_power = np.log(start)
    stop_power = np.log(stop)
    values = np.logspace(start_power, stop_power, num_points, base=np.e)
    return values


def PhotonToFlux(energy, phot):
    elow = energy[:-1]
    ehigh = energy[1:]
    En = (elow + ehigh) / 2
    dEn = ehigh - elow

    Flux_E = phot * En / dEn

    return En, Flux_E

def LogPlot(energy, photon, ax, col):
    eline = 6.4
    En, flux = PhotonToFlux(energy, photon)

    E_F_E = En*flux

    ax.loglog(En, E_F_E / E_F_E[En > eline][0], color=col)





N_BIN =  1000
energy = LogSpace(0.001, 1000.0, N_BIN+1)


inc = 70
spin = 0.998
alpha= 3.
gamma = 2.
logXi = 3.1



ray_xill = RayXill(gamma=gamma, alpha=alpha, Ecut=300.,
                   logXi=logXi, spin=spin, incl=inc, xtable_path="xillver-a-Ec5.fits")

phot_rayxill = ray_xill.Run(8, energy, cache=1)

line = RayLine(alpha=alpha, Eline=6.4, spin=spin, incl=inc)

phot_line = line.Run(8, energy, cache=1)


LampostGeom(spin=0.9, hlp=5)


fig, ax = plt.subplots()



LogPlot(energy, phot_rayxill, ax, "r")
# LogPlot(energy, phot_line, ax, "r")



ax.set_xlim(1, 50)

plt.show()
