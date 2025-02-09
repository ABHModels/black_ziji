import os
import sys
import numpy as np

# ----------------------------------------------------------------------------------------------------
# INPUT PARAMETERS
# ----------------------------------------------------------------------------------------------------

path_to_xillver_file = "xillver-a-Ec5.fits"
gamma = 1.8
afe = 1
logxi = 1.
ecut = 300
incl = 70

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------



xillver_command = "python3 xillver/get_xillver_spectrum.py %f %f %f %f %f %s"%(gamma, afe, logxi, ecut, incl, path_to_xillver_file)
os.system(xillver_command)


xill_spec_filename = "data/xill_spec_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_incl_%.2f.dat" % (gamma, afe, logxi, ecut, incl)

