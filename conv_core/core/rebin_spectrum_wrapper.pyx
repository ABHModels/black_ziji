# rebin_spectrum_wrapper.pyx

cimport numpy as np
import numpy as np
ctypedef np.float64_t dtype_t




cdef extern from "rebin_spec.hpp":
    void rebin_spectrum(const double* ener, double* flu, int nbins, const double* ener0, const double* flu0, int nbins0)

def rebin_spectrum_py(np.ndarray[dtype_t, ndim=1] ener, 
                      np.ndarray[dtype_t, ndim=1] ener0,
                      np.ndarray[dtype_t, ndim=1] flu0):
    

    cdef int nbins = ener.shape[0] - 1
    cdef int nbins0 = ener0.shape[0] - 1
    cdef np.ndarray[dtype_t, ndim=1] flu = np.zeros(nbins)
    
    rebin_spectrum(&ener[0], &flu[0], nbins, &ener0[0], &flu0[0], nbins0)

    return flu
 
