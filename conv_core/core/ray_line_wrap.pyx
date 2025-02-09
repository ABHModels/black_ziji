# rebin_spectrum_wrapper.pyx

import numpy as np
cimport numpy as np
ctypedef np.float64_t dtype_t

cdef extern from "ray_line.h":
    void RayLine(double* radi_em, double* gred,double* dSsc, int Nph, double alp, double E_line, double *Eobs, double* spec_out, int Nspec, int Npar)
    void RayConv(double* radi_em, double* gred,double* dSsc, int Nph, 
			 double alp, 
			 double *Eobs, double* phot_ref, double* phot_out, int Nspec, 
			 int Npar)

    void RayXill(double *radi_em, double *gred, double *cosem, double *dSsc,
        int Nph, double alp, double *Eobs, double *phot_ref_2d,
        double *phot_out, int Nspec, int Npar)

def Ray_py(np.ndarray[dtype_t, ndim=1] radi_em, 
                      np.ndarray[dtype_t, ndim=1] gred, 
                      np.ndarray[dtype_t, ndim=1] dSsc,
						alp,
                                                Eline,
                      np.ndarray[dtype_t, ndim=1] Eobs,
                      Npar):


    cdef int Nph = radi_em.shape[0] - 1
    cdef int Nspec = Eobs.shape[0] - 1

    cdef np.ndarray[dtype_t, ndim=1] spec_out = np.zeros(Nspec)


    RayLine(<dtype_t*>radi_em.data, <dtype_t*>gred.data, <dtype_t*>dSsc.data, 
            Nph, alp, Eline, <dtype_t*>Eobs.data, <dtype_t*>spec_out.data, Nspec, Npar)

    return spec_out


def RayConv_py(np.ndarray[dtype_t, ndim=1] radi_em, 
                      np.ndarray[dtype_t, ndim=1] gred, 
                      np.ndarray[dtype_t, ndim=1] dSsc,
						alp,
                      np.ndarray[dtype_t, ndim=1] Eobs,
                      np.ndarray[dtype_t, ndim=1] phot_ref,
                      Npar):

    cdef int Nph = radi_em.shape[0] - 1
    cdef int Nspec = Eobs.shape[0] - 1
    cdef np.ndarray[dtype_t, ndim=1] phot_out = np.zeros(Nspec)

    RayConv(<dtype_t*>radi_em.data, <dtype_t*>gred.data, <dtype_t*>dSsc.data, 
            Nph, alp, <dtype_t*>Eobs.data, <dtype_t*>phot_ref.data, <dtype_t*>phot_out.data, Nspec, Npar)

    return phot_out



def RayXill_py(np.ndarray[dtype_t, ndim=1] radi_em, 
                      np.ndarray[dtype_t, ndim=1] gred,
                      np.ndarray[dtype_t, ndim=1] cosem, 
                      np.ndarray[dtype_t, ndim=1] dSsc,
						alp,
                      np.ndarray[dtype_t, ndim=1] Eobs,
                      np.ndarray[dtype_t, ndim=1] phot_ref_2d,
                      Npar):

    cdef int Nph = radi_em.shape[0] - 1
    cdef int Nspec = Eobs.shape[0] - 1
    cdef np.ndarray[dtype_t, ndim=1] phot_out = np.zeros(Nspec)

    RayXill(<dtype_t*>radi_em.data, <dtype_t*>gred.data, <dtype_t*>cosem.data, <dtype_t*>dSsc.data, Nph, alp, 
            <dtype_t*>Eobs.data, <dtype_t*>phot_ref_2d.data, <dtype_t*>phot_out.data, Nspec, Npar)

    return phot_out