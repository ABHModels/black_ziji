#include <math.h>
#include <omp.h>
#include <stdlib.h>

void findMaxMin(const double *arr, int N, double *max, double *min) {
  if (N <= 0)
    return; // Handle empty array case

  *max = *min = arr[0]; // Initialize max and min with the first element
  for (int i = 1; i < N; i++) {
    if (arr[i] > *max) {
      *max = arr[i]; // Update max if a larger element is found
    }
    if (arr[i] < *min) {
      *min = arr[i]; // Update min if a smaller element is found
    }
  }
}

int binary_search(const double *arr, int n, double val) {

  if (n <= 1)
    return -1;

  int klo = 0;
  int khi = n - 1;
  int k;
  while ((khi - klo) > 1) {
    k = (khi + klo) / 2;
    if (arr[k] > val) {
      khi = k;
    } else {
      klo = k;
    }
  }
  return klo;
}

void RayLine(double *radi_em, double *gred, double *dSsc, int Nph, double alp,
             double E_line, double *Eobs, double *spec_out, int Nspec,
             int Npar) {

  double Intens;
  double gmax, gmin;
  int imin, imax;

  double Eline = E_line;

  for (int j = 0; j < Nspec; j++) {
    spec_out[j] = 0.;
  }

  findMaxMin(gred, Nph, &gmax, &gmin);
  imin = binary_search(Eobs, Nspec, gmin * Eline);
  imax = binary_search(Eobs, Nspec, gmax * Eline) + 1;

  int i, j;
  omp_set_num_threads(Npar);
#pragma omp parallel private(i, j, Intens)

#pragma omp for reduction(+ : spec_out[ : Nspec])
  for (i = 0; i < Nph; i++) {
    Intens = gred[i] * gred[i] * gred[i] * pow(radi_em[i], -alp) * dSsc[i];

    for (j = imin; j < imax; j++) {
      if ((gred[i] * Eline <= Eobs[j + 1]) && (gred[i] * Eline >= Eobs[j])) {
        spec_out[j] = spec_out[j] + Intens;
      }
    }
  }
}

void RayConv(double *radi_em, double *gred, double *dSsc, int Nph, double alp,
             double *Eobs, double *phot_ref, double *phot_out, int Nspec,
             int Npar) {

  double Trans;
  double logEmin, dlogE, Eh, El, dN_l, dN_h; // Eemis,logEmax,
  // int imin, imax;

  for (int j = 0; j < Nspec; j++) {
    phot_out[j] = 0.;
    phot_ref[j] /= (Eobs[j + 1] - Eobs[j]); // photon dN to number flux dN/dE
  }

  logEmin = log(Eobs[0]);
  // logEmax = log(Eobs[Nspec]);
  dlogE = log(Eobs[1]) - logEmin;

  int i, j, kem;
  omp_set_num_threads(Npar);
#pragma omp parallel private(i, j, kem, Eh, El, Trans, dN_l, dN_h)

#pragma omp for reduction(+ : phot_out[ : Nspec])
  for (i = 0; i < Nph; i++) {
    Trans = gred[i] * gred[i] * pow(radi_em[i], -alp) * dSsc[i];

    for (j = 0; j < Nspec; j++) {
      Eh = Eobs[j + 1] * gred[i];
      El = Eobs[j] * gred[i];

      if ((El < Eobs[Nspec]) && (Eh > Eobs[0])) {
        kem = (log(Eh) - logEmin) / dlogE;

        if (kem == 0) {
          dN_l = 0.;
          dN_h = Trans * phot_ref[j] * (Eh - Eobs[0]);

          phot_out[0] = phot_out[0] + dN_h;
        }

        else if (kem == Nspec) {
          dN_l = Trans * phot_ref[j] * (Eobs[Nspec] - El);
          dN_h = 0.;

          phot_out[Nspec - 1] = phot_out[Nspec - 1] + dN_l;
        }

        else {
          dN_l = Trans * phot_ref[j] * (Eobs[kem] - El);
          dN_h = Trans * phot_ref[j] * (Eh - Eobs[kem]);

          phot_out[kem - 1] = phot_out[kem - 1] + dN_l;
          phot_out[kem] = phot_out[kem] + dN_h;
        }
      }
    }
  }

  // for(int j=0; j<Nspec; j++)
  // {
  // 	phot_out[j] *= (Eobs[j+1] - Eobs[j]);//photon dN to number flux dN/dE
  // }
}

// phot_ref_2d  -- 100 * Nspec
void RayXill(double *radi_em, double *gred, double *cosem, double *dSsc,
             int Nph, double alp, double *Eobs, double *phot_ref_2d,
             double *phot_out, int Nspec, int Npar) {

  double incl_xill[] = {18.19487, 31.78833, 41.40962, 49.4584, 56.63298, 
                      63.25631, 69.51268, 75.52248, 81.37307, 87.13402};

  int N_incem_zone = 10;
  double Trans;
  double logEmin, dlogE, Eh, El, dN_l, dN_h; // Eemis,logEmax,
  double frac_inc, spec_inc, emis_ang;
  int linc;
  // int imin, imax;

  for (int j = 0; j < Nspec; j++) {
    phot_out[j] = 0.;
  }

  for (int i_em = 0; i_em < N_incem_zone; i_em++) {
    for (int j_en = 0; j_en < Nspec; j_en++) {
      phot_ref_2d[i_em * Nspec + j_en] /=
          (Eobs[j_en + 1] - Eobs[j_en]); // photon dN to number flux dN/dE
    }
  }

  logEmin = log(Eobs[0]);
  // logEmax = log(Eobs[Nspec]);
  dlogE = log(Eobs[1]) - logEmin;

  int i, j, kem;
  omp_set_num_threads(Npar);
#pragma omp parallel private(i, j, kem, Eh, El, Trans, dN_l, dN_h, linc,       \
                                 frac_inc, spec_inc, emis_ang)

#pragma omp for reduction(+ : phot_out[ : Nspec])
  for (i = 0; i < Nph; i++) {
    Trans = gred[i] * gred[i] * pow(radi_em[i], -alp) * dSsc[i];
    // linc = cosem[i] * (N_incem_zone - 1);
    // frac_inc = cosem[i] * (N_incem_zone - 1) - 1.0 * linc;
    emis_ang = acos(cosem[i])*180/M_PI;
    linc = binary_search(incl_xill, 10, emis_ang);
    frac_inc = (emis_ang - incl_xill[linc])/(incl_xill[linc+1] - incl_xill[linc]);


    for (j = 0; j < Nspec; j++) {
      if(emis_ang<incl_xill[0])
        spec_inc = phot_ref_2d[j];
      else if (emis_ang>incl_xill[9])
        spec_inc = phot_ref_2d[9 * Nspec + j];
      else
        spec_inc = (1. - frac_inc) * phot_ref_2d[linc * Nspec + j] +
                 frac_inc * phot_ref_2d[(linc + 1) * Nspec + j];

      // spec_inc = phot_ref_2d[linc * Nspec + j];

      Eh = Eobs[j + 1] * gred[i];
      El = Eobs[j] * gred[i];

      if ((El < Eobs[Nspec]) && (Eh > Eobs[0])) {
        kem = (log(Eh) - logEmin) / dlogE;

        if (kem == 0) {
          dN_l = 0.;
          dN_h = Trans * spec_inc * (Eh - Eobs[0]);

          phot_out[0] = phot_out[0] + dN_h;
        }

        else if (kem == Nspec) {
          dN_l = Trans * spec_inc * (Eobs[Nspec] - El);
          dN_h = 0.;

          phot_out[Nspec - 1] = phot_out[Nspec - 1] + dN_l;
        }

        else {
          dN_l = Trans * spec_inc * (Eobs[kem] - El);
          dN_h = Trans * spec_inc * (Eh - Eobs[kem]);

          phot_out[kem - 1] = phot_out[kem - 1] + dN_l;
          phot_out[kem] = phot_out[kem] + dN_h;
        }
      }
    }
  }

  // for(int j=0; j<Nspec; j++)
  // {
  // 	phot_out[j] *= (Eobs[j+1] - Eobs[j]);//photon dN to number flux dN/dE
  // }
}