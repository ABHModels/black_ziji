
// nbins = n_en - 1 = n_flu
void rebin_spectrum(double *ener, double *flu, int nbins, double *ener0, double *flu0, int nbins0) {

  int ii;
  int jj;
  int imin = 0;
  int imax = 0;

  for (ii = 0; ii < nbins; ii++) {

    flu[ii] = 0.0;

    /* check of the bin is outside the given energy range */
    if ((ener0[0] <= ener[ii + 1]) && (ener0[nbins0] >= ener[ii])) {

      /* need to make sure we are in the correct bin */
      while (ener0[imin] <= ener[ii] && imin <= nbins0) {
        imin++;
      }
      // need to set it back, as we just crossed to the next bin
      if (imin > 0) {
        imin--;
      }
      while ((ener0[imax] <= ener[ii + 1] && imax < nbins0)) {
        imax++;
      }
      if (imax > 0) {
        imax--;
      }

      double elo = ener[ii];
      double ehi = ener[ii + 1];
      if (elo < ener0[imin]) elo = ener0[imin];
      if (ehi > ener0[imax + 1]) ehi = ener0[imax + 1];

      if (imax == imin) {
        flu[ii] = (ehi - elo) / (ener0[imin + 1] - ener0[imin]) * flu0[imin];
      } else {

        double dmin = (ener0[imin + 1] - elo) / (ener0[imin + 1] - ener0[imin]);
        double dmax = (ehi - ener0[imax]) / (ener0[imax + 1] - ener0[imax]);

        flu[ii] += flu0[imin] * dmin + flu0[imax] * dmax;

        for (jj = imin + 1; jj <= imax - 1; jj++) {
          flu[ii] += flu0[jj];
        }

      }

    }

  }
}