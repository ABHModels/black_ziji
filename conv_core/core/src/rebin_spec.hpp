
#include <stddef.h> // for size_t

// rebin_spectrum:
//   ener[0..nbins]                : new bin edges
//   rebinned_photon_counts[0..nbins-1]
//                                 : output counts per new bin
//   ener0[0..nbins0]              : original bin edges
//   orig_photon_counts[0..nbins0-1]
//                                 : input counts per original bin
void rebin_spectrum(const double *restrict ener,
                    double *restrict rebinned_photon_counts, size_t nbins,
                    const double *restrict ener0,
                    const double *restrict orig_photon_counts, size_t nbins0) {
  size_t j = 0; // pointer into original bins

  for (size_t i = 0; i < nbins; ++i) {
    double elo = ener[i];
    double ehi = ener[i + 1];
    double sum = 0.0;

    // skip old bins that end before the new-bin low edge
    while (j < nbins0 && ener0[j + 1] <= elo) {
      ++j;
    }

    // accumulate overlaps with original bins k = j … until start ≥ ehi
    for (size_t k = j; k < nbins0 && ener0[k] < ehi; ++k) {
      // inline min/max
      double overlap_lo = (elo > ener0[k] ? elo : ener0[k]);
      double overlap_hi = (ehi < ener0[k + 1] ? ehi : ener0[k + 1]);

      if (overlap_hi > overlap_lo) {
        // compute bin-width on the fly
        double width0 = ener0[k + 1] - ener0[k];
        sum += (overlap_hi - overlap_lo) / width0 * orig_photon_counts[k];
      }
    }

    rebinned_photon_counts[i] = sum;
  }
}
