
/*
taper_knm()
  Taper high degree knm coefficients to try to reduce ringing.
We use a cosine taper set to 1 at nmin, and going to 0 at
nmax + 1. Don't make it go to 0 at nmax or the spectrum will
suddenly drop many orders of magnitude between nmax-1 and nmax
*/

static int
taper_knm(const size_t nmin, gsl_matrix * K, green_workspace * green_p)
{
  const size_t nmax = green_p->nmax;
  const double fac = M_PI / (2.0 * (nmax + 1 - nmin));
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, green_p->mmax);
      int m;
      double wn = 1.0;

      /* compute taper weight */
      if (n > nmin)
        {
          double val = cos(fac * (n - nmin));
          wn = val * val;
        }

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, green_p);
          gsl_vector_view row = gsl_matrix_row(K, cidx);

          gsl_vector_scale(&row.vector, wn);
        }
    }

  return GSL_SUCCESS;
}
