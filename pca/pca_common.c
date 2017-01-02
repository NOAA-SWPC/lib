
/*
taper_knm()
  Taper high degree knm coefficients to try to reduce ringing.
We use a cosine taper set to 1 at nmin, and going to 0 at
nmax
*/

static int
taper_knm(const size_t nmin, gsl_matrix * K, green_workspace * green_p)
{
  const size_t nmax = green_p->nmax;
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, green_p->mmax);
      int m;
      double wn = 1.0;

      /* compute taper weight */
      if (n > nmin)
        wn = cos((n - nmin) * M_PI / (double)nmax);

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, green_p);
          gsl_vector_view row = gsl_matrix_row(K, cidx);

          gsl_vector_scale(&row.vector, wn);
        }
    }

  return GSL_SUCCESS;
}
