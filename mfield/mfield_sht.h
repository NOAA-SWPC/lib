/*
 * mfield_sht.h
 */

#ifndef INCLUDED_mfield_sht_h
#define INCLUDED_mfield_sht_h

#include <gsl/gsl_vector.h>

#include <nfft3.h>

typedef struct
{
  size_t N; /* maximum spherical harmonic degree */
  size_t M; /* number of data points for evaluation */

  nfsft_plan plan;
} mfield_sht_workspace;

mfield_sht_workspace *mfield_sht_alloc(const size_t N, const size_t M);
void mfield_sht_free(mfield_sht_workspace *w);
int mfield_sht_calc(const gsl_vector *g, const double r,
                    const gsl_vector *theta,
                    const gsl_vector *phi, gsl_vector *Y, gsl_vector *Z,
                    mfield_sht_workspace *w);

#endif /* INCLUDED_mfield_sht_h */
