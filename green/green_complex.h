/*
 * green_complex.h
 */

#ifndef INCLUDED_green_complex_h
#define INCLUDED_green_complex_h

#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

typedef struct
{
  size_t nmax;          /* maximum spherical harmonic degree */
  size_t mmax;          /* maximum spherical harmonic order */
  size_t nnm;           /* number of total Green's functions for nmax */
  double R;             /* reference radius (km) */
  double *Pnm;          /* associated Legendre functions */
  double *dPnm;         /* derivatives of associated Legendre functions */
  complex double *Ynm;  /* Pnm * exp(i m phi) */
  complex double *dYnm; /* dPnm * exp(i m phi) */
} green_complex_workspace;

/*
 * Prototypes
 */

green_complex_workspace *green_complex_alloc(const size_t nmax, const size_t mmax, const double R);
void green_complex_free(green_complex_workspace *w);
int green_complex_int(const double r, const double theta, const double phi,
                      complex double *X, complex double *Y, complex double *Z,
                      green_complex_workspace *w);
int green_complex_ext(const double r, const double theta, const double phi,
                      complex double *X, complex double *Y, complex double *Z,
                      green_complex_workspace *w);
int green_complex_calc_int(const size_t nmax, const size_t mmax, const double R,
                           const double r, const double theta,
                           const complex double *Ynm, const complex double *dYnm,
                           complex double *X, complex double *Y, complex double *Z);
int green_complex_calc_ext(const size_t nmax, const size_t mmax, const double R,
                           const double r, const double theta,
                           const complex double *Ynm, const complex double *dYnm,
                           complex double *X, complex double *Y, complex double *Z);
inline size_t green_complex_nmidx(const size_t n, const int m, const green_complex_workspace *w);
size_t green_complex_nnm(const green_complex_workspace *w);
int green_complex_Ynm(const double theta, const double phi, green_complex_workspace *w);

#endif /* INCLUDED_green_complex_h */
