/*
 * green.h
 */

#ifndef INCLUDED_green_h
#define INCLUDED_green_h

typedef struct
{
  size_t nmax;     /* maximum spherical harmonic degree */
  size_t nnm;      /* number of total Green's functions for nmax */
  double R;        /* reference radius (km) */
  double *cosmphi; /* array of cos(m phi) values */
  double *sinmphi; /* array of sin(m phi) values */
  double *Plm;     /* associated Legendre functions */
  double *dPlm;    /* derivatives of associated Legendre functions */
} green_workspace;

/*
 * Prototypes
 */

green_workspace *green_alloc(const size_t nmax);
void green_free(green_workspace *w);
int green_calc_int(const double r, const double theta, const double phi,
                   double *X, double *Y, double *Z, green_workspace *w);
int green_calc_ext(const double r, const double theta, const double phi,
                   double *X, double *Y, double *Z, green_workspace *w);
size_t green_nmidx(const size_t n, const int m);

#endif /* INCLUDED_green_h */
