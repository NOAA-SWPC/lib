/*
 * green.h
 */

#ifndef INCLUDED_green_h
#define INCLUDED_green_h

typedef struct
{
  size_t nmax;     /* maximum spherical harmonic degree */
  size_t nnm;      /* number of (n,m) coefficients */

  double *dX;      /* X basis functions for internal field */
  double *dY;      /* Y basis functions for internal field */
  double *dZ;      /* Z basis functions for internal field */

  double *dX_ext;  /* X basis functions for external field */
  double *dY_ext;  /* Y basis functions for external field */
  double *dZ_ext;  /* Z basis functions for external field */

  double *Plm;     /* legendre functions */
  double *dPlm;    /* legendre derivatives */

  double *cosmphi; /* cos(m phi) */
  double *sinmphi; /* sin(m phi) */
} green_workspace;

/*
 * Prototypes
 */

green_workspace *green_alloc(const size_t nmax);
void green_free(green_workspace *w);
size_t green_nnm(const green_workspace *w);
size_t green_nmax(const green_workspace *w);
int green_calc(const double r, const double theta, const double phi,
               const double R, green_workspace *w);
int green_eval(const double r, const double theta, const double phi,
               const double R, const double *coeffs,
               double B[3], green_workspace *w);
int green_calc_ext(const double r, const double theta, const double phi,
                   const double R, green_workspace *w);
size_t green_nmidx(const size_t n, const int m);

#endif /* INCLUDED_green_h */
