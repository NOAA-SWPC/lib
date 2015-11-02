/*
 * mfield_green.h
 */

#ifndef INCLUDED_mfield_green_h
#define INCLUDED_mfield_green_h

typedef struct
{
  size_t nmax;      /* maximum spherical harmonic degree */
  size_t nnm;       /* total number of (n,m) coefficients */
  double *cosmphi;  /* array of cos(m phi) */
  double *sinmphi;  /* array of sin(m phi) */
  double *Plm;      /* associated legendres */
  double *dPlm;     /* associated legendre derivatives */
  double *dX;       /* internal basis functions for X */
  double *dY;       /* internal basis functions for Y */
  double *dZ;       /* internal basis functions for Z */
  double *dX_ext;   /* external basis functions for X */
  double *dY_ext;   /* external basis functions for Y */
  double *dZ_ext;   /* external basis functions for Z */
  double R;         /* reference radius (km) */
} mfield_green_workspace;

/*
 * Prototypes
 */

mfield_green_workspace *mfield_green_alloc(const size_t nmax, const double R);
void mfield_green_free(mfield_green_workspace *w);
int mfield_green_calc(const double r, const double theta, const double phi,
                      mfield_green_workspace *w);
int mfield_green_ext(const double r, const double theta, const double phi,
                     mfield_green_workspace *w);
size_t mfield_green_nmidx(const size_t n, const int m);
double mfield_green_X(const size_t idx, const mfield_green_workspace *w);
double mfield_green_Y(const size_t idx, const mfield_green_workspace *w);
double mfield_green_Z(const size_t idx, const mfield_green_workspace *w);

#endif /* INCLUDED_mfield_green_h */
