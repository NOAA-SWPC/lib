/*
 * mfield_eval.h
 */

#ifndef INCLUDED_mfield_eval_h
#define INCLUDED_mfield_eval_h

typedef struct
{
  size_t nmax;      /* maximum spherical harmonic degree */
  double R;         /* geomagnetic reference radius (km) */
  double epoch;     /* t0 epoch (years) */

  size_t nnm;       /* number of (n,m) pairs */
  size_t p;         /* total number of model coefficients */

  double *c;        /* model coefficients */
  size_t sv_offset; /* offset into 'c' of SV coefficients */
  size_t sa_offset; /* offset into 'c' of SA coefficients */

  double *cosmphi;  /* array of cos(m phi) */
  double *sinmphi;  /* array of sin(m phi) */

  double *Plm;      /* associated legendres */
  double *dPlm;     /* associated legendre derivatives */

  double *dX;       /* basis functions for X */
  double *dY;       /* basis functions for Y */
  double *dZ;       /* basis functions for Z */
} mfield_eval_workspace;

/*
 * Prototypes
 */

mfield_eval_workspace *mfield_eval_alloc(const size_t nmax, const double epoch);
void mfield_eval_free(mfield_eval_workspace *w);
int mfield_eval(const double t, const double r, const double theta,
                const double phi, double B[4], mfield_eval_workspace *w);
int mfield_eval_ext(const double t, const double r, const double theta,
                    const double phi, const double E_st, const double I_st,
                    double B[4], mfield_eval_workspace *w);
int mfield_eval_dBdt(const double t, const double r, const double theta,
                     const double phi, double dBdt[4],
                     mfield_eval_workspace *w);
int mfield_eval_g_ext(const double t, const double r, const double theta, const double phi,
                      const double E_st, const double I_st,
                      const double *g, const double *dg,
                      double B[4], mfield_eval_workspace *w);
mfield_eval_workspace *mfield_eval_read(const char *filename);
int mfield_eval_green(const double r, const double theta, const double phi,
                      mfield_eval_workspace *w);
size_t mfield_eval_nmidx(const size_t n, const int m);

#endif /* INCLUDED_mfield_eval_h */
