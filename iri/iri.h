/*
 * iri.h
 * Patrick Alken
 */

#ifndef INCLUDED_iri_h
#define INCLUDED_iri_h

#include <indices/indices.h>

#define IRI_NR_OUTF         20
#define IRI_NC_OUTF         1000
#define IRI_N_OARR          100
#define IRI_N_FLAGS         50 /* JF flags */

#define IRI_FIDX(i, j)           ((j)*IRI_NR_OUTF + (i))

typedef struct
{
  double Ne;    /* electron density in m^{-3} */
  double Tn;    /* neutral temperature in K */
  double Ti;    /* ion temperature in K */
  double Te;    /* electron temperature in K */
  double n_Op;  /* O+ ion density in m^{-3} */
  double n_Hp;  /* H+ ion density in m^{-3} */
  double n_HEp; /* HE+ ion density in m^{-3} */
  double n_O2p; /* O2+ ion density in m^{-3} */
  double n_NOp; /* NO+ ion density in m^{-3} */
  double n_Np;  /* N+ ion density in m^{-3} */
  double n_clus; /* cluster ion density in m^{-3} */
} iri_result;

typedef struct
{
  /* output arrays */
  float outf[IRI_NR_OUTF * IRI_NC_OUTF];
  float oarr[IRI_N_OARR];

  int jf[IRI_N_FLAGS]; /* input flags */
  int jmag; /* 0 = geographic, 1 = geomagnetic coords */

  /* array containing IRI parameters for each altitude */
  iri_result *iri_results;

  size_t nalt;          /* size of iri_results array */

  double f107_override; /* override F10.7 value */

  f107_workspace *f107_workspace_p;
} iri_workspace;

/*
 * Prototypes
 */

iri_workspace *iri_alloc(size_t nalt, const char *f107_datadir);
void iri_free(iri_workspace *w);
void iri_f107_override(double f107, iri_workspace *w);
int iri_calc(double theta, double phi, time_t t, double altmin,
             double altstp, size_t nalt, iri_workspace *w);
iri_result *iri_get_result(size_t idx, iri_workspace *w);

/* IRI prototypes */

void initialize_(void);
void iri_sub_(int *jf, int *jmag, float *alati, float *along,
              int *iyyyy, int *mmdd, float *dhour, float *hbeg,
              float *hend, float *hstp, float *outf, float *oarr);

#endif /* INCLUDED_iri_h */
