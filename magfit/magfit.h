/*
 * magfit.h
 */

#ifndef INCLUDED_magfit_h
#define INCLUDED_magfit_h

#include <satdata/satdata.h>

#include "track.h"

/* mu_0 in units of: nT / (kA km^{-1}) */
#define MAGFIT_MU_0                       (400.0 * M_PI)

/* general flags */
#define MAGFIT_FLG_FIT_X                  (1 << 0) /* fit X component */
#define MAGFIT_FLG_FIT_Y                  (1 << 1) /* fit Y component */
#define MAGFIT_FLG_FIT_Z                  (1 << 2) /* fit Z component */

/* SECS flags */
#define MAGFIT_SECS_FLG_FIT_DF            (1 << 0) /* fit divergence-free SECS */
#define MAGFIT_SECS_FLG_FIT_CF            (1 << 1) /* fit curl-free SECS */

/* recommended lmax for 1D SECS */
#define MAGFIT_SECS_LMAX                  200

typedef struct
{
  /* SECS parameters */
  double lat_spacing1d; /* latitude pole spacing for 1D SECS (degrees) */
  double lat_spacing2d; /* latitude pole spacing for 2D SECS (degrees) */
  double lat_min;       /* minimum geocentric latitude for SECS poles (degrees) */
  double lat_max;       /* maximum geocentric latitude for SECS poles (degrees) */
  double lon_spacing;   /* longitude pole spacing for 2D SECS (degrees) */
  double lon_min;       /* minimum longitude for SECS poles (degrees) */
  double lon_max;       /* maximum longitude for SECS poles (degrees) */
  double R;             /* reference radius (km) */
  size_t lmax;          /* maximum spherical harmonic degree for 1D SECS expansion */
  size_t secs_flags;    /* MAGFIT_SECS_xxx flags */

  /* PCA parameters */
  size_t pca_modes;     /* number of PCA modes to use */

  size_t flags;         /* MAGFIT_FLG_xxx */
  double qdmax;         /* latitude range for fitting data: [-qdmax,qdmax] */
} magfit_parameters;

typedef struct
{
  const char *name;
  void * (*alloc)(const void * params);
  int (*reset) (void * state);
  size_t (*ncoeff) (void * state);
  int (*add_datum) (const double t, const double r, const double theta, const double phi, const double qdlat,
                    const double B[3], void * state);
  int (*fit) (double * rnorm, double * snorm, void * state);
  int (*eval_B) (const double r, const double theta, const double phi, double B[3], void * state);
  int (*eval_J) (const double r, const double theta, const double phi, double J[3], void * state);
  double (*eval_chi) (const double theta, const double phi, void * state);
  void (*free) (void * state);
} magfit_type;

typedef struct
{
  const magfit_type *type;
  magfit_parameters params;
  void *state;
} magfit_workspace;

const magfit_type * magfit_secs1d;
const magfit_type * magfit_secs2d;
const magfit_type * magfit_pca;
const magfit_type * magfit_gaussint;
const magfit_type * magfit_rc;

/*
 * Prototypes
 */

magfit_workspace *magfit_alloc(const magfit_type * T, const magfit_parameters * params);
void magfit_free(magfit_workspace *w);
magfit_parameters magfit_default_parameters(void);
int magfit_reset(magfit_workspace *w);
size_t magfit_add_track(track_data *tptr, const satdata_mag *data, magfit_workspace *w);
int magfit_fit(double * rnorm, double * snorm, magfit_workspace *w);
int magfit_apply_track(track_data *tptr, satdata_mag *data, magfit_workspace *w);
int magfit_eval_B(const double r, const double theta, const double phi, double B[3], magfit_workspace *w);
int magfit_eval_J(const double r, const double theta, const double phi, double J[3], magfit_workspace *w);
int magfit_print_track(const int header, FILE *fp, const track_data *tptr, const satdata_mag *data,
                       magfit_workspace *w);
int magfit_print_rms(const int header, FILE *fp, const double lon0, const track_data *tptr,
                     const satdata_mag *data, magfit_workspace *w);
int magfit_print_map(FILE *fp, const double r, magfit_workspace *w);

#endif /* INCLUDED_magfit_h */
