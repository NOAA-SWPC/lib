/*
 * magdata.h
 */

#ifndef INCLUDED_magdata_h
#define INCLUDED_magdata_h

#include <satdata/satdata.h>

#include "track_weight.h"

/* data flags */
#define MAGDATA_FLG_X                 (1 << 0)  /* X measurement available */
#define MAGDATA_FLG_Y                 (1 << 1)  /* Y measurement available */
#define MAGDATA_FLG_Z                 (1 << 2)  /* Z measurement available */
#define MAGDATA_FLG_F                 (1 << 3)  /* F measurement available */
#define MAGDATA_FLG_GRAD_NS           (1 << 4)  /* along-track difference (gradient) data available */
#define MAGDATA_FLG_DX                (1 << 5)  /* DX measurement available */
#define MAGDATA_FLG_DY                (1 << 6)  /* DY measurement available */
#define MAGDATA_FLG_DZ                (1 << 7)  /* DZ measurement available */
#define MAGDATA_FLG_DF                (1 << 8)  /* DF measurement available */
#define MAGDATA_FLG_TRACK_START       (1 << 9)  /* start of track in data structure */
#define MAGDATA_FLG_FIT_MF            (1 << 10) /* fit main field to this data point */
#define MAGDATA_FLG_FIT_EULER         (1 << 11) /* fit Euler angles to this data point */
#define MAGDATA_FLG_DISCARD           (1 << 12) /* discard this data point */

/* global flags */
#define MAGDATA_GLOBFLG_EULER         (1 << 0)  /* fit Euler angles to this dataset */

/* k_b * mu_0 in units of: nT^2 cm^3 / K */
#define MAGDATA_KB_MU0                (1.73497445090703e-05)

/* check if data point is discarded */
#define MAGDATA_Discarded(x)          ( (x) & MAGDATA_FLG_DISCARD )

/* vector measurement available */
#define MAGDATA_ExistVector(x)        (!MAGDATA_Discarded(x) && \
                                       ((x) & MAGDATA_FLG_X) && \
                                       ((x) & MAGDATA_FLG_Y) && \
                                       ((x) & MAGDATA_FLG_Z))

/* scalar measurement available */
#define MAGDATA_ExistScalar(x)        (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_F))

/* check if fitting Euler angles to this data point */
#define MAGDATA_FitEuler(x)           (MAGDATA_ExistVector(x) && \
                                       ((x) & MAGDATA_FLG_FIT_EULER))

typedef struct
{
  double t;              /* timestamp (CDF_EPOCH) */
  double r;              /* radius (km) */
  double theta;          /* geocentric colatitude (radians) */
  double phi;            /* longitude (radians) in [-pi,pi] */
  double qdlat;          /* QD latitude (degrees) */
  double B_nec[3];       /* magnetic field vector in NEC frame (nT) */
  double B_vfm[3];       /* magnetic field vector in VFM frame (nT) */
  double B_model[3];     /* main/crustal/external field model in NEC frame (nT) */
  double F;              /* magnetic field scalar measurement (nT) */
  double q[4];           /* quaternions */
  double ne;             /* electron density cm^{-3} */
  int satdir;            /* +1 north, -1 south */

  double r_ns;           /* radius for along-track point (km) */
  double theta_ns;       /* geocentric (or QD) colatitude for along-track point (radians) */
  double phi_ns;         /* longitude for along-track point (radians) */
  double qdlat_ns;       /* QD latitude for along-track point (degrees) */
  double B_nec_ns[3];    /* along-track magnetic field vector in NEC frame (nT) */
  double B_model_ns[3];  /* along-track main/crustal/external field model in NEC frame (nT) */
  double F_ns;           /* along-track magnetic field scalar measurement (nT) */

  size_t flags;          /* MAGDATA_FLG_xxx flags */
} magdata_datum;

typedef struct
{
  double *t;           /* timestamp (CDF_EPOCH) */
  double *ts;          /* scaled timestamp for MF modeling (dimensionless), not written to output file */
  double *r;           /* radius (km) */
  double *theta;       /* geocentric colatitude (radians) */
  double *phi;         /* longitude (radians) */
  double *qdlat;       /* QD latitude (degrees) */
  double *Bx_nec;      /* NEC X measurement (nT) */
  double *By_nec;      /* NEC Y measurement (nT) */
  double *Bz_nec;      /* NEC Z measurement (nT) */
  double *Bx_vfm;      /* VFM X measurement (nT) */
  double *By_vfm;      /* VFM Y measurement (nT) */
  double *Bz_vfm;      /* VFM Z measurement (nT) */
  double *Bx_model;    /* NEC X main/crustal/external field (nT) */
  double *By_model;    /* NEC Y main/crustal/external field (nT) */
  double *Bz_model;    /* NEC Z main/crustal/external field (nT) */
  double *F;           /* scalar measurement (nT) */
  double *q;           /* quaternions */
  double *weights;     /* spatial weights */
  int *satdir;         /* +1 north, -1 south */

  double *ne;          /* electron density cm^{-3} */

  double *r_ns;        /* radius for along-track point (km) */
  double *theta_ns;    /* geocentric colatitude for along-track point (radians) */
  double *phi_ns;      /* longitude for along-track point (radians) */
  double *qdlat_ns;    /* QD latitude for along-track point (degrees) */
  double *Bx_nec_ns;   /* NEC along-track X measurement (nT) */
  double *By_nec_ns;   /* NEC along-track Y measurement (nT) */
  double *Bz_nec_ns;   /* NEC along-track Z measurement (nT) */
  double *Bx_model_ns; /* NEC along-track X main/crustal/external field (nT) */
  double *By_model_ns; /* NEC along-track Y main/crustal/external field (nT) */
  double *Bz_model_ns; /* NEC along-track Z main/crustal/external field (nT) */
  double *F_ns;        /* along-track scalar measurement (nT) */

  size_t *flags;       /* MAGDATA_FLG_xxx flags */

  size_t n;            /* total number of data */
  size_t ntot;         /* total array allocation */
  double R;            /* Earth radius (km) */
  double rmin;         /* minimum radius stored (km) */
  double rmax;         /* maximum radius stored (km) */
  size_t nx;           /* number of X measurements */
  size_t ny;           /* number of Y measurements */
  size_t nz;           /* number of Z measurements */
  size_t nf;           /* number of F measurements */
  size_t ndx;          /* number of DX measurements */
  size_t ndy;          /* number of DY measurements */
  size_t ndz;          /* number of DZ measurements */
  size_t ndf;          /* number of DF measurements */
  size_t nvec;         /* number of vector measurements */
  size_t ngrad;        /* number of along-track difference measurements */
  size_t nres;         /* number of total residuals */

  size_t global_flags; /* MAGDATA_GLOBFLG_xxx flags applying to all data */

  track_weight_workspace *weight_workspace_p;
} magdata;

/*
 * Prototypes
 */

magdata *magdata_alloc(const size_t n, const double R);
magdata *magdata_realloc(const size_t n, const double R, magdata *data);
void magdata_free(magdata *data);
int magdata_datum_init(magdata_datum *datum);
int magdata_add(const magdata_datum *datum, magdata *data);
int magdata_init(magdata *data);
int magdata_unit_weights(magdata *data);
int magdata_calc(magdata *data);
int magdata_print(const char *filename, const magdata *data);
int magdata_map(const char *filename, const magdata *data);
int magdata_residual(const size_t idx, double B[4], const magdata *data);
int magdata_residual_ns(const size_t idx, double B[4], const magdata *data);
int magdata_residual_dB_ns(const size_t idx, double dB[4], const magdata *data);
int magdata_t(double *t0, double *t1, const magdata *data);
int magdata_write(const char *filename, magdata *data);
magdata *magdata_read(const char *filename, magdata *data);
size_t magdata_ndiscard(const magdata *data);
size_t magdata_neuler(const magdata *data);
int magdata_clear(magdata *data);
int magdata_flag_t(const double t0, const double t1, magdata *data);

#endif /* INCLUDED_magdata_h */
