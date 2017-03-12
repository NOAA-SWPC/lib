/*
 * mfield_data.h
 */

#ifndef INCLUDED_mfield_data_h
#define INCLUDED_mfield_data_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include "magdata.h"

/* selectively fit different components if needed for testing */
#define MFIELD_FIT_X          1
#define MFIELD_FIT_Y          1
#define MFIELD_FIT_Z          1
#define MFIELD_FIT_F          0
#define MFIELD_FIT_DX_NS      1
#define MFIELD_FIT_DY_NS      1
#define MFIELD_FIT_DZ_NS      1
#define MFIELD_FIT_DX_EW      1
#define MFIELD_FIT_DY_EW      1
#define MFIELD_FIT_DZ_EW      1

typedef struct
{
  size_t nsources;   /* number of data sources (satellites) */
  double epoch;      /* model epoch (decimal years) */
  magdata **mdata;

  double *t0;        /* array of size nsources for first time of each satellite */
  double *t1;        /* array of size nsources for last time of each satellite */

  double t_mu;       /* mean of timestamps (years) */
  double t_sigma;    /* stddev of timestamps (years) */
  double t0_data;    /* timestamp of first data point (CDF_EPOCH) */
  double t1_data;    /* timestamp of last data point (CDF_EPOCH) */

  gsl_rstat_workspace *rstat_workspace_p;
} mfield_data_workspace;

/*
 * Prototypes
 */

mfield_data_workspace *mfield_data_alloc(const size_t nsources,
                                         const double epoch);
void mfield_data_free(mfield_data_workspace *w);
int mfield_data_copy(const size_t sat_idx, satdata_mag *data,
                     const size_t flags, mfield_data_workspace *w);
size_t mfield_data_filter_time(const double tmin, const double tmax,
                               mfield_data_workspace *w);
size_t mfield_data_filter_euler(mfield_data_workspace *w);
size_t mfield_data_filter_comp(mfield_data_workspace *w);
int mfield_data_init(mfield_data_workspace *w);
double mfield_data_epoch(mfield_data_workspace *w);
int mfield_data_map(const char *dir_prefix, const mfield_data_workspace *w);
int mfield_data_print(const char *dir_prefix, const mfield_data_workspace *w);
magdata *mfield_data_ptr(const size_t idx, const mfield_data_workspace *w);
int mfield_data_t(double *t0, double *t1, const magdata *data);

#endif /* INCLUDED_mfield_data_h */
