/*
 * julia.h
 */

#ifndef INCLUDED_julia_h
#define INCLUDED_julia_h

#include <time.h>
#include <gsl/gsl_math.h>

#define JULIA_MAX_DATA         10000000
#define JULIA_MAX_BUFFER       2048

#define JULIA_LON              (-76.87)
#define JULIA_GEODETIC_LAT     (-11.95)
#define JULIA_GEOCENTRIC_LAT   (-11.872280)

#define JULIA_LON_RAD          (JULIA_LON * M_PI / 180.0)

#define JULIA_AVG_IDX_FILE "/nfs/satmag/palken/JULIA/julia_avg.idx"
#define JULIA_IDX_FILE     "/nfs/satmag/palken/JULIA/julia.idx"

typedef struct
{
  time_t t;       /* timestamp */
  double v_zonal; /* zonal drift (m/s) */
  double v_vert;  /* vertical drift (m/s) */
} julia_datum;

typedef struct
{
  julia_datum *array;
  size_t n;    /* total data stored */
  size_t ntot; /* total data allocated */
} julia_data;

julia_data *julia_alloc(const size_t n);
void julia_free(julia_data *data);
int julia_read_avg_data(const char *filename, julia_data *data);
julia_data *julia_read_idx(const char *idx_filename, julia_data *data);
julia_data *julia_read_avg_idx(const char *idx_filename, julia_data *data);
int julia_sort(julia_data *data);
int julia_search(const time_t t, const double window,
                 const julia_data *data, size_t *index);

#endif /* INCLUDED_julia_h */
