/*
 * magdata_list.h
 */

#ifndef INCLUDED_magdata_list_h
#define INCLUDED_magdata_list_h

#include <gsl/gsl_math.h>

#include "magdata.h"

#define MAGDATA_LIST_IDX_X                   0
#define MAGDATA_LIST_IDX_Y                   1
#define MAGDATA_LIST_IDX_Z                   2
#define MAGDATA_LIST_IDX_F                   3
#define MAGDATA_LIST_IDX_DX_NS               4
#define MAGDATA_LIST_IDX_DY_NS               5
#define MAGDATA_LIST_IDX_DZ_NS               6
#define MAGDATA_LIST_IDX_DF_NS               7
#define MAGDATA_LIST_IDX_DX_EW               8
#define MAGDATA_LIST_IDX_DY_EW               9
#define MAGDATA_LIST_IDX_DZ_EW               10
#define MAGDATA_LIST_IDX_DF_EW               11

#define MAGDATA_LIST_IDX_TOTAL               12 /* total data */
#define MAGDATA_LIST_IDX_VECTOR              13 /* total vector data */
#define MAGDATA_LIST_IDX_VGRAD               14 /* total vector gradient data */
#define MAGDATA_LIST_IDX_END                 15

typedef struct
{
  size_t n;       /* number of data sources (satellites, observatories) */
  size_t n_tot;   /* number of data sources allocated */
  magdata **mdata;
} magdata_list;

/*
 * Prototypes
 */

magdata_list *magdata_list_alloc(const size_t nsources);
void magdata_list_free(magdata_list *w);
size_t magdata_list_add(const char *filename, magdata_list *w);
size_t magdata_list_filter_time(const double tmin, const double tmax,
                                magdata_list *w);
size_t magdata_list_filter_euler(magdata_list *w);
int magdata_list_map(const char *dir_prefix, const magdata_list *w);
int magdata_list_print(const char *dir_prefix, const magdata_list *w);
magdata *magdata_list_ptr(const size_t idx, const magdata_list *w);
int magdata_list_rminmax(const magdata_list * list, double * rmin, double * rmax);
int magdata_list_count(const magdata_list * list, size_t count[]);

#endif /* INCLUDED_magdata_list_h */
