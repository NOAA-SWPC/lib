/*
 * magdata_list.h
 */

#ifndef INCLUDED_magdata_list_h
#define INCLUDED_magdata_list_h

#include <gsl/gsl_math.h>

#include "magdata.h"

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
int magdata_list_print(const char *dir_prefix, const gsl_vector *wts_spatial, const magdata_list *w);
magdata *magdata_list_ptr(const size_t idx, const magdata_list *w);

#endif /* INCLUDED_magdata_list_h */
