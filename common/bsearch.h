/*
 * bsearch.h
 */

#ifndef INCLUDED_bsearch_h
#define INCLUDED_bsearch_h

#include <time.h>

/*
 * Prototypes
 */

size_t bsearch_double(const double x_array[], double x,
                      size_t index_lo, size_t index_hi);
size_t bsearch_desc_double(const double x_array[], double x,
                           size_t index_lo, size_t index_hi);
size_t bsearch_timet(const time_t x_array[], const time_t x,
                     size_t index_lo, size_t index_hi);

#endif /* INCLUDED_bsearch_h */
