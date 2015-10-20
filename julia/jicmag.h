/*
 * jicmag.h
 */

#ifndef INCLUDED_jicmag_h
#define INCLUDED_jicmag_h

#include <time.h>

#define JICMAG_MAX_DATA     1000000
#define JICMAG_MAX_BUFFER   2048

#define JIC_LONGITUDE       (-76.87)

typedef struct
{
  time_t t; /* timestamp */
  double H; /* horizontal field (nT) */
} jicmag_datum;

typedef struct
{
  jicmag_datum *array;
  size_t n;    /* total data stored */
  size_t ntot; /* total data allocated */
} jicmag_data;

jicmag_data *jicmag_alloc(const size_t n);
void jicmag_free(jicmag_data *data);
jicmag_data *jicmag_read_idx(const char *idx_filename);
int jicmag_sort(jicmag_data *data);

#endif /* INCLUDED_jicmag_h */
