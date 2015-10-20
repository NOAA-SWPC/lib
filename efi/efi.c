/*
 * efi.c
 * Patrick Alken
 *
 * This module reads Swarm EFI ascii data files
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <satdata/satdata.h>

#include "efi.h"

static int efi_compare_sort(const void *a, const void *b);
static int efi_compare_search(const void *a, const void *b);

/* XXX Global */
double EFI_TIME_WINDOW = 0.0;

/*
efi_read_data()
  Read data from a EFI data file

Inputs: filename - data file
        w        - kp workspace
*/

int
efi_read_data(const char *filename, efi_data *data)
{
  int s = 0;
  FILE *fp;
  char buf[EFI_MAX_BUFFER];
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "efi_read_data: fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 1;
    }

  /* use GMT time */
  putenv("TZ=GMT");

  n = data->n;

  while (fgets(buf, sizeof(buf), fp) != NULL)
    {
      int c;
      double t;
      time_t t0 = -2208988800; /* Jan 1 1900 */
      double glat, phi, galt, Ephi;

      /* EFI timestamps are seconds since Jan 1 1900 */
      c = sscanf(buf, "%lf %lf %lf %lf %lf",
                 &t,
                 &glat,
                 &phi,
                 &galt,
                 &Ephi);
      if (c < 5)
        continue;

      data->array[n].t = (time_t) t + t0;
      data->array[n].glat = glat;
      data->array[n].longitude = phi;
      data->array[n].galt = galt;
      data->array[n].E_phi = Ephi;

      if (++n >= data->ntot)
        {
          fprintf(stderr, "efi_read_data: error: ntot too small\n");
          exit(1);
        }
    }

  fclose(fp);

  data->n = n;

  return s;
} /* efi_read_data() */

/*
efi_read_idx()
  Read all available EFI data

Inputs: idx_filename - index file
        data         - (output) append data here
                       or set to NULL for new
                       allocation
*/

efi_data *
efi_read_idx(const char *idx_filename, efi_data *data)
{
  int s;
  FILE *fp;
  char buffer[EFI_MAX_BUFFER];

  fp = fopen(idx_filename, "r");
  if (!fp)
    {
      fprintf(stderr, "efi_read_idx: unable to open index file %s: %s\n",
              idx_filename, strerror(errno));
      return 0;
    }

  if (data == NULL)
    data = efi_alloc(EFI_MAX_DATA);

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      buffer[strlen(buffer) - 1] = '\0';
      s = efi_read_data(buffer, data);
      if (s)
        break;
    }

  fclose(fp);

  return data;
} /* efi_read_idx() */

efi_data *
efi_alloc(const size_t n)
{
  efi_data *data;

  data = calloc(1, sizeof(efi_data));
  if (!data)
    return 0;

  data->array = malloc(n * sizeof(efi_datum));

  data->ntot = n;
  data->n = 0;

  return data;
} /* efi_alloc() */

void
efi_free(efi_data *data)
{
  if (data->array)
    free(data->array);

  free(data);
} /* efi_free() */

int
efi_sort(efi_data *data)
{
  qsort(data->array, data->n, sizeof(efi_datum), efi_compare_sort);
  return 0;
} /* efi_sort() */

/*
efi_search()
  Search EFI data for a measurement within a specified window of a
given time

Inputs: t      - timestamp to search for
        window - allowed time window
        data   - EFI data
        index  - (output) index into array of data if found

Return: 0 if found, -1 if not found
*/

int
efi_search(const time_t t, const double window, const efi_data *data,
           size_t *index)
{
  int s = 0;
  efi_datum *ptr;

  /* XXX Global variable */
  EFI_TIME_WINDOW = window;

  ptr = bsearch(&t, data->array, data->n, sizeof(efi_datum),
                efi_compare_search);
  if (ptr)
    {
      size_t idx = ptr - data->array;
      assert(ptr->t == data->array[idx].t);
      *index = idx;
    }
  else
    s = -1; /* not found */

  return s;
}

static int
efi_compare_sort(const void *a, const void *b)
{
  efi_datum *da = (efi_datum *) a;
  efi_datum *db = (efi_datum *) b;
  const time_t t1 = da->t;
  const time_t t2 = db->t;

  return (t1 - t2);
} /* efi_compare_sort() */

static int
efi_compare_search(const void *a, const void *b)
{
  const time_t t1 = *(const time_t *) a;
  efi_datum *db = (efi_datum *) b;
  const time_t t2 = db->t;
  double diff = (t2 - t1) / 60.0;

  if (fabs(diff) < EFI_TIME_WINDOW)
    return 0; /* found */
  else if (t1 < t2)
    return -1;
  else
    return 1;
} /* efi_compare_search() */
