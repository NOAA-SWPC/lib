/*
 * efi.h
 */

#ifndef INCLUDED_efi_h
#define INCLUDED_efi_h

#include <time.h>
#include <gsl/gsl_math.h>

#define EFI_MAX_DATA     500000
#define EFI_MAX_BUFFER   2048

#define EFI_DATA_FILE    "/nfs/satmag_work/palken/corr/efi/efi.idx"

typedef struct
{
  time_t t;         /* timestamp */
  double longitude; /* longitude (degrees) */
  double glat;      /* geodetic latitude (degrees) */
  double galt;      /* geodetic altitude (km) */
  double E_phi;     /* eastward electric field (mV/m) */
} efi_datum;

typedef struct
{
  efi_datum *array;
  size_t n;    /* total data stored */
  size_t ntot; /* total data allocated */
} efi_data;

efi_data *efi_alloc(const size_t n);
void efi_free(efi_data *data);
efi_data *efi_read_idx(const char *idx_filename, efi_data *data);
int efi_sort(efi_data *data);
int efi_search(const time_t t, const double window, const efi_data *data,
               size_t *index);

#endif /* INCLUDED_efi_h */
