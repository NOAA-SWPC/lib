/*
 * tiegcm_read.c
 *
 * This module contains routines to read tiegcm NetCDF files
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <netcdf.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "tiegcm.h"

/*
tiegcm_read()
  Read tiegcm NetCDF data file

Inputs: filename - NetCDF file
        data     - (input/output) where to store newly read data
                   If NULL, this data structure is allocated to store
                   new data. If not NULL, new data is appended to
                   pre-existing data

Return: pointer to data structure containing data
*/

tiegcm_data *
tiegcm_read(const char *filename, tiegcm_data *data)
{
  int status;
  size_t cur_idx, i;
  int ncid;
  int yearid, doyid, utid;
  size_t nrec, ntot, ncur;

  status = nc_open(filename, NC_NOWRITE, &ncid);
  if (status)
    {
      fprintf(stderr, "tiegcm_read: unable to open %s: %s\n",
              filename, nc_strerror(status));
      return data;
    }

  status = 0;
  status += nc_inq_varid(ncid, "year", &yearid);
  status += nc_inq_varid(ncid, "doy", &doyid);
  status += nc_inq_varid(ncid, "ut", &utid);

  if (status)
    {
      fprintf(stderr, "tiegcm_read: error reading variable ids: %s\n",
              nc_strerror(status));
      return data;
    }

  /* since records start at 0 */
  nrec = 142;

  /* allocate more data if needed */
  ntot = data ? data->ntot : 0;
  ncur = data ? data->n : 0;
  if (ntot < ncur + nrec)
    data = tiegcm_realloc(nrec + ntot, data);

  if (!data)
    return 0;

  cur_idx = data->n;

  status = 0;
  status += nc_get_var_double(ncid, yearid, &(data->year[cur_idx]));
  status += nc_get_var_double(ncid, doyid, &(data->doy[cur_idx]));
  status += nc_get_var_double(ncid, utid, &(data->ut[cur_idx]));

  if (status)
    {
      fprintf(stderr, "tiegcm_read: error reading variables: %s\n",
              nc_strerror(status));
      return data;
    }

  data->n += nrec;

  /* compute timestamps */
  putenv("TZ=GMT");

  for (i = cur_idx; i < data->n; ++i)
    {
      int iyear = (int) data->year[i];
      int idoy = (int) data->doy[i];
      int month, day, hour, min, sec;
      struct tm tminfo;

      doy2md(iyear, idoy, &month, &day);

      hour = (int) data->ut[i];
      min = (int) ((data->ut[i] - hour) * 60.0);
      sec = 0;

      tminfo.tm_sec = sec;
      tminfo.tm_min = min;
      tminfo.tm_hour = hour;
      tminfo.tm_mday = day;
      tminfo.tm_mon = month - 1;
      tminfo.tm_year = iyear - 1900;
      tminfo.tm_isdst = 0;

      data->t[i] = mktime(&tminfo);
    }

  nc_close(ncid);

  return data;
}
