/*
 * tiegcm3d_read.c
 *
 * This module contains routines to read tiegcm NetCDF files
 *
 * Required fields in TIEGCM files:
 *
 * Dimensions:
 *
 * time
 * glon
 * glat
 *
 * Variables:
 * glon
 * glat
 * year
 * doy
 * ut
 * Bn_grd
 * Be_grd
 * Bu_grd
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
#include <gsl/gsl_vector.h>

#include <common/common.h>

#include "tiegcm3d.h"

/*
tiegcm3d_read()
  Read tiegcm NetCDF data file

Inputs: filename - NetCDF file
        data     - (input/output) where to store newly read data
                   If NULL, this data structure is allocated to store
                   new data. If not NULL, new data is appended to
                   pre-existing data

Return: pointer to data structure containing data
*/

tiegcm3d_data *
tiegcm3d_read(const char *filename, tiegcm3d_data *data)
{
  int status, status2;
  size_t cur_idx, i;
  int ncid;
  int timeid, hid, glonid, glatid;
  int doyid, utid, hvid, glonvid, glatvid, Jtid, Jpid;
  size_t ntot, ncur;
  size_t nt, nr, nlon, nlat;

  status = nc_open(filename, NC_NOWRITE, &ncid);
  if (status)
    {
      fprintf(stderr, "tiegcm3d_read: unable to open %s: %s\n",
              filename, nc_strerror(status));
      return data;
    }

  status = 0;
  status += nc_inq_dimid(ncid, "time", &timeid);
  status += nc_inq_dimid(ncid, "hgt_sph", &hid);

  /* when generating netcdf from matlab I had to rename glat/glon to nlat/nlon */
  status2 = nc_inq_dimid(ncid, "glon_sph", &glonid);
  if (status2)
    status2 = nc_inq_dimid(ncid, "nlon", &glonid);

  status2 = nc_inq_dimid(ncid, "glat_sph", &glatid);
  if (status2)
    status2 = nc_inq_dimid(ncid, "nlat", &glatid);

  status += status2;

  if (status)
    {
      fprintf(stderr, "tiegcm3d_read: error reading dimension ids: %s\n",
              nc_strerror(status));
      return data;
    }

  /* find dimension lengths */
  status = 0;
  status += nc_inq_dimlen(ncid, timeid, &nt);
  status += nc_inq_dimlen(ncid, hid, &nr);
  status += nc_inq_dimlen(ncid, glonid, &nlon);
  status += nc_inq_dimlen(ncid, glatid, &nlat);

  if (status)
    {
      fprintf(stderr, "tiegcm3d_read: error reading dimension lengths: %s\n",
              nc_strerror(status));
      return data;
    }

  status = 0;
  status += nc_inq_varid(ncid, "doy", &doyid);
  status += nc_inq_varid(ncid, "ut", &utid);
  status += nc_inq_varid(ncid, "hgt_sph", &hvid);
  status += nc_inq_varid(ncid, "glon_sph", &glonvid);
  status += nc_inq_varid(ncid, "glat_sph", &glatvid);
  status += nc_inq_varid(ncid, "Keast", &Jpid);
  status += nc_inq_varid(ncid, "Ksouth", &Jtid);

  if (status)
    {
      fprintf(stderr, "tiegcm3d_read: error reading variable ids: %s\n",
              nc_strerror(status));
      return data;
    }

  /* allocate more data if needed */
  ntot = data ? data->ntot : 0;
  ncur = data ? data->nt : 0;
  if (ntot < ncur + nt)
    data = tiegcm3d_realloc(nt + ntot, nr, nlon, nlat, data);

  if (!data)
    return 0;

  cur_idx = data->nt;

  if (nt > data->nt_max)
    {
      fprintf(stderr, "tiegcm3d_read: error: nt_max too small [%zu]\n", data->nt_max);
      return 0;
    }

  status = 0;
  status += nc_get_var_double(ncid, doyid, &(data->doy[cur_idx]));
  status += nc_get_var_double(ncid, utid, &(data->ut[cur_idx]));
  status += nc_get_var_double(ncid, hvid, data->r);
  status += nc_get_var_double(ncid, glonvid, data->glon);
  status += nc_get_var_double(ncid, glatvid, data->glat);
  status += nc_get_var_double(ncid, Jtid, data->Jt);
  status += nc_get_var_double(ncid, Jpid, data->Jp);

  if (status)
    {
      fprintf(stderr, "tiegcm3d_read: error reading variables: %s\n",
              nc_strerror(status));
      return data;
    }

  data->nt += nt;

  /* convert heights to km radii */
  for (i = 0; i < nr; ++i)
    {
      double height = data->r[i];
      data->r[i] = 6371.2 + height * 1.0e-3;
    }

  /* compute timestamps */
  putenv("TZ=GMT");

  for (i = cur_idx; i < data->nt; ++i)
    {
      int iyear = 2009;
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
