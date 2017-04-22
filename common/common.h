/* common.h
 * 
 * Copyright (C) 2006, 2007 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_common_h
#define INCLUDED_common_h

#include <time.h>
#include <sys/time.h>

#define R_EARTH_KM        (6371.2)
#define R_EARTH_M         (R_EARTH_KM * 1.0e3)

#define CIDX2(i,Ni,j,Nj)      ((j) + (Nj)*(i))
#define CIDX3(i,Ni,j,Nj,k,Nk) ((k) + (Nk)*((j) + (Nj)*(i)))

#define FIDX3(i,Ni,j,Nj,k,Nk) ((i) + (Ni)*((j) + (Nj)*(k)))

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

/*
 * Prototypes
 */

time_t fday2timet(double fday);
double time2fday(time_t t);
time_t date2timet(int sec, int min, int hour, int day, int month, int year);
double time2lunar(time_t t);
double get_season(time_t t);
double get_year(time_t t);
int get_doy(time_t t);
double get_localtime(time_t t, double longitude);
double get_ut(time_t t);
double get_zenith(double fday, double latitude, double longitude);
double get_fday(int year, double longitude, double local_time,
                int season);
int doy2md(int year, int doy, int *month, int *day);
int is_leap_year(int year);
double lon_fix(double phi);
double wrap360(double x);
double wrap180(double x);
double wrap2pi(double x);
double wrappi(double x);
double linint(double x, double a, double b, double fa, double fb);
double **array2d_alloc(size_t size_x, size_t size_y);
void array2d_free(size_t size_x, double **a);
int sphcross(double A[3], double B[3], double C[3]);
int sph_basis(const double theta, const double phi,
              double rhat[3], double that[3], double phat[3]);
void sph2ecef_vec(const double theta, const double phi, const double V_sph[3], double V_ecef[3]);
void ecef2sph_vec(const double theta, const double phi, const double V_ecef[3], double V_sph[3]);
double vec_dot(const double a[3], const double b[3]);
double vec_norm(const double v[3]);
int vec_unit(const double v[3], double unit[3]);
double time_diff(struct timeval a, struct timeval b);
int progress_bar(FILE *fp, const double progress, const size_t bar_width);

extern int putenv(char *string);

#endif /* INCLUDED_common_h */
