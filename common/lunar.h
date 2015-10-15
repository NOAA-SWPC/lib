/*
 * lunar.h
 * Patrick Alken
 */

#ifndef INCLUDED_lunar_h
#define INCLUDED_lunar_h

#include <time.h>

typedef struct
{
  double phase; /* lunar phase (radians) new = 0, full = pi */
  double age;   /* moon's age from new (days) */
  double eclat; /* ecliptic latitude (radians) */
  double eclon; /* ecliptic longitude (radians) */
  double dist;  /* distance in units of earth radii */
} lunardata;

double lunartime(time_t t, double longitude);
double lunartime_r(time_t t, double longitude);
double lunarphase(time_t t);
int lunarcalc(int year, int month, int day, lunardata *data);

#endif /* INCLUDED_lunar_h */
