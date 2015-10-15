/*
 * julian.h
 */

#ifndef INCLUDED_julian_h
#define INCLUDED_julian_h

#include <time.h>

/*
 * Prototypes
 */

double date2julian(const int year, const int month, const int day,
                   const int hour, const int min, const double sec);
double timet2julian(const time_t t);
double julian2GMST(double jd);
double julian2GAST(double jd);

#endif /* INCLUDED_julian_h */
