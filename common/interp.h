/*
 * interp.h
 */

#ifndef INCLUDED_interp_h
#define INCLUDED_interp_h

#include "bsearch.h"

/*
 * Prototypes
 */

double interp1d(double a, double b, double fa, double fb, double x);
double interp_xy(const double x_array[], const double y_array[], const size_t n,
                 const double x);
double interp_cosine_1d(double a, double b, double fa, double fb, double x);
double interp2d(double x, double y, double gridpts[4][2],
                double funcvals[4]);
double interp3d(double x, double y, double z, double gridpts[8][3],
                double funcvals[8]);

#endif /* INCLUDED_interp_h */
