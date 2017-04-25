/*
 * interp.c
 *
 * Interpolation routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "bsearch.h"
#include "interp.h"

/*
interp1d()
  1D interpolation: find f(x) using endpoints f(a) and f(b)
*/

double
interp1d(double a, double b, double fa, double fb, double x)
{
  double fx;

  fx = (b - x) / (b - a) * fa + (x - a) / (b - a) * fb;

  return fx;
}

double
interp_xy(const double x_array[], const double y_array[], const size_t n,
          const double x)
{
  size_t idx = bsearch_double(x_array, x, 0, n - 1);
  double y = interp1d(x_array[idx], x_array[idx + 1], y_array[idx], y_array[idx + 1], x);

  return y;
}

double
interp_mu2(double a, double b, double x)
{
  return ((x - a) / (b - a));
}

/* interpolate f(x) between f(a) and f(b) using cosines */
double
interp_cosine_1d(double a, double b, double fa, double fb, double x)
{
  double mu = interp_mu2(a, b, x);
  double mu2 = (1.0 - cos(mu * M_PI)) / 2.0;

  return (fa * (1.0 - mu2) + fb * mu2);
} /* interp_cosine_1d() */

int equals( double a, double b, double tolerance )
{
    return ( a == b ) ||
      ( ( a <= ( b + tolerance ) ) &&
        ( a >= ( b - tolerance ) ) );
}

double cross2( double x0, double y0, double x1, double y1 )
{
    return x0*y1 - y0*x1;
}

int in_range( double val, double range_min, double range_max, double tol )
{
    return ((val+tol) >= range_min) && ((val-tol) <= range_max);
}

/* Returns number of solutions found.  If there is one valid solution, it will be put in s and t */
int inverseBilerp( double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double* sout, double* tout, double* s2out, double* t2out )
{
    int t_valid, t2_valid;

    double a  = cross2( x0-x, y0-y, x0-x2, y0-y2 );
    double b1 = cross2( x0-x, y0-y, x1-x3, y1-y3 );
    double b2 = cross2( x1-x, y1-y, x0-x2, y0-y2 );
    double c  = cross2( x1-x, y1-y, x1-x3, y1-y3 );
    double b  = 0.5 * (b1 + b2);

    double s, s2, t, t2;

    double am2bpc = a-2*b+c;
    /* this is how many valid s values we have */
    int num_valid_s = 0;

    if ( equals( am2bpc, 0, 1e-10 ) )
    {
        if ( equals( a-c, 0, 1e-10 ) )
        {
                /* Looks like the input is a line */
                /* You could set s=0.5 and solve for t if you wanted to */
                return 0;
        }
        s = a / (a-c);
        if ( in_range( s, 0, 1, 1e-10 ) )
                num_valid_s = 1;
    }
    else
    {
        double sqrtbsqmac = sqrt( b*b - a*c );
        s  = ((a-b) - sqrtbsqmac) / am2bpc;
        s2 = ((a-b) + sqrtbsqmac) / am2bpc;
        num_valid_s = 0;
        if ( in_range( s, 0, 1, 1e-10 ) )
        {
                num_valid_s++;
                if ( in_range( s2, 0, 1, 1e-10 ) )
                        num_valid_s++;
        }
        else
        {
                if ( in_range( s2, 0, 1, 1e-10 ) )
                {
                        num_valid_s++;
                        s = s2;
                }
        }
    }

    if ( num_valid_s == 0 )
        return 0;

    t_valid = 0;
    if ( num_valid_s >= 1 )
    {
        double tdenom_x = (1-s)*(x0-x2) + s*(x1-x3);
        double tdenom_y = (1-s)*(y0-y2) + s*(y1-y3);
        t_valid = 1;
        if ( equals( tdenom_x, 0, 1e-10 ) && equals( tdenom_y, 0, 1e-10 ) )
        {
                t_valid = 0;
        }
        else
        {
                /* Choose the more robust denominator */
                if ( fabs( tdenom_x ) > fabs( tdenom_y ) )
                {
                        t = ( (1-s)*(x0-x) + s*(x1-x) ) / ( tdenom_x );
                }
                else
                {
                        t = ( (1-s)*(y0-y) + s*(y1-y) ) / ( tdenom_y );
                }
                if ( !in_range( t, 0, 1, 1e-10 ) )
                        t_valid = 0;
        }
    }

    /* Same thing for s2 and t2 */
    t2_valid = 0;
    if ( num_valid_s == 2 )
    {
        double tdenom_x = (1-s2)*(x0-x2) + s2*(x1-x3);
        double tdenom_y = (1-s2)*(y0-y2) + s2*(y1-y3);
        t2_valid = 1;
        if ( equals( tdenom_x, 0, 1e-10 ) && equals( tdenom_y, 0, 1e-10 ) )
        {
                t2_valid = 0;
        }
        else
        {
                /* Choose the more robust denominator */
                if ( fabs( tdenom_x ) > fabs( tdenom_y ) )
                {
                        t2 = ( (1-s2)*(x0-x) + s2*(x1-x) ) / ( tdenom_x );
                }
                else
                {
                        t2 = ( (1-s2)*(y0-y) + s2*(y1-y) ) / ( tdenom_y );
                }
                if ( !in_range( t2, 0, 1, 1e-10 ) )
                        t2_valid = 0;
        }
    }

    /* Final cleanup */
    if ( t2_valid && !t_valid )
    {
        s = s2;
        t = t2;
        t_valid = t2_valid;
        t2_valid = 0;
    }

    /* Output */
    if ( t_valid )
    {
        *sout = s;
        *tout = t;
    }

    if ( t2_valid )
    {
        *s2out = s2;
        *t2out = t2;
    }

    return t_valid + t2_valid;
}

void bilerp( double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double s, double t, double* x, double* y )
{
    *x = t*(s*x3+(1-s)*x2) + (1-t)*(s*x1+(1-s)*x0);
    *y = t*(s*y3+(1-s)*y2) + (1-t)*(s*y1+(1-s)*y0);
}

/*
interp2d()
  2D interpolation: find f(x, y) where (x, y) lies in the square
defined by gridpoints (a,b,c,d), where:

gridpts[0] = a = (x1,y1)         c----------d
gridpts[1] = b = (x2,y2)         | P(x,y)   |
gridpts[2] = c = (x3,y3)         |          |
gridpts[3] = d = (x4,y4)         a----------b

The grid is not necessarily equally spaced (y1 != y2 for ex), but it
is assumed that (x,y) lies inside the (rough) square shown above.

with function values given:

funcvals[0] = f(a)
funcvals[1] = f(b)
funcvals[2] = f(c)
funcvals[3] = f(d)
*/

double
interp2d(double x, double y, double gridpts[4][2], double funcvals[4])
{
  double f;
  double s, t;
  double s2, t2;
  int num_st;

  /* check inputs */
  if (x < gridpts[0][0] || x < gridpts[2][0] ||
      x > gridpts[1][0] || x > gridpts[3][0])
    {
      fprintf(stderr, "interp2d: x = %e does not lie in the square of a = %e, b = %e, c = %e, d = %e\n",
              x,
              gridpts[0][0],
              gridpts[1][0],
              gridpts[2][0],
              gridpts[3][0]);
      return 0.0;
    }

  num_st = inverseBilerp(gridpts[0][0], gridpts[0][1],
                         gridpts[1][0], gridpts[1][1],
                         gridpts[2][0], gridpts[2][1],
                         gridpts[3][0], gridpts[3][1],
                         x, y,
                         &s, &t,
                         &s2, &t2);

  if (num_st != 1)
    {
      fprintf(stderr, "interp2d: failed to find s, t\n");
      return 0.0;
    }

  f = t * (s * funcvals[3] + (1.0 - s) * funcvals[2]) +
      (1.0 - t) * (s * funcvals[1] + (1.0 - s) * funcvals[0]);

  return f;
} /* interp2d() */

/*
interp3d()
  3D interpolation: find f(x, y, z) where (x, y, z) lies in the cube
defined by gridpoints (a,b,c,d,e,f,g,h), where:

gridpts[0] = a = (x1,y1,z1)
gridpts[1] = b = (x2,y2,z1)
gridpts[2] = c = (x3,y3,z1)
gridpts[3] = d = (x4,y4,z1)
gridpts[4] = e = (x5,y5,z2)
gridpts[5] = f = (x6,y6,z2)
gridpts[6] = g = (x7,y7,z2)
gridpts[7] = h = (x8,y8,z2)

(a,b,c,d) define the bottom square (lower z = z1) using the same conventions
as interp2d(). (e,f,g,h) define the top square (z = z2)

The grid is not necessarily equally spaced in the (x, y) coordinates,
but is assumed to be uniform in z

with function values given:

funcvals[0] = f(a)
funcvals[1] = f(b)
funcvals[2] = f(c)
funcvals[3] = f(d)
funcvals[4] = f(e)
funcvals[5] = f(f)
funcvals[6] = f(g)
funcvals[7] = f(h)
*/

double
interp3d(double x, double y, double z, double gridpts[8][3],
         double funcvals[8])
{
  double f;
  double f1, f2;
  double s, t;
  double s2, t2;
  int num_st;
  double z1, z2;

  num_st = inverseBilerp(gridpts[0][0], gridpts[0][1],
                         gridpts[1][0], gridpts[1][1],
                         gridpts[2][0], gridpts[2][1],
                         gridpts[3][0], gridpts[3][1],
                         x, y,
                         &s, &t,
                         &s2, &t2);

  if (num_st != 1)
    {
      fprintf(stderr, "interp2d: failed to find s, t\n");
      return 0.0;
    }

  f1 = t * (s * funcvals[3] + (1.0 - s) * funcvals[2]) +
       (1.0 - t) * (s * funcvals[1] + (1.0 - s) * funcvals[0]);

  num_st = inverseBilerp(gridpts[4][0], gridpts[4][1],
                         gridpts[5][0], gridpts[5][1],
                         gridpts[6][0], gridpts[6][1],
                         gridpts[7][0], gridpts[7][1],
                         x, y,
                         &s, &t,
                         &s2, &t2);

  if (num_st != 1)
    {
      fprintf(stderr, "interp2d: failed to find s, t\n");
      return 0.0;
    }

  f2 = t * (s * funcvals[7] + (1.0 - s) * funcvals[6]) +
       (1.0 - t) * (s * funcvals[5] + (1.0 - s) * funcvals[4]);

  /* assume uniform spacing in z coordinate */
  z1 = gridpts[0][2];
  z2 = gridpts[4][2];
  f = interp1d(z1, z2, f1, f2, z);

  return f;
} /* interp3d() */
