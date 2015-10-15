#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void get_ext(double g1,double  g2,double  g3,double  phi,double  theta,double  rrel,
             double *x, double *y, double *z) 

                                /* ext deg 1 has no radial dependence */
{
  *x = g2 * cos(theta) * cos(phi) - g1 * sin(theta) + g3 * cos(theta) * sin(phi);
  *y = g2 * sin(phi) - g3 * cos(phi);
  *z = g1 * cos(theta) + g2 * sin(theta) * cos(phi) + g3 * sin(theta) * sin(phi);
}



void get_int(double g1,double  g2,double  g3,double  phi,double  theta,double  rrel, 
             double *x, double *y, double *z) 
{
  double rm3;
  rm3 = 1 / (rrel*rrel*rrel);
  *x = rm3 * (g2  * cos(theta) * cos(phi) - g1 * sin(theta) + g3 * cos(theta) * sin(phi));
  *y = rm3 * (g2  * sin(phi) - g3 * cos(phi));
  *z = -2 * rm3 * (g1 * cos(theta) + g2 * sin(theta) * cos(phi) + g3 * sin(theta) * sin(phi));
}










