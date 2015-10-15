/*
 * coord.h
 * Patrick Alken
 */

#ifndef INCLUDED_coord_h
#define INCLUDED_coord_h

#include "Gdef.h"

/*
 * Prototypes
 */

void trans(int action, double fday, double lon_rad, double lat_rad,
           double *long_rad_s, double *lat_rad_s);
void sphere2cart(double phi_rad, double lat_rad, double r, double X[3]);
void sphere2cart_vec(double phi_rad,double lat_rad,double vlat,double vphi,
                     double vminusr, double *vx,double *vy,double *vz);
void cart2sphere_vec(double phi_rad,double lat_rad,double vx,double vy,
                     double vz, double *vlat,double *vphi,double *vminusr);
void my_delaz(double lat1, double lon1, double lat2, double lon2,
              double *delta, double *az);
void trans_vec(int action,double fday,double lon_rad,double lat_rad,
               double vlat,double vphi,double vminusr,
               double *lon_rad_s,double *lat_rad_s,double *vlat_s,
               double *vphi_s,double *vminusr_s);

#endif /* INCLUDED_coord_h */
