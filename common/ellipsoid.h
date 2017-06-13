/*
 * ellipsoid.h
 */

#ifndef INCLUDED_ellipsoid_h
#define INCLUDED_ellipsoid_h

/* semi-major and semi-minor axes from WGS84 */
#define WGS84_A        (6378.137)
#define WGS84_B        (6356.7523142)

/*
 * mu distance parameter specifying Earth surface -
 * mu = atanh(WGS84_B / WGS84_A)
 */
#define WGS84_MU       (3.19471282344095)

/*
 * Prototypes
 */

int ellipsoid_basis(const double r[3], double e_mu[3], double e_nu[3],
                    double e_phi[3]);
int ellipsoid_basis_mu(const double r[3], const double mu, double e_mu[3],
                       double e_nu[3], double e_phi[3]);
int ellipsoid_nec2ell(const double r_ECEF[3], const double B_nec[3], double B_ell[3]);

#endif /* INCLUDED_ellipsoid_h */
