/*
 * eph_data.h
 */

#ifndef INCLUDED_eph_data_h
#define INCLUDED_eph_data_h

#define EPH_MAX_DATA         35000000

#define EPH_DATA_FLG_ECI     (1 << 0) /* ephemeris is ECI */
#define EPH_DATA_FLG_ECEF    (1 << 1) /* ephemeris is ECEF */

/* ephemeris data */
typedef struct
{
  double t[EPH_MAX_DATA];         /* timestamp (CDF_EPOCH) */
  double X[EPH_MAX_DATA];         /* X in km */
  double Y[EPH_MAX_DATA];         /* Y in km */
  double Z[EPH_MAX_DATA];         /* Z in km */
  double VX[EPH_MAX_DATA];        /* V_x in km/s */
  double VY[EPH_MAX_DATA];        /* V_y in km/s */
  double VZ[EPH_MAX_DATA];        /* V_z in km/s */
  double latitude[EPH_MAX_DATA];  /* geocentric latitude in deg */
  double longitude[EPH_MAX_DATA]; /* geocentric longitude in deg */
  size_t n;                       /* number of points stored */
  size_t flags;                   /* ECI or ECEF */
} eph_data;

/*
 * Prototypes
 */

void eph_data_free(eph_data *data);
eph_data *eph_data_read_bowman(const char *filename);
eph_data *eph_data_read_tena(const char *filename);

#endif /* INCLUDED_eph_data_h */
