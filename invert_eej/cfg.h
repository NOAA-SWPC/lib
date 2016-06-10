/*
 * cfg.h
 *
 * Configuration file parameters
 */

#ifndef INCLUDED_cfg_h
#define INCLUDED_cfg_h

#include <libconfig.h>

#define CFG_STRING      (1 << 0)
#define CFG_INT         (1 << 1)
#define CFG_DOUBLE      (1 << 2)
#define CFG_OPTIONAL    (1 << 3) /* optional parameter */

typedef struct
{
  const char *prev_input_filename; /* EEF MAGxLR_1B previous day input file (CDF) */
  const char *input_filename;  /* EEF MAGxLR_1B input file (CDF) */
  const char *output_filename; /* EEF output file (CDF) */
  const char *log_dir;         /* log directory prefix */
  const char *test_output_dir; /* test output directory */
  const char *kp_data_file;    /* KP data file */
  const char *f107_data_file;  /* F10.7 data file */
  const char *dst_data_file;   /* DST data file */
  const char *igrf_cof_file;   /* IGRF coefficient file */
  const char *lith_cof_file;   /* Lithospheric field coef file */
  const char *core_cof_file;   /* Core field coefficient file */
  const char *pmf_cof_file;    /* location of PMF coefficient file */
  const char *pmf_sm_file;     /* location of POMME SM file */
  const char *pmf_gsm_file;    /* location of POMME GSM file */
  int pomme_residual_deg;      /* POMME main field degree for residuals */
  int pde_int_deg;             /* POMME main field degree for PDE solver */
  double lt_min;               /* minimum allowed local time */
  double lt_max;               /* maximum allowed local time */
  double max_lat_gap;          /* maximum allowed latitude gap (degrees) */
  int dst_num_days;            /* number of days of Dst data for Est/Ist */
  double curr_altitude;        /* altitude of current shell (km) */
  int ncurr;                   /* number of line currents */
  double prof_lat_min;         /* minimum latitude for profile inversion */
  double prof_lat_max;         /* maximum latitude for profile inversion */
  double prof_lat_step;        /* latitude step for profile inversion */
  int nr;                      /* radial grid points */
  int ntheta;                  /* theta grid points */
  double r_earth;              /* radius of earth in km */
  double alt_min;              /* minimum altitude for PDE solver (km) */
  double alt_max;              /* maximum altitude for PDE solver (km) */
  double theta_min;            /* minimum colatitude for PDE solver (deg) */
  double theta_max;            /* maximum colatitude for PDE solver (deg) */
  double kp_max;               /* maximum allowed kp */
  double lon_min;              /* minimum longitude (degrees) */
  double lon_max;              /* maximum longitude (degrees) */
  int sq_nmax_int;             /* internal spherical harmonic degree for Sq filter */
  int sq_mmax_int;             /* internal spherical harmonic order for Sq filter */
  int sq_nmax_ext;             /* external spherical harmonic degree for Sq filter */
  int sq_mmax_ext;             /* external spherical harmonic order for Sq filter */
} cfg_parameters;

typedef struct
{
  const char *name;
  void *location;
  size_t flags;
} cfg_struct;

typedef struct
{
  cfg_parameters params;
  config_t config;
} cfg_workspace;

/*
 * Prototypes
 */
cfg_workspace *cfg_alloc(const char *config_filename);
void cfg_free(cfg_workspace *w);

extern cfg_parameters cfg_params;

#endif /* INCLUDED_cfg_h */
