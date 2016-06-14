/*
 * main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <libconfig.h>

#include <satdata/satdata.h>

#include "cfg.h"
#include "common.h"
#include "mag.h"
#include "track.h"

static int
fill_parameters(mag_params *params)
{
  int s = 0;

  if (cfg_params.core_cof_file != NULL)
    params->core_file = (char *) cfg_params.core_cof_file;
  if (cfg_params.lith_cof_file != NULL)
    params->lith_file = (char *) cfg_params.lith_cof_file;
  if (cfg_params.kp_data_file != NULL)
    params->kp_file = (char *) cfg_params.kp_data_file;
  if (cfg_params.f107_data_file != NULL)
    params->f107_file = (char *) cfg_params.f107_data_file;
  if (cfg_params.dst_data_file != NULL)
    params->dst_file = (char *) cfg_params.dst_data_file;
  if (cfg_params.log_dir != NULL)
    params->log_dir = (char *) cfg_params.log_dir;
  if (cfg_params.lt_min >= 0.0)
    params->lt_min = cfg_params.lt_min;
  if (cfg_params.lt_max >= 0.0)
    params->lt_max = cfg_params.lt_max;
  if (cfg_params.max_lat_gap >= 0.0)
    params->dlat_max = cfg_params.max_lat_gap;
  if (cfg_params.r_earth >= 0.0)
    params->r_earth = cfg_params.r_earth;
  if (cfg_params.curr_altitude >= 0.0)
    params->curr_altitude = cfg_params.curr_altitude;
  if (cfg_params.kp_max >= 0.0)
    params->kp_max = cfg_params.kp_max;
  if (cfg_params.lon_min >= -1000.0)
    params->lon_min = cfg_params.lon_min;
  if (cfg_params.lon_max >= -1000.0)
    params->lon_max = cfg_params.lon_max;
  if (cfg_params.ncurr >= 0)
    params->ncurr = (size_t) cfg_params.ncurr;
  if (cfg_params.sq_nmax_int >= 0)
    params->sq_nmax_int = (size_t) cfg_params.sq_nmax_int;
  if (cfg_params.sq_mmax_int >= 0)
    params->sq_mmax_int = (size_t) cfg_params.sq_mmax_int;
  if (cfg_params.sq_nmax_ext >= 0)
    params->sq_nmax_ext = (size_t) cfg_params.sq_nmax_ext;
  if (cfg_params.sq_mmax_ext >= 0)
    params->sq_mmax_ext = (size_t) cfg_params.sq_mmax_ext;
  if (cfg_params.calc_field_models >= 0)
    params->calc_field_models = cfg_params.calc_field_models;
  if (cfg_params.main_nmax_int >= 0)
    params->main_nmax_int = cfg_params.main_nmax_int;

  return s;
}

static int
parse_config_file(const char *filename, mag_params *params)
{
  int s;
  config_t cfg;
  config_setting_t *setting;
  double fval;
  char *strptr;

  config_init(&cfg);

  s = config_read_file(&cfg, filename);
  if (!s)
    {
      fprintf(stderr, "parse_config_file: %s:%d - %s\n",
              config_error_file(&cfg),
              config_error_line(&cfg),
              config_error_text(&cfg));
      config_destroy(&cfg);
      return -1;
    }

  if (config_lookup_float(&cfg, "qdlat_max", &fval))
    params->qdlat_max = fval;

  config_destroy(&cfg);

  return 0;
}

static void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --champ_file | -c file              - CHAMP index file\n");
  fprintf(stderr, "\t --swarm_file | -s file              - Swarm index file\n");
  fprintf(stderr, "\t --output_file | -o file             - output file\n");
  fprintf(stderr, "\t --log_dir | -l dir                  - log directory\n");
  fprintf(stderr, "\t --lt_min | -a lt_min                - minimum local time (hours)\n");
  fprintf(stderr, "\t --lt_max | -b lt_max                - maximum local time (hours)\n");
  fprintf(stderr, "\t --lon_min | -d lon_min              - minimum longitude (degrees)\n");
  fprintf(stderr, "\t --lon_max | -e lon_max              - maximum longitude (degrees)\n");
  fprintf(stderr, "\t --season_min | -f season_min        - minimum season (doy)\n");
  fprintf(stderr, "\t --season_max | -g season_max        - maximum season (doy)\n");
  fprintf(stderr, "\t --season_min2 | -h season_min2      - minimum season 2 (doy)\n");
  fprintf(stderr, "\t --season_max2 | -j season_max2      - maximum season 2 (doy)\n");
  fprintf(stderr, "\t --sq_nmax_int | -t nmax             - nmax for Sq internal\n");
  fprintf(stderr, "\t --sq_mmax_int | -u mmax             - mmax for Sq internal\n");
  fprintf(stderr, "\t --sq_nmax_ext | -v nmax             - nmax for Sq external\n");
  fprintf(stderr, "\t --sq_mmax_ext | -w mmax             - mmax for Sq external\n");
  fprintf(stderr, "\t --curr_alt | -k curr_alt            - altitude of line current shell (km)\n");
  fprintf(stderr, "\t --qdmax | -q qdmax                  - maximum QD latitude for line currents (deg)\n");
  fprintf(stderr, "\t --profiles_only | -p                - compute magnetic/current profiles only (no EEF)\n");
  fprintf(stderr, "\t --vector | -z                       - use vector data instead of scalar\n");
}

int
main(int argc, char *argv[])
{
  satdata_mag *data = NULL;
  struct timeval tv0, tv1;
  mag_workspace *mag_workspace_p;
  track_workspace *track_workspace_p;
  cfg_workspace *config_workspace_p = NULL;
  mag_params params;

  params.kp_max = 20.0;
  params.year = -1; /* filled in below */
  params.log_dir = "log";
  params.output_file = NULL;
  params.curr_altitude = 110.0;
  params.ncurr = 161;
  params.qdlat_max = 20.0;
  params.lt_min = 6.0;
  params.lt_max = 18.0;
  params.lon_min = -200.0;
  params.lon_max = 200.0;
  params.season_min = 0.0;
  params.season_max = 367.0;
  params.season_min2 = -1.0;
  params.season_max2 = -1.0;
  params.profiles_only = 0;
  params.use_vector = 0;
  params.dlat_max = 2.0;
  params.r_earth = 6371.2;

  /* spherical harmonic degrees for Sq filter */
  params.sq_nmax_int = 12;
  params.sq_mmax_int = 0;
  params.sq_nmax_ext = 1;
  params.sq_mmax_ext = 1;

  params.kp_file = KP_IDX_FILE;
  params.f107_file = F107_IDX_FILE;
  params.dst_file = DST_IDX_FILE;

  params.calc_field_models = 0;
  params.core_file = NULL;
  params.lith_file = NULL;
  params.main_nmax_int = 15;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "champ_file", required_argument, NULL, 'c' },
          { "swarm_file", required_argument, NULL, 's' },
          { "output_file", required_argument, NULL, 'o' },
          { "log_dir", required_argument, NULL, 'l' },
          { "lt_min", required_argument, NULL, 'a' },
          { "lt_max", required_argument, NULL, 'b' },
          { "lon_min", required_argument, NULL, 'd' },
          { "lon_max", required_argument, NULL, 'e' },
          { "season_min", required_argument, NULL, 'f' },
          { "season_max", required_argument, NULL, 'g' },
          { "season_min2", required_argument, NULL, 'h' },
          { "season_max2", required_argument, NULL, 'j' },
          { "curr_alt", required_argument, NULL, 'k' },
          { "ncurr", required_argument, NULL, 'm' },
          { "qdmax", required_argument, NULL, 'q' },
          { "sq_nmax_int", required_argument, NULL, 't' },
          { "sq_mmax_int", required_argument, NULL, 'u' },
          { "sq_nmax_ext", required_argument, NULL, 'v' },
          { "sq_mmax_ext", required_argument, NULL, 'w' },
          { "vector", required_argument, NULL, 'z' },
          { "kp_max", required_argument, NULL, 'A' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:d:e:f:g:h:j:k:l:m:o:pq:s:t:u:v:w:zA:C:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            params.lt_min = atof(optarg);
            break;

          case 'b':
            params.lt_max = atof(optarg);
            break;

          case 'c':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_champ_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            if (!data)
              exit(1);
            fprintf(stderr, "done (%zu points read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            /* check for instrument flags */
            {
              size_t nflag;

              fprintf(stderr, "main: filtering for instrument flags...");
              nflag = satdata_champ_filter_instrument(1, SATDATA_FLG_ONESC, data);
              fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
                      nflag, data->n, (double)nflag / (double)data->n * 100.0);
            }

            break;

          case 'd':
            params.lon_min = atof(optarg);
            break;

          case 'e':
            params.lon_max = atof(optarg);
            break;

          case 'f':
            params.season_min = atof(optarg);
            break;

          case 'g':
            params.season_max = atof(optarg);
            break;

          case 'h':
            params.season_min2 = atof(optarg);
            break;

          case 'j':
            params.season_max2 = atof(optarg);
            break;

          case 'k':
            params.curr_altitude = atof(optarg);
            break;

          case 's':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_swarm_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            if (!data)
              exit(1);
            fprintf(stderr, "done (%zu points read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            {
              size_t nflag;

              fprintf(stderr, "main: filtering instrument flags...");
              nflag = satdata_swarm_filter_instrument(1, data);
              fprintf(stderr, "done (%zu data flagged)\n", nflag);
            }

            break;

          case 'o':
            params.output_file = optarg;
            break;

          case 'l':
            params.log_dir = optarg;
            break;

          case 'm':
            params.ncurr = (size_t) atoi(optarg);
            break;

          case 'p':
            params.profiles_only = 1;
            break;

          case 'q':
            params.qdlat_max = atof(optarg);
            break;

          case 't':
            params.sq_nmax_int = (size_t) atoi(optarg);
            break;

          case 'u':
            params.sq_mmax_int = (size_t) atoi(optarg);
            break;

          case 'v':
            params.sq_nmax_ext = (size_t) atoi(optarg);
            break;

          case 'w':
            params.sq_mmax_ext = (size_t) atoi(optarg);
            break;

          case 'z':
            params.use_vector = 1;
            break;

          case 'A':
            params.kp_max = atof(optarg);
            break;

          case 'C':
            fprintf(stderr, "main: reading config file %s...", optarg);
            config_workspace_p = cfg_alloc(optarg);
            fprintf(stderr, "done\n");
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  if (config_workspace_p)
    fill_parameters(&params);

  fprintf(stderr, "main: maximum allowed kp:      %.1f\n", params.kp_max);
  fprintf(stderr, "main: line current shell:      %.1f [km]\n", params.curr_altitude);
  fprintf(stderr, "main: number of line currents: %zu\n", params.ncurr);
  fprintf(stderr, "main: LT min:                  %.1f [h]\n", params.lt_min);
  fprintf(stderr, "main: LT max:                  %.1f [h]\n", params.lt_max);
  fprintf(stderr, "main: longitude min:           %.1f [deg]\n", params.lon_min);
  fprintf(stderr, "main: longitude max:           %.1f [deg]\n", params.lon_max);
  fprintf(stderr, "main: season min:              %.1f [deg]\n", params.season_min);
  fprintf(stderr, "main: season max:              %.1f [deg]\n", params.season_max);
  fprintf(stderr, "main: season min 2:            %.1f [deg]\n", params.season_min2);
  fprintf(stderr, "main: season max 2:            %.1f [deg]\n", params.season_max2);
  fprintf(stderr, "main: Sq internal nmax         %zu\n", params.sq_nmax_int);
  fprintf(stderr, "main: Sq internal mmax         %zu\n", params.sq_mmax_int);
  fprintf(stderr, "main: Sq external nmax         %zu\n", params.sq_nmax_ext);
  fprintf(stderr, "main: Sq external mmax         %zu\n", params.sq_mmax_ext);

  track_workspace_p = track_alloc();

  track_init(data, NULL, track_workspace_p);

  params.year = (int) satdata_epoch2year(data->t[0]);
  mag_workspace_p = mag_alloc(&params);

  mag_proc(&params, track_workspace_p, data, mag_workspace_p);

  mag_free(mag_workspace_p);
  satdata_mag_free(data);
  track_free(track_workspace_p);

  if (config_workspace_p)
    cfg_free(config_workspace_p);

  return 0;
} /* main() */
