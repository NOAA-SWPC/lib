/*
 * main.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <satdata/satdata.h>

#include "common.h"
#include "mag.h"
#include "track.h"

void
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
  fprintf(stderr, "\t --curr_alt | -k curr_alt            - altitude of line current shell (km)\n");
  fprintf(stderr, "\t --qdmax | -q qdmax                  - maximum QD latitude for line currents (deg)\n");
  fprintf(stderr, "\t --profiles_only | -p                - compute magnetic/current profiles only (no EEF)\n");
}

int
main(int argc, char *argv[])
{
  /*
   * rms thresholds; only need to check scalar field; set threshold high to process data
   * during strong storms: 17 March 2015 storm has scalar rms of up to 120 nT
   */
  const double thresh_swarm[] = { -1.0, -1.0, -1.0, 120.0 };
  const double thresh_champ[] = { -1.0, -1.0, -1.0, 120.0 };

  const double *thresh = NULL;
  satdata_mag *data = NULL;
  struct timeval tv0, tv1;
  mag_workspace *mag_workspace_p;
  track_workspace *track_workspace_p;
  mag_params params;

  params.year = -1; /* filled in below */
  params.log_dir = "log";
  params.output_file = NULL;
  params.curr_altitude = 110.0;
  params.ncurr = 81;
  params.qdlat_max = 20.0;
  params.lt_min = MAG_LT_MIN;
  params.lt_max = MAG_LT_MAX;
  params.lon_min = -200.0;
  params.lon_max = 200.0;
  params.season_min = 0.0;
  params.season_max = 367.0;
  params.season_min2 = -1.0;
  params.season_max2 = -1.0;
  params.profiles_only = 0;

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
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:d:e:h:j:k:l:m:o:ps:", long_options, &option_index);
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
            thresh = thresh_champ;
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
            thresh = thresh_swarm;
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

  fprintf(stderr, "main: maximum allowed kp:      %.1f\n", MAG_MAX_KP);
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

  track_workspace_p = track_alloc();

  track_init(data, NULL, track_workspace_p);

  /* flag tracks with high rms */
  {
    size_t nrms = track_flag_rms("rms.dat", thresh, data, track_workspace_p);
    fprintf(stderr, "main: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
            nrms, data->n, (double) nrms / (double) data->n * 100.0);
  }

  params.year = (int) satdata_epoch2year(data->t[0]);
  mag_workspace_p = mag_alloc(&params);

  mag_proc(&params, track_workspace_p, data, mag_workspace_p);

  mag_free(mag_workspace_p);
  satdata_mag_free(data);
  track_free(track_workspace_p);

  return 0;
} /* main() */
