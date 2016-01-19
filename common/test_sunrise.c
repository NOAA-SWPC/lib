/*
 * test_zenith.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "solarpos.h"
#include "solpos00.h"
#include "tcalc.h"

int
main(int argc, char *argv[])
{
  long retval;
  struct posdata pd;
  int c;
  double lat, lon;
  time_t ts;
  double sr1, ss1, sr2, ss2;
  double fday;
  double epsr, epss;
  double epsr_max, epss_max;
  solarpos_workspace *spw;
  int hr1, mr1, secr1;
  int hs1, ms1, secs1;

  lat = -10.0;
  lon = -180.0;

  spw = solarpos_alloc();

  while ((c = getopt(argc, argv, "y:h:m:s:a:o:n:d:")) != (-1))
    {
      switch (c)
        {
          case 'a':
            lat = atoi(optarg);
            break;

          case 'o':
            lon = atoi(optarg);
            break;

          case '?':
          default:
            printf("usage: %s [-a lat] [-o lon]\n", argv[0]);
            exit(1);
        }
    }

  S_init(&pd);

  pd.longitude = (float) lon;
  pd.latitude = (float) lat;
  pd.timezone = 0.0;
  pd.year = 1970;
  pd.daynum = 1;
  pd.hour = 19;
  pd.minute = 0;
  pd.second = 0;

  pd.temp = 10.0;
  pd.press = 1013.0;
  pd.aspect = 180.0;

  retval = S_solpos(&pd);
  S_decode(retval, &pd);

  printf("sunrise = %f, sunset = %f\n",
         pd.sretr / 60.0,
         pd.ssetr / 60.0);

  epsr_max = -1.0;
  epss_max = -1.0;

  for (lat = -10.0; lat < 10.0; lat += 0.5) {
  for (lon = -180.0; lon < 180.0; lon += 1.0) {

  for (ts = 0; ts < 2000000000; ts += 86400)
    {
      solarpos_calc_sunrs(ts, lat*M_PI/180.0, lon*M_PI/180.0,
                          &sr1, &ss1, spw);

      sr1 /= 60.0;
      ss1 /= 60.0;

      fday = time2fday(ts);
      sr2 = tcalc_sunrise(fday, lat*M_PI/180.0, lon*M_PI/180.0);
      ss2 = tcalc_sunset(fday, lat*M_PI/180.0, lon*M_PI/180.0);

      sr2 /= 60.0;
      ss2 /= 60.0;

      epsr = fabs(sr1 - sr2);
      epss = fabs(ss1 - ss2);
      epsr *= 60.0;
      epss *= 60.0;

      if (epsr > epsr_max)
        epsr_max = epsr;
      if (epss > epss_max)
        epss_max = epss;

#if 0
      hr1 = (int) sr1;
      mr1 = (int) ((sr1 - hr1) * 60.0);
      secr1 = (int) ((sr1 - hr1 - mr1/60.0) * 3600.0);

      hs1 = (int) ss1;
      ms1 = (int) ((ss1 - hs1) * 60.0);
      secs1 = (int) ((ss1 - hs1 - ms1/60.0) * 3600.0);

      printf("%f %f %ld %02d:%02d:%02d %02d:%02d:%02d %02d:%02d %02d:%02d %f %f\n",
             lat,
             lon,
             ts,
             hr1,
             mr1,
             secr1,
             hs1,
             ms1,
             secs1,
             (int)sr2,
             (int) ((sr2 - (int)sr2) * 60.0),
             (int)ss2,
             (int) ((ss2 - (int)ss2) * 60.0),
             epsr,
             epss);
#endif
    }

  }
  }

  fprintf(stderr, "epsr_max = %f mins (%f hrs)\n",
          epsr_max, epsr_max / 60.0);
  fprintf(stderr, "epss_max = %f mins (%f hrs)\n",
          epss_max, epss_max / 60.0);

  solarpos_free(spw);

  return 0;
} /* main() */
