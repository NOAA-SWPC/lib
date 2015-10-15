/*
 * test_values.c
 *
 * Generate test values for a model coefficient file
 *
 * Usage: ./test_values -e epoch_year
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "geo.h"

#include "msynth.h"

/* output lots of points for all years for general testing */
#define TEST_ALL       0

/* output less points for GUI testing */
#define TEST_GUI       1

#define NMAX           740
#define START_YEAR     1980

/*
print_test_values()
  Print test values for model from given time/position

Inputs: header - print header
        altd   - geodetic altitude (km)
        latd   - geodetic latitude (deg)
        lond   - geodetic longitude (deg)
        t      - time (decimal years)
        w      - workspace
*/

static void
print_test_values(const int header, const double altd, const double latd,
                  const double lond, const double t,
                  msynth_workspace *w)
{
  double r, lat; /* geocentric position */
  double theta, phi;
  double B[4], dBdt[4];
  double D, I; /* declination/inclination (degrees) */
  double X, Y, Z, H, F;
  double dD, dI;
  double dX, dY, dZ, dH, dF;

  if (header)
    {
      size_t i = 1;

      printf("# Field %zu: Date\n", i++);
      printf("# Field %zu: WGS84 Altitude (km)\n", i++);
      printf("# Field %zu: geodetic latitude (deg)\n", i++);
      printf("# Field %zu: geodetic longitude (deg)\n", i++);
      printf("# Field %zu: declination (deg)\n", i++);
      printf("# Field %zu: inclination (deg)\n", i++);
      printf("# Field %zu: H (nT)\n", i++);
      printf("# Field %zu: X (nT)\n", i++);
      printf("# Field %zu: Y (nT)\n", i++);
      printf("# Field %zu: Z (nT)\n", i++);
      printf("# Field %zu: F (nT)\n", i++);
      printf("# Field %zu: dD/dt (deg/year)\n", i++);
      printf("# Field %zu: dI/dt (deg/year)\n", i++);
      printf("# Field %zu: dH/dt (nT/year)\n", i++);
      printf("# Field %zu: dX/dt (nT/year)\n", i++);
      printf("# Field %zu: dY/dt (nT/year)\n", i++);
      printf("# Field %zu: dZ/dt (nT/year)\n", i++);
      printf("# Field %zu: dF/dt (nT/year)\n", i++);
      fflush(stdout);
      return;
    }

  geodetic2geo(latd * M_PI / 180.0, altd, &lat, &r);

  theta = M_PI / 2.0 - lat;
  phi = lond * M_PI / 180.0;

  msynth_eval(t, r, theta, phi, B, w);
  msynth_eval_dBdt(t, r, theta, phi, dBdt, w);

  /* rotate from geocentric back to ellipsoidal frame (eq 16-17 of WMM report) */
  {
    double sd = sin(lat - latd*M_PI/180.0);
    double cd = cos(lat - latd*M_PI/180.0);
    double Bx = B[0];
    double dBx = dBdt[0];

    B[0] = Bx * cd - B[2] * sd;
    B[2] = Bx * sd + B[2] * cd;

    dBdt[0] = dBx * cd - dBdt[2] * sd;
    dBdt[2] = dBx * sd + dBdt[2] * cd;
  }

  X = B[0];
  Y = B[1];
  Z = B[2];
  H = gsl_hypot(X, Y);
  F = gsl_hypot(H, Z);
  D = atan2(Y, X) * 180.0 / M_PI; /* convert to deg */
  I = atan2(Z, H) * 180.0 / M_PI; /* convert to deg */

  dX = dBdt[0];
  dY = dBdt[1];
  dZ = dBdt[2];
  dH = (X*dX + Y*dY) / H;
  dF = (X*dX + Y*dY + Z*dZ) / F;
  dI = (H*dZ - Z*dH) / (F*F) * 3437.74677; /* convert to minutes */
  dD = (X*dY - Y*dX) / (H*H) * 3437.74677; /* convert to minutes */

#if TEST_ALL
  /* 2 decimal place for HDGM testing */
  printf("%6.1f %6g %5g %5g %6.1f %6.1f %8.1f %8.1f %8.1f %8.1f %8.1f %4.1f %4.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
         t,
         altd,
         latd,
         lond,
         D,
         I,
         H,
         X,
         Y,
         Z,
         F,
         dD / 60.0,
         dI / 60.0,
         dH,
         dX,
         dY,
         dZ,
         dF);
#elif TEST_GUI
  /* 0 decimal places for HDGM gui testing */
  printf("%6.1f %6g %5g %5g %6.2f %6.2f %8.0f %8.0f %8.0f %8.0f %8.0f %4.2f %4.2f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
         t,
         altd,
         latd,
         lond,
         D,
         I,
         H,
         X,
         Y,
         Z,
         F,
         dD / 60.0,
         dI / 60.0,
         dH,
         dX,
         dY,
         dZ,
         dF);
#elif 0
  /* 12 decimal places for debugging */
  printf("%6.1f %6g %5g %5g %6.2f %6.2f %8.12f %8.12f %8.12f %8.12f %8.12f %4.1f %4.1f %5.12f %5.12f %5.12f %5.12f %5.12f\n",
         t,
         altd,
         latd,
         lond,
         D,
         I,
         H,
         X,
         Y,
         Z,
         F,
         dD,
         dI,
         dH,
         dX,
         dY,
         dZ,
         dF);
#endif
  fflush(stdout);
} /* print_test_values() */

int
main(int argc, char *argv[])
{
  int c;
  msynth_workspace *msynth_p = NULL;
  char *outfile = NULL;
  double t0 = -1.0, t1;

  while ((c = getopt(argc, argv, "c:w:o:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            fprintf(stderr, "main: reading coefficient file %s...", optarg);
            msynth_p = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'w':
            fprintf(stderr, "main: reading WMM coefficient file %s...", optarg);
            msynth_p = msynth_wmm_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!msynth_p)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-w wmm_file] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (t0 < 0.0)
    t0 = msynth_p->epochs[0];

  t1 = t0 + 5.0;

  {
    double t;

    /* print header */
    print_test_values(1, 0.0, 0.0, 0.0, t0, msynth_p);

#if TEST_ALL
    for (t = t0; t < t1; t += 0.1)
      {
        print_test_values(0, 0.0, 80.0, 120.0, t, msynth_p);
        print_test_values(0, 0.0, 0.0, 120.0, t, msynth_p);
        print_test_values(0, 0.0, -80.0, 120.0, t, msynth_p);
        print_test_values(0, 0.0, 80.0, 0.0, t, msynth_p);
        print_test_values(0, 0.0, 0.0, 0.0, t, msynth_p);
        print_test_values(0, 0.0, -80.0, 0.0, t, msynth_p);
        print_test_values(0, 0.0, 80.0, -120.0, t, msynth_p);
        print_test_values(0, 0.0, 0.0, -120.0, t, msynth_p);
        print_test_values(0, 0.0, -80.0, -120.0, t, msynth_p);

        print_test_values(0, 305.0, 80.0, 120.0, t, msynth_p);
        print_test_values(0, 305.0, 0.0, 120.0, t, msynth_p);
        print_test_values(0, 305.0, -80.0, 120.0, t, msynth_p);
        print_test_values(0, 305.0, 80.0, 0.0, t, msynth_p);
        print_test_values(0, 305.0, 0.0, 0.0, t, msynth_p);
        print_test_values(0, 305.0, -80.0, 0.0, t, msynth_p);
        print_test_values(0, 305.0, 80.0, -120.0, t, msynth_p);
        print_test_values(0, 305.0, 0.0, -120.0, t, msynth_p);
        print_test_values(0, 305.0, -80.0, -120.0, t, msynth_p);
      }
#elif TEST_GUI
    {
      gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

      gsl_rng_set(r, time(NULL));

      for (t = t0 + 0.1; t < t1; t += 0.6)
        {
          double alt = floor(gsl_rng_uniform(r) * 300); /* random alt in [0,300] km */

          print_test_values(0, alt, 80.0, 120.0, t, msynth_p);
          print_test_values(0, alt, 0.0, -10.0, t, msynth_p);
          print_test_values(0, alt, -80.0, 120.0, t, msynth_p);
        }

      gsl_rng_free(r);
    }
#endif
  }

  msynth_free(msynth_p);

  return 0;
}
