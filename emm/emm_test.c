/*
 * emm_test.c
 *
 * Generate test values for EMM
 *
 * Usage: ./emm_test -e epoch_year [-c msynth_coef_file] [-n nmax]
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
#define TEST_ALL       1

/* output less points for GUI testing */
#define TEST_GUI       0

#define NMAX           740
#define START_YEAR     1900

#define indx(n,m)     ( (n) * ((n) + 1) / 2 + (m) )

/*
print_test_values()
  Print test values for EMM model from given time/position

Inputs: header - print header
        altd   - geodetic altitude (km)
        latd   - geodetic latitude (deg)
        lond   - geodetic longitude (deg)
        t      - time (decimal years)
*/

static void
print_test_values(const int header, const double altd, const double latd,
                  const double lond, const double t,
                  msynth_workspace *msynth_core, msynth_workspace *msynth_crust)
{
  double r, lat; /* geocentric position */
  double theta, phi;
  double B[4], dBdt[4], B_crust[4];
  double D, I; /* declination/inclination (degrees) */
  double X, Y, Z, H, F;
  double dD, dI;
  double dX, dY, dZ, dH, dF;
  size_t i;

  if (header)
    {
      i = 1;
      printf("# Field %zu: Date\n", i++);
      printf("# Field %zu: altitude (km)\n", i++);
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

  msynth_eval(t, r, theta, phi, B, msynth_core);
  msynth_eval_dBdt(t, r, theta, phi, dBdt, msynth_core);
  msynth_eval(t, r, theta, phi, B_crust, msynth_crust);

  /* add main and crustal fields */
  for (i = 0; i < 3; ++i)
    B[i] += B_crust[i];

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
  /* 1 decimal place for EMM testing */
  printf("%5.1f %6.1f %6.1f %6.1f %7.2f %6.2f %9.1f %9.1f %9.1f %9.1f %9.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n",
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
  /* 0 decimal places for EMM gui testing */
#if 0
  printf("%6.1f %6g %5g %5g %6.2f %6.2f %8.0f %8.0f %8.0f %8.0f %8.0f %4.2f %4.2f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
#else
  /* tab-delimited for Nir */
  printf("%.1f\t%g\t%g\t%g\t%.2f\t%.2f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.2f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
#endif
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
  msynth_workspace *msynth_core = NULL;
  msynth_workspace *msynth_crust;
  size_t nmax = NMAX;
  int year = -1;
  double epoch;

  while ((c = getopt(argc, argv, "e:n:")) != (-1))
    {
      switch (c)
        {
          case 'e':
            year = atoi(optarg);
            break;

          case 'n':
            nmax = (size_t) atoi(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-e epoch_year> [-n nmax]\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (year < 0)
    {
      fprintf(stderr, "Usage: %s <-e epoch_year> [-n nmax] msynth1.txt msynth2.txt ...\n", argv[0]);
      exit(1);
    }

  epoch = (double) year;

  if (optind < argc)
    {
      size_t nmodels = argc - optind;

      fprintf(stderr, "main: allocating workspace for %zu models...",
              nmodels);
      msynth_core = msynth_alloc(15, nmodels, NULL);
      fprintf(stderr, "done\n");

      for (c = optind; c < argc; ++c)
        {
          fprintf(stderr, "main: reading %s...", argv[c]);
          msynth_core = msynth_read2(argv[c], msynth_core);
          fprintf(stderr, "done\n");
        }

      /* recompute SV as differences of main field snapshots */
      fprintf(stderr, "main: computing SV coefficients...");
      msynth_calc_sv(msynth_core);
      fprintf(stderr, "done\n");
    }
  else
    {
      fprintf(stderr, "main: error: no main field models specified\n");
      exit(1);
    }

  fprintf(stderr, "main: reading crustal field from %s...", MSYNTH_EMM_FILE);
  msynth_crust = msynth_emm_read(MSYNTH_EMM_FILE);
  fprintf(stderr, "done\n");

  msynth_set(1, nmax, msynth_core);

  {
#if 1
    const double test_lat[3] = { 78.3, 0.0, -85.9 };
    const double test_lon[3] = { 123.7, 0.0, -116.5 };
#else
    const double test_lat[3] = { 80.0, 0.0, -80.0 };
    const double test_lon[3] = { 120.0, 0.0, -120.0 };
#endif
    double t;
    double alt1 = 0.0;
    double alt2 = 100.0;

    /* print header */
    print_test_values(1, 0.0, 0.0, 0.0, epoch, msynth_core, msynth_crust);

#if TEST_ALL
    for (t = epoch; t < epoch + 0.9; t += 0.1)
      {
        print_test_values(0, alt1, test_lat[0], test_lon[0], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[1], test_lon[0], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[2], test_lon[0], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[0], test_lon[1], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[1], test_lon[1], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[2], test_lon[1], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[0], test_lon[2], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[1], test_lon[2], t, msynth_core, msynth_crust);
        print_test_values(0, alt1, test_lat[2], test_lon[2], t, msynth_core, msynth_crust);

        print_test_values(0, alt2, test_lat[0], test_lon[0], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[1], test_lon[0], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[2], test_lon[0], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[0], test_lon[1], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[1], test_lon[1], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[2], test_lon[1], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[0], test_lon[2], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[1], test_lon[2], t, msynth_core, msynth_crust);
        print_test_values(0, alt2, test_lat[2], test_lon[2], t, msynth_core, msynth_crust);
      }
#elif TEST_GUI
    {
      gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

      /*gsl_rng_set(r, time(NULL));*/

      for (t = epoch + 0.1; t < epoch + 0.9; t += 0.6)
        {
          double alt = floor(gsl_rng_uniform(r) * 501); /* random altitude in [0,500] km */

          print_test_values(0, alt, test_lat[0], test_lon[0], t, msynth_core, msynth_crust);
        }

      gsl_rng_free(r);
    }
#endif
  }

  msynth_free(msynth_core);
  msynth_free(msynth_crust);

  return 0;
}
