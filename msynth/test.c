/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "common.h"
#include "msynth.h"

int
test_igrf(void)
{
  int s = 0;
  msynth_workspace *w = msynth_igrf_read(MSYNTH_IGRF_FILE);
  const double theta = M_PI / 2.0 - 20.0 * M_PI / 180.0;
  const double phi = -10.0 * M_PI / 180.0;
  const double r = w->R + 100.0;
  const double tol = 1.0e-12;
  double B[4];

  /* test values taken from Matlab IGRF implementation */

  msynth_eval(1900.100200152207, r, theta, phi, B, w);
  gsl_test_rel(B[0], 2.768833125295744e+04, tol, "1900 Bx");
  gsl_test_rel(B[1], -7.851025953873755e+03, tol, "1900 By");
  gsl_test_rel(B[2], 1.977217624607306e+04, tol, "1900 Bz");

  msynth_eval(1950.689241248097, r, theta, phi, B, w);
  gsl_test_rel(B[0], 2.921582815253614e+04, tol, "1950 Bx");
  gsl_test_rel(B[1], -6.207002781334470e+03, tol, "1950 By");
  gsl_test_rel(B[2], 1.549704479525150e+04, tol, "1950 Bz");

  msynth_eval(1993.270063165906, r, theta, phi, B, w);
  gsl_test_rel(B[0], 3.072569520786891e+04, tol, "1993 Bx");
  gsl_test_rel(B[1], -3.612575058057932e+03, tol, "1993 By");
  gsl_test_rel(B[2], 1.251608986114605e+04, tol, "1993 Bz");

  msynth_eval(2004.336037492410, r, theta, phi, B, w);
  gsl_test_rel(B[0], 3.099014890382297e+04, tol, "2004 Bx");
  gsl_test_rel(B[1], -2.867576485753952e+03, tol, "2004 By");
  gsl_test_rel(B[2], 1.207786099143877e+04, tol, "2004 Bz");

  msynth_eval(2008.046192319368, r, theta, phi, B, w);
  gsl_test_rel(B[0], 31119.0088326705481, tol, "2008 Bx");
  gsl_test_rel(B[1], -2649.4736155154269, tol, "2008 By");
  gsl_test_rel(B[2], 11924.6187095175737, tol, "2008 Bz");

  msynth_eval(2014.428093797565, r, theta, phi, B, w);
  gsl_test_rel(B[0], 31332.8470307078642, tol, "2014 Bx");
  gsl_test_rel(B[1], -2255.57980121877199, tol, "2014 By");
  gsl_test_rel(B[2], 11645.0908870569947, tol, "2014 Bz");

  msynth_free(w);

  return s;
}

int
test_wmm(void)
{
  int s = 0;
  msynth_workspace *w = msynth_wmm_read(MSYNTH_WMM_FILE);

  /*
   * Enter: geodetic lat = 30.0, lon = -65.2, alt = 0 in WMM software
   * to get these geocentric values
   */
  const double theta = M_PI / 2.0 - 29.833634633909668 * M_PI / 180.0;
  const double phi = -65.2 * M_PI / 180.0;
  const double r = 6372.7793888732622;

  const double tol = 1.0e-12;
  double B[4];

  msynth_eval(2010.0, r, theta, phi, B, w);
  gsl_test_rel(B[0], 24180.278568789774, tol, "2010 Bx");
  gsl_test_rel(B[1], -6306.7436574298436, tol, "2010 By");
  gsl_test_rel(B[2], 37197.167838020032, tol, "2010 Bz");

  msynth_eval(2011.35, r, theta, phi, B, w);
  gsl_test_rel(B[0], 24214.650050182827, tol, "2011 Bx");
  gsl_test_rel(B[1], -6300.4659913459882, tol, "2011 By");
  gsl_test_rel(B[2], 36975.024407508579, tol, "2011 Bz");

  msynth_eval(2013.882, r, theta, phi, B, w);
  gsl_test_rel(B[0], 24279.115673062242, tol, "2013 Bx");
  gsl_test_rel(B[1], -6288.6918798465049, tol, "2013 By");
  gsl_test_rel(B[2], 36558.382062282581, tol, "2013 Bz");

  msynth_free(w);

  return s;
}

int
test_chaos(void)
{
  int s = 0;
  msynth_workspace *w = msynth_chaos_read(MSYNTH_CHAOS_FILE);

  msynth_free(w);

  return s;
}

int
main(int argc, char *argv[])
{
  fprintf(stderr, "testing IGRF...");
  test_igrf();
  fprintf(stderr, "done\n");

  fprintf(stderr, "testing WMM...");
  test_wmm();
  fprintf(stderr, "done\n");

  fprintf(stderr, "testing CHAOS...");
  test_chaos();
  fprintf(stderr, "done\n");

  exit (gsl_test_summary());

  return 0;
} /* main() */
