/*
 * mfield_eval_main
 *
 * Example program for calling DMSP-MAG-1
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mfield_eval.h"

int
main(int argc, char *argv[])
{
  mfield_eval_workspace *w;
  const double m_pi = 3.14159265358979;
  double r, theta, phi, t;
  double E_st, I_st;
  double B[4];

  w = mfield_eval_read("dmsp_mag_1.txt");
  if (!w)
    exit(1);

  r = w->R + 300.0;
  theta = 45.0 * m_pi / 180.0;
  phi = 20.0 * m_pi / 180.0;
  t = 2013.2;

  mfield_eval(t, r, theta, phi, B, w);

  fprintf(stderr, "internal field:\n");
  fprintf(stderr, "X = %.2f [nT]\n", B[0]);
  fprintf(stderr, "Y = %.2f [nT]\n", B[1]);
  fprintf(stderr, "Z = %.2f [nT]\n", B[2]);
  fprintf(stderr, "F = %.2f [nT]\n", B[3]);

  E_st = -10.0;
  I_st = -3.5;
  mfield_eval_ext(t, r, theta, phi, E_st, I_st, B, w);

  fprintf(stderr, "external field:\n");
  fprintf(stderr, "X = %.2f [nT]\n", B[0]);
  fprintf(stderr, "Y = %.2f [nT]\n", B[1]);
  fprintf(stderr, "Z = %.2f [nT]\n", B[2]);
  fprintf(stderr, "F = %.2f [nT]\n", B[3]);

  mfield_eval_free(w);

  return 0;
} /* main() */
