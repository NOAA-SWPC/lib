/*
 * mfield_covar.c
 *
 * Rewrite binary coefficient in ascii format
 *
 * usage: mfield_covar -c binary_coef_file -t decimal_year
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "oct.h"
#include "mfield.h"
#include "msynth.h"

/* compute b^T A b */
static double
compute_sigmasq(const gsl_vector *b, const gsl_matrix *A, gsl_vector *work)
{
  double result;

  /* compute work = A b */
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, b, 0.0, work);

  /* compute result = b . work */
  gsl_blas_ddot(b, work, &result);

  return result;
}

/*
compute_covar()
  Compute global rms error in X,Y,Z using covariance matrix
of main field coefficients

Inputs: t        - time at which to compute rms error (decimal years)
        C        - covariance matrix of gnm coefficients
        rms      - (output)
                   rms[0] = sigmasq X rms over globe [nT^2]
                   rms[1] = sigmasq Y rms over globe [nT^2]
                   rms[2] = sigmasq Z rms over globe [nT^2]
                   rms[3] = sigmasq H rms over globe [nT^2]
                   rms[4] = sigmasq F rms over globe [nT^2]
                   rms[5] = sigmasq D rms over globe [rad^2]
                   rms[6] = sigmasq I rms over globe [rad^2]
                   rms[7] = sigmasq GV rms over globe [rad^2]
        mfield_p - mfield workspace
        msynth_p - msynth workspace

Notes:
1) rms are computed over the globe at Earth's surface
*/

static int
compute_covar(const double t, const gsl_matrix *C,
              double rms[3],
              mfield_workspace *mfield_p,
              msynth_workspace *msynth_p)
{
  int s = 0;
  const size_t nmax = 12;
  const size_t nnm = (nmax + 1) * (nmax + 1) - 1;
  const size_t p = nnm;
  const double t0 = mfield_p->epoch;
  const double dt = t - t0;
  gsl_matrix_const_view SC0 = gsl_matrix_const_submatrix(C, 0, 0, p, p);
  gsl_matrix_const_view SR0 = gsl_matrix_const_submatrix(C, mfield_p->sv_offset, mfield_p->sv_offset, p, p);
  gsl_matrix_const_view SC0R0 = gsl_matrix_const_submatrix(C, mfield_p->sv_offset, 0, p, p);
  gsl_matrix *SC = gsl_matrix_alloc(p, p);
  gsl_vector *work = gsl_vector_alloc(p);
  mfield_green_workspace *green_p = mfield_p->green_workspace_p;
  gsl_vector_view dX = gsl_vector_view_array(green_p->dX, p);
  gsl_vector_view dY = gsl_vector_view_array(green_p->dY, p);
  gsl_vector_view dZ = gsl_vector_view_array(green_p->dZ, p);
  double r, theta, phi;
  size_t n = 0;
  size_t ngv = 0;
  size_t i, j;

  for (i = 0; i < p; ++i)
    {
      for (j = 0; j < p; ++j)
        {
          double sc0 = gsl_matrix_get(&SC0.matrix, i, j);
          double sr0 = gsl_matrix_get(&SR0.matrix, i, j);
          double sc0r0 = gsl_matrix_get(&SC0R0.matrix, i, j);

          /* Sigma_C from eq (5) of Nikos paper */
          gsl_matrix_set(SC, i, j, sc0 + dt * dt * sr0 + 2.0 * dt * sc0r0);
        }
    }

  r = mfield_p->R;
  theta = M_PI / 2.0;
  phi = 0.0;

  for (i = 0; i < 8; ++i)
    rms[i] = 0.0;

  for (theta = 0.01; theta < M_PI; theta += 5.0 * M_PI / 180.0)
    {
      double lat = 90.0 - theta * 180.0 / M_PI;

      for (phi = 0.0; phi < 2.0 * M_PI; phi += 5.0 * M_PI / 180.0)
        {
          double SSX, SSY, SSZ, SSH, SSF, SSD, SSI;
          double B[4];

          /* compute main field for this point */
          msynth_eval(t, r, theta, phi, B, msynth_p);

          /* compute Green's functions for this point */
          mfield_green_calc(r, theta, phi, green_p);

          /* compute b^T SC b for X,Y,Z */
          SSX = compute_sigmasq(&dX.vector, SC, work);
          SSY = compute_sigmasq(&dY.vector, SC, work);
          SSZ = compute_sigmasq(&dZ.vector, SC, work);

          /* compute rms H, F, D, I */
          {
            double H = gsl_hypot(B[0], B[1]);
            double term1 = B[0] * (1.0 + (B[1]/B[0])*(B[1]/B[0]));
            double term2 = H * (1.0 + (B[2]/H)*(B[2]/H));
            double Dterm = 1.0 / (term1 * term1);
            double Iterm = 1.0 / (term2 * term2);

            SSH = (B[0] / H) * (B[0] / H) * SSX +
                  (B[1] / H) * (B[1] / H) * SSY;
            SSF = (B[0] / B[3]) * (B[0] / B[3]) * SSX +
                  (B[1] / B[3]) * (B[1] / B[3]) * SSY +
                  (B[2] / B[3]) * (B[2] / B[3]) * SSZ;
            SSD = (B[1] / B[0]) * (B[1] / B[0]) * Dterm * SSX +
                  Dterm * SSY;
            SSI = (B[2] / H) * (B[2] / H) * Iterm * SSH +
                  Iterm * SSZ;
          }

          rms[0] += SSX * SSX;
          rms[1] += SSY * SSY;
          rms[2] += SSZ * SSZ;
          rms[3] += SSH * SSH;
          rms[4] += SSF * SSF;
          rms[5] += SSD * SSD;
          rms[6] += SSI * SSI;
          ++n;

          if (lat > 55.0 || lat < -55.0)
            {
              rms[7] += SSD * SSD;
              ++ngv;
            }
        }
    }

  for (i = 0; i < 7; ++i)
    rms[i] = sqrt(rms[i] / (double)n);

  rms[7] = sqrt(rms[7] / (double)ngv);

  gsl_matrix_free(SC);
  gsl_vector_free(work);

  return s;
}

int
main(int argc, char *argv[])
{
  int c;
  char *binary_coef_file = NULL;
  char *ascii_coef_file = NULL;
  char *matrix_file = NULL;
  mfield_workspace *mfield_workspace_p;
  msynth_workspace *msynth_workspace_p;
  double t = -1.0;

  while ((c = getopt(argc, argv, "c:a:t:m:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            binary_coef_file = optarg;
            break;

          case 'a':
            ascii_coef_file = optarg;
            break;

          case 't':
            t = atof(optarg);
            break;

          case 'm':
            matrix_file = optarg;
            break;

          default:
            break;
        }
    }

  if (!binary_coef_file || !ascii_coef_file)
    {
      fprintf(stderr, "usage: %s <-c binary_coef_file> <-a ascii_coef_file> [-t decimal_year]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: loading coefficients from %s...", ascii_coef_file);
  msynth_workspace_p = msynth_read(ascii_coef_file);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: loading covariance matrix from %s...", binary_coef_file);
  mfield_workspace_p = mfield_read(binary_coef_file);
  fprintf(stderr, "done\n");

  if (matrix_file)
    {
      gsl_matrix *covar = mfield_workspace_p->covar;
      const size_t n = covar->size1;
      size_t i, j;
      FILE *fp;

      fprintf(stderr, "main: writing covariance matrix to %s...", matrix_file);

      fp = fopen(matrix_file, "w");

      /*print_octave(mfield_workspace_p->covar, matrix_file);*/

      /* format requested by Nikos */
      for (j = 0; j < n; ++j)
        {
          for (i = 0; i <= j; ++i)
            {
              double cij = gsl_matrix_get(covar, i, j);

              if (i != 0 && i % 4 == 0)
                fprintf(fp, "\n");

              fprintf(fp, "%25.15E ", cij);
            }
          fprintf(fp, "\n");
        }

      fclose(fp);

      fprintf(stderr, "done\n");
    }

  if (t < 0.0)
    t = mfield_workspace_p->epoch;

  fprintf(stderr, "main: time = %f\n", t);

  {
    double rms[7];

    fprintf(stderr, "main: computing rms errors of field components...");
    compute_covar(t, mfield_workspace_p->covar, rms,
                  mfield_workspace_p, msynth_workspace_p);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: rms sigma X = %f [nT]\n", sqrt(rms[0]));
    fprintf(stderr, "main: rms sigma Y = %f [nT]\n", sqrt(rms[1]));
    fprintf(stderr, "main: rms sigma Z = %f [nT]\n", sqrt(rms[2]));
    fprintf(stderr, "main: rms sigma H = %f [nT]\n", sqrt(rms[3]));
    fprintf(stderr, "main: rms sigma F = %f [nT]\n", sqrt(rms[4]));
    fprintf(stderr, "main: rms sigma D = %f [deg]\n", sqrt(rms[5]) * 180.0 / M_PI);
    fprintf(stderr, "main: rms sigma I = %f [deg]\n", sqrt(rms[6]) * 180.0 / M_PI);
    fprintf(stderr, "main: rms sigma GV = %f [deg]\n", sqrt(rms[7]) * 180.0 / M_PI);
  }

  mfield_free(mfield_workspace_p);
  msynth_free(msynth_workspace_p);

  return 0;
}
