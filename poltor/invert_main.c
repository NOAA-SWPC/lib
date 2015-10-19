/*
 * invert_main.c
 *
 * Invert vector magnetic field measurements for toroidal currents
 * flowing in a thin shell at 110km altitude
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <flow/flow.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "oct.h"
#include "track.h"

#include "poltor.h"

#define POLTOR_FIT_X               1 /* fit X component */
#define POLTOR_FIT_Y               1 /* fit Y component */
#define POLTOR_FIT_Z               1 /* fit Z component */
#define POLTOR_FIT_DX              1 /* fit DX component */
#define POLTOR_FIT_DY              1 /* fit DY component */
#define POLTOR_FIT_DZ              1 /* fit DZ component */

/* define to use QD coordinates in spherical harmonic transforms */
#define POLTOR_QD_HARMONICS        1

/*
set_flags()
  The preproc program will have calculated the magdata struct,
but won't set the individual X/Y/Z fitting flags, since that should
be done here. However it would set the MAGDATA_FLG_DZ_NS flag if
gradient information is available for each point.

So go through data, and based on the FIT_xxx settings, re-assign
magdata flags
*/

int
set_flags(magdata *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      size_t flags = 0;

      /* only fit X/Y at mid/low latitudes */
      if (data->qdlat[i] >= -60.0 && data->qdlat[i] <= 60.0)
        {
#if POLTOR_FIT_X
          flags |= MAGDATA_FLG_X;
#endif
#if POLTOR_FIT_Y
          flags |= MAGDATA_FLG_Y;
#endif
        }

#if POLTOR_FIT_Z
      flags |= MAGDATA_FLG_Z;
#endif

      if (data->flags[i] & MAGDATA_FLG_DZ_NS)
        {
          /* vector along-track gradient is available for this point */

          /* only fit X/Y at mid/low latitudes */
          if (data->qdlat[i] >= -60.0 && data->qdlat[i] <= 60.0)
            {
#if POLTOR_FIT_DX
              flags |= MAGDATA_FLG_DX_NS;
#endif
#if POLTOR_FIT_DY
              flags |= MAGDATA_FLG_DY_NS;
#endif
            }

          /* fit B_z at all latitudes */
#if POLTOR_FIT_DZ
          flags |= MAGDATA_FLG_DZ_NS;
#endif
        }

      data->flags[i] = flags;
    }

  return s;
} /* set_flags() */

int
print_correlation(const char *filename, poltor_workspace *w)
{
  int s = 0;
  const size_t p = w->p;
  const size_t nmax = GSL_MIN(w->nmax_int, w->nmax_sh);
  const size_t mmax = GSL_MIN(w->mmax_int, w->mmax_sh);
  gsl_matrix_complex *B;
  size_t n;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_correlation: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  B = gsl_matrix_complex_calloc(p, p);

  /* compute correlation matrix */
  s = lls_complex_correlation(B, w->lls_workspace_p);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax);
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t idx1 = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          size_t idx2 = poltor_jnmidx(0, n, m, w);
          gsl_complex z = gsl_matrix_complex_get(B, idx1, idx2);

          fprintf(fp, "%3zu %3d %e\n",
                  n,
                  m,
                  gsl_complex_abs(z));
        }

      fprintf(fp, "\n\n");
    }

  gsl_matrix_complex_free(B);

  fclose(fp);

  return s;
} /* print_correlation() */

int
print_coefficients(poltor_workspace *w)
{
#if POLTOR_SYNTH_DATA

  /* print and check synthetic coefficients */
  poltor_synth_print(w);

#else

  const double b = w->b;
  const double R = w->R;
  const size_t nmax_int = GSL_MIN(3, w->nmax_int);
  const size_t mmax_int = GSL_MIN(3, w->mmax_int);
  const size_t nmax_ext = GSL_MIN(2, w->nmax_ext);
  const size_t mmax_ext = GSL_MIN(2, w->mmax_ext);
  const size_t nmax_sh = GSL_MIN(2, w->nmax_sh);
  const size_t mmax_sh = GSL_MIN(2, w->mmax_sh);
  const size_t nmax_tor = GSL_MIN(3, w->nmax_tor);
  const size_t mmax_tor = GSL_MIN(3, w->mmax_tor);
  size_t n;

  /* print internal poloidal coefficients */
  fprintf(stderr, "Internal poloidal coefficients:\n");
  for (n = 1; n <= nmax_int; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_int);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double gnm = GSL_REAL(coef);
          double qnm = -(2.0*n + 1.0) / (double)n * gnm *
                       pow(R / b, n + 3.0);

          fprintf(stderr, "g(%2zu,%2d) = %12g q(%2zu,%2d) = %12g [nT]\n",
                  n, m, gnm,
                  n, m, qnm);
        }
    }

  /* print external poloidal coefficients */
  fprintf(stderr, "External poloidal coefficients:\n");
  for (n = 1; n <= nmax_ext; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_ext);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double knm = GSL_REAL(coef);

          fprintf(stderr, "k(%2zu,%2d) = %12g [nT]\n",
                  n,
                  m,
                  knm);
        }
    }

  /* print shell poloidal coefficients */
  {
    char buf[2048];
    char *bufptr = buf;
    int offset;
    size_t j;

    fprintf(stderr, "Shell poloidal coefficients:\n");
    for (n = 1; n <= nmax_sh; ++n)
      {
        int ni = (int) GSL_MIN(n, mmax_sh);
        int m;

        /*
         * only need to print positive m coefficients since
         * c_{n,-m} = c_{nm}
         */
        for (m = 0; m <= ni; ++m)
          {
            bufptr = buf;
            for (j = 0; j <= w->shell_J; ++j)
              {
                size_t cidx = poltor_jnmidx(j, n, m, w);
                gsl_complex coef = poltor_get(cidx, w);
                double qnm = GSL_REAL(coef);

                sprintf(bufptr, "%12g%s %n",
                        qnm,
                        (j < w->shell_J) ? "," : "",
                        &offset);
                bufptr += offset;
              }

            fprintf(stderr, "q(%2zu,%2d) = %s [nT]\n", n, m, buf);
          }
      }
  }

  /* print toroidal coefficients */
  fprintf(stderr, "Shell toroidal coefficients:\n");
  for (n = 1; n <= nmax_tor; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_tor);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double phinm = GSL_REAL(coef);

          fprintf(stderr, "phi(%2zu,%2d) = %12g [nT]\n",
                  n,
                  m,
                  phinm);
        }
    }

#endif /* POLTOR_SYNTH_DATA */

  return 0;
} /* print_coefficients() */

int
print_Lcurve(const char *filename, poltor_workspace *w)
{
  int s = 0;
  FILE *fp;
  const size_t p = w->p;
  double rnorm, Lnorm;
  gsl_vector_complex_view v = gsl_vector_complex_subvector(w->rhs, 0, p);
  size_t i;

  fp = fopen(filename, "a");
  if (!fp)
    {
      fprintf(stderr, "print_Lcurve: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  /* construct A and b, and calculate chi^2 = ||b - A c||^2 */
  poltor_build_ls(0, w);
  rnorm = sqrt(w->chisq);

  /* compute v = L c; L is stored in w->L by poltor_solve() */
  for (i = 0; i < p; ++i)
    {
      gsl_complex ci = gsl_vector_complex_get(w->c, i);
      double li = gsl_vector_get(w->L, i);
      gsl_complex val = gsl_complex_mul_real(ci, li);

      gsl_vector_complex_set(&v.vector, i, val);
    }

  /* compute || L c || */
  Lnorm = gsl_blas_dznrm2(&v.vector);

  fprintf(fp, "%.12e %.12e %.6e %.6e %.6e\n",
          log(rnorm),
          log(Lnorm),
          w->alpha_int,
          w->alpha_sh,
          w->alpha_tor);

  printcv_octave(w->residuals, "r");
  printcv_octave(w->c, "c");
  printv_octave(w->L, "L");

  fclose(fp);

  return s;
} /* print_Lcurve() */

int
print_residuals(const char *filename, poltor_workspace *w)
{
  size_t i, j;
  FILE *fp;
  magdata *data = w->data;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: X observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: DX observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: DY observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: DZ observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: X residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: DX residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: DY residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: DZ residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal DX model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal DY model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal DZ model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal DX model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal DY model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal DZ model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal DX model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal DY model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal DZ model (nT)\n", i++);
  fprintf(fp, "# Field %zu: vector along-track gradient available (1 or 0)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double B_obs[4], B_model[3], B_int[4], B_sh[4], B_ext[4], B_tor[4], B_res[3];
      double dB_model[3], B_int2[4], B_sh2[4], B_ext2[4], B_tor2[4];
      double dB_obs[4] = { 0.0, 0.0, 0.0, 0.0 };
      double dB_res[3] = { 0.0, 0.0, 0.0 };
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lt = get_localtime(unix_time, data->phi[i]);
      double theta = poltor_theta(i, w);

      /* store observation vector */
      magdata_residual(i, B_obs, data);

      /* compute individual magnetic field models at this point */
      poltor_eval_B_all(data->r[i], theta, data->phi[i], B_int, B_ext, B_sh, B_tor, w);

      if (data->flags[i] & MAGDATA_FLG_DZ_NS)
        {
          double theta_ns = poltor_theta_ns(i, w);

          magdata_residual_dB_ns(i, dB_obs, data);

          /* evaluate model at along-track point */
          poltor_eval_B_all(data->r_ns[i], theta_ns, data->phi_ns[i], B_int2, B_ext2, B_sh2, B_tor2, w);
        }

      /* compute total model and residual vector */
      for (j = 0; j < 3; ++j)
        {
          B_model[j] = B_int[j] + B_ext[j] + B_sh[j] + B_tor[j];
          B_res[j] = B_obs[j] - B_model[j];

          if (data->flags[i] & MAGDATA_FLG_DZ_NS)
            {
              dB_model[j] = B_int2[j] + B_ext2[j] + B_sh2[j] + B_tor2[j] - B_model[j];
              dB_res[j] = dB_obs[j] - dB_model[j];
            }
        }

      fprintf(fp, "%f %.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d\n",
              satdata_epoch2year(data->t[i]),
              lt,
              data->r[i] - w->R,
              data->phi[i] * 180.0 / M_PI,
              90.0 - data->theta[i] * 180.0 / M_PI,
              data->qdlat[i],
              B_obs[0],
              B_obs[1],
              B_obs[2],
              data->satdir[i] * dB_obs[0],
              data->satdir[i] * dB_obs[1],
              data->satdir[i] * dB_obs[2],
              B_res[0],
              B_res[1],
              B_res[2],
              data->satdir[i] * dB_res[0],
              data->satdir[i] * dB_res[1],
              data->satdir[i] * dB_res[2],
              B_int[0],
              B_int[1],
              B_int[2],
              B_sh[0],
              B_sh[1],
              B_sh[2],
              B_tor[0],
              B_tor[1],
              B_tor[2],
              data->satdir[i] * (B_int2[0] - B_int[0]),
              data->satdir[i] * (B_int2[1] - B_int[1]),
              data->satdir[i] * (B_int2[2] - B_int[2]),
              data->satdir[i] * (B_sh2[0] - B_sh[0]),
              data->satdir[i] * (B_sh2[1] - B_sh[1]),
              data->satdir[i] * (B_sh2[2] - B_sh[2]),
              data->satdir[i] * (B_tor2[0] - B_tor[0]),
              data->satdir[i] * (B_tor2[1] - B_tor[1]),
              data->satdir[i] * (B_tor2[2] - B_tor[2]),
              data->flags[i] & MAGDATA_FLG_DZ_NS ? 1 : 0);
    }

  fclose(fp);

  return 0;
} /* print_residuals() */

int
print_chisq(const char *filename, poltor_workspace *w)
{
  FILE *fp;
  poltor_parameters *params = &(w->params);

  fp = fopen(filename, "a");
  if (!fp)
    {
      fprintf(stderr, "print_chisq: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fprintf(fp, "%zu %zu %zu %zu %zu %zu %e\n",
          params->nmax_int,
          params->mmax_int,
          params->nmax_tor,
          params->mmax_tor,
          params->nmax_sh,
          params->mmax_sh,
          w->chisq / w->dof);

  fclose(fp);

  return 0;
} /* print_chisq() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options] input1.dat input2.dat ...\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --nmax_int | -n nmax                - internal nmax\n");
  fprintf(stderr, "\t --mmax_int | -m mmax                - internal mmax\n");
  fprintf(stderr, "\t --nmax_tor | -a nmax                - toroidal nmax\n");
  fprintf(stderr, "\t --mmax_tor | -b mmax                - toroidal mmax\n");
  fprintf(stderr, "\t --nmax_sh | -e nmax                 - shell nmax\n");
  fprintf(stderr, "\t --mmax_sh | -f mmax                 - shell mmax\n");
  fprintf(stderr, "\t --alpha_int | -c alpha_int          - damping parameter for internal coefficients\n");
  fprintf(stderr, "\t --alpha_sh | -d alpha_sh            - damping parameter for shell poloidal coefficients\n");
  fprintf(stderr, "\t --alpha_tor | -j alpha_tor          - damping parameter for shell toroidal coefficients\n");
  fprintf(stderr, "\t --residual_file | -r file           - residual output file\n");
  fprintf(stderr, "\t --coef_file | -o file               - coefficient output file\n");
  fprintf(stderr, "\t --chisq_file | -p file              - chi^2 output file (for L-curve analysis)\n");
  fprintf(stderr, "\t --universal_time | -t decimal UT    - Universal time (hours)\n");
  fprintf(stderr, "\t --lls_file | -l lls_file            - LS system file (matrix and rhs)\n");
  fprintf(stderr, "\t --lcurve_file | -k lcurve_file      - output file for L-curve data\n");
  fprintf(stderr, "\t --maxit | -q maxit                  - number of robust iterations\n");
  fprintf(stderr, "\t --print_data | -u                   - print data file\n");
}

int
main(int argc, char *argv[])
{
  size_t nmax_int = 60;
  size_t mmax_int = 6;
  size_t nmax_ext = 0;
  size_t mmax_ext = 0;
  size_t nmax_sh = 60;
  size_t mmax_sh = 5;
  size_t nmax_tor = 60;
  size_t mmax_tor = 5;
  double alpha_int = 1.0;
  double alpha_sh = 1.0;
  double alpha_tor = 1.0;
  size_t robust_maxit = 5;
  const double R = R_EARTH_KM;
  const double b = R + 110.0;   /* radius of internal current shell (Sq+EEJ) */
  const double d = R + 350.0;   /* radius of current shell for gravity/diamag */
  double universal_time = 11.0; /* UT in hours for data selection */
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *spectrum_file = "poltor.s";
  char *corr_file = "corr.dat";
  char *residual_file = NULL;
  char *output_file = NULL;
  char *chisq_file = NULL;
  char *lls_file = NULL;
  char *Lcurve_file = NULL;
  magdata *mdata = NULL;
  poltor_workspace *poltor_p;
  poltor_parameters params;
  struct timeval tv0, tv1;
  int print_data = 0;

#if POLTOR_SYNTH_DATA
  nmax_int = 30;
  mmax_int = 10;
  nmax_ext = 2;
  mmax_ext = 2;
  nmax_sh = 20;
  mmax_sh = 10;
  nmax_tor = 30;
  mmax_tor = 10;
#endif

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "nmax_int", required_argument, NULL, 'n' },
          { "mmax_int", required_argument, NULL, 'm' },
          { "nmax_tor", required_argument, NULL, 'a' },
          { "mmax_tor", required_argument, NULL, 'b' },
          { "nmax_sh", required_argument, NULL, 'e' },
          { "mmax_sh", required_argument, NULL, 'f' },
          { "nmax_ext", required_argument, NULL, 'g' },
          { "mmax_ext", required_argument, NULL, 'h' },
          { "residual_file", required_argument, NULL, 'r' },
          { "output_file", required_argument, NULL, 'o' },
          { "chisq_file", required_argument, NULL, 'p' },
          { "universal_time", required_argument, NULL, 't' },
          { "lls_file", required_argument, NULL, 'l' },
          { "lcurve_file", required_argument, NULL, 'k' },
          { "alpha_int", required_argument, NULL, 'c' },
          { "alpha_sh", required_argument, NULL, 'd' },
          { "alpha_tor", required_argument, NULL, 'j' },
          { "maxit", required_argument, NULL, 'q' },
          { "print_data", no_argument, NULL, 'u' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:d:e:f:g:h:j:k:l:m:n:o:p:q:r:t:u", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'n':
            nmax_int = (size_t) atoi(optarg);
            break;

          case 'm':
            mmax_int = (size_t) atoi(optarg);
            break;

          case 'a':
            nmax_tor = (size_t) atoi(optarg);
            break;

          case 'b':
            mmax_tor = (size_t) atoi(optarg);
            break;

          case 'e':
            nmax_sh = (size_t) atoi(optarg);
            break;

          case 'f':
            mmax_sh = (size_t) atoi(optarg);
            break;

          case 'g':
            nmax_ext = (size_t) atoi(optarg);
            break;

          case 'h':
            mmax_ext = (size_t) atoi(optarg);
            break;

          case 'c':
            alpha_int = atof(optarg);
            break;

          case 'd':
            alpha_sh = atof(optarg);
            break;

          case 'j':
            alpha_tor = atof(optarg);
            break;

          case 'r':
            residual_file = optarg;
            break;

          case 'k':
            Lcurve_file = optarg;
            break;

          case 'o':
            output_file = optarg;
            break;

          case 't':
            universal_time = atof(optarg);
            break;

          case 'p':
            chisq_file = optarg;
            break;

          case 'l':
            lls_file = optarg;
            break;

          case 'q':
            robust_maxit = (size_t) atoi(optarg);
            break;

          case 'u':
            print_data = 1;
            break;

          default:
            break;
        }
    }

  while (optind < argc)
    {
      fprintf(stderr, "main: reading %s...", argv[optind]);
      gettimeofday(&tv0, NULL);
      mdata = magdata_read(argv[optind], mdata);
      gettimeofday(&tv1, NULL);

      if (!mdata)
        exit(1);

      fprintf(stderr, "done (%zu data total, %g seconds)\n",
              mdata->n, time_diff(tv0, tv1));

      ++optind;
    }

  if (!mdata)
    {
      print_help(argv);
      exit(1);
    }

  mmax_int = GSL_MIN(mmax_int, nmax_int);
  mmax_ext = GSL_MIN(mmax_ext, nmax_ext);
  mmax_sh = GSL_MIN(mmax_sh, nmax_sh);
  mmax_tor = GSL_MIN(mmax_tor, nmax_tor);

  fprintf(stderr, "main: universal time = %.1f\n", universal_time);

  fprintf(stderr, "main: nmax_int  = %zu\n", nmax_int);
  fprintf(stderr, "main: mmax_int  = %zu\n", mmax_int);
  fprintf(stderr, "main: nmax_ext  = %zu\n", nmax_ext);
  fprintf(stderr, "main: mmax_ext  = %zu\n", mmax_ext);
  fprintf(stderr, "main: nmax_sh   = %zu\n", nmax_sh);
  fprintf(stderr, "main: mmax_sh   = %zu\n", mmax_sh);
  fprintf(stderr, "main: nmax_tor  = %zu\n", nmax_tor);
  fprintf(stderr, "main: mmax_tor  = %zu\n", mmax_tor);
  fprintf(stderr, "main: alpha_int = %g\n", alpha_int);
  fprintf(stderr, "main: alpha_sh  = %g\n", alpha_sh);
  fprintf(stderr, "main: alpha_tor = %g\n", alpha_tor);

  if (residual_file)
    fprintf(stderr, "main: residual file = %s\n", residual_file);

  if (Lcurve_file)
    fprintf(stderr, "main: L-curve file  = %s\n", Lcurve_file);

  /*
   * re-compute flags for fitting components / gradient, etc;
   * must be called before magdata_init()
   */
  set_flags(mdata);

  fprintf(stderr, "main: initializing spatial weighting histogram...");
  gettimeofday(&tv0, NULL);
  magdata_init(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* re-compute weights, nvec, nres based on flags update */
  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if POLTOR_SYNTH_DATA
  fprintf(stderr, "main: setting unit spatial weights...");
  magdata_unit_weights(mdata);
  fprintf(stderr, "done\n");
#endif

  fprintf(stderr, "main: print_data = %d\n", print_data);
  if (print_data)
    {
      fprintf(stderr, "main: writing data to %s...", data_file);
      magdata_print(data_file, mdata);
      fprintf(stderr, "done\n");

      fprintf(stderr, "main: writing data map to %s...", datamap_file);
      magdata_map(datamap_file, mdata);
      fprintf(stderr, "done\n");
    }

  fprintf(stderr, "main: satellite rmin = %.1f (%.1f) [km]\n",
          mdata->rmin, mdata->rmin - mdata->R);
  fprintf(stderr, "main: satellite rmax = %.1f (%.1f) [km]\n",
          mdata->rmax, mdata->rmax - mdata->R);

  params.R = R;
  params.b = b;
  params.d = d;
  params.rmin = GSL_MAX(mdata->rmin, mdata->R + 250.0);
  params.rmax = GSL_MIN(mdata->rmax, mdata->R + 450.0);
  params.nmax_int = nmax_int;
  params.mmax_int = mmax_int;
  params.nmax_ext = nmax_ext;
  params.mmax_ext = mmax_ext;
  params.nmax_sh = nmax_sh;
  params.mmax_sh = mmax_sh;
  params.nmax_tor = nmax_tor;
  params.mmax_tor = mmax_tor;
  params.shell_J = 0;
  params.data = mdata;
  params.alpha_int = alpha_int;
  params.alpha_sh = alpha_sh;
  params.alpha_tor = alpha_tor;

#if POLTOR_QD_HARMONICS
  params.flags = POLTOR_FLG_QD_HARMONICS;
#else
  params.flags = 0;
#endif

  poltor_p = poltor_alloc(&params);

  fprintf(stderr, "main: poltor rmin = %.1f (%.1f) [km]\n",
          params.rmin, params.rmin - mdata->R);
  fprintf(stderr, "main: poltor rmax = %.1f (%.1f) [km]\n",
          params.rmax, params.rmax - mdata->R);

#if POLTOR_SYNTH_DATA
  fprintf(stderr, "main: replacing with synthetic data...");
  gettimeofday(&tv0, NULL);
  poltor_synth(poltor_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

  if (lls_file)
    {
      /* use previously computed LS system from file */
      fprintf(stderr, "main: loading LS system from %s...", lls_file);
      lls_complex_load(lls_file, poltor_p->lls_workspace_p);
      fprintf(stderr, "done\n");

      /* solve LS system */
      poltor_solve(poltor_p);
    }
  else
    {
      size_t maxiter = robust_maxit;
      size_t iter = 0;
      char buf[2048];

#if POLTOR_SYNTH_DATA
      maxiter = 1;
#endif

      while (iter++ < maxiter)
        {
          fprintf(stderr, "main: ROBUST ITERATION %zu/%zu\n", iter, maxiter);

          /* build LS system */
          poltor_calc(poltor_p);

          /* solve LS system */
          poltor_solve(poltor_p);

          sprintf(buf, "%s.iter%zu", spectrum_file, iter);
          fprintf(stderr, "main: printing spectrum to %s...", buf);
          poltor_print_spectrum(buf, poltor_p);
          fprintf(stderr, "done\n");
        }
    }

  print_coefficients(poltor_p);

  fprintf(stderr, "main: printing correlation data to %s...", corr_file);
  print_correlation(corr_file, poltor_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: printing spectrum to %s...", spectrum_file);
  poltor_print_spectrum(spectrum_file, poltor_p);
  fprintf(stderr, "done\n");

  if (Lcurve_file)
    {
      fprintf(stderr, "main: writing L-curve data to %s...", Lcurve_file);
      print_Lcurve(Lcurve_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (output_file)
    {
      fprintf(stderr, "main: writing output coefficients to %s...", output_file);
      poltor_write(output_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (residual_file)
    {
      fprintf(stderr, "main: printing residuals to %s...", residual_file);
      print_residuals(residual_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (chisq_file)
    {
      fprintf(stderr, "main: printing chisq/dof to %s...", chisq_file);
      print_chisq(chisq_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
  poltor_free(poltor_p);

  return 0;
} /* main() */
