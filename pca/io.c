/*
 * io.c
 *
 * I/O routines for TIEGCM project
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "io.h"
#include "tiegcm.h"

static int matlab_dlmwrite_complex(FILE * fp, const gsl_matrix_complex * A);

int
pca_write_data(const char *filename, const size_t nmax, const size_t mmax, const tiegcm_data *data)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fwrite(&nmax, sizeof(size_t), 1, fp);
  fwrite(&mmax, sizeof(size_t), 1, fp);
  fwrite(&(data->nt), sizeof(size_t), 1, fp);
  fwrite(data->ut, sizeof(double), data->nt, fp);

  fclose(fp);

  return 0;
}

int
pca_read_data(const char *filename, size_t *nmax, size_t *mmax, size_t *nt, double *ut)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "pca_read_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fread(nmax, sizeof(size_t), 1, fp);
  fread(mmax, sizeof(size_t), 1, fp);

  if (nt && ut)
    {
      fread(nt, sizeof(size_t), 1, fp);
      fread(ut, sizeof(double), *nt, fp);
    }

  fclose(fp);

  return 0;
}

int
pca_write_vector(const char *filename, const gsl_vector * v)
{
  int s;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_vector: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fwrite(&(v->size), sizeof(size_t), 1, fp);
  s = gsl_vector_fwrite(fp, v);

  fclose(fp);

  return s;
}

gsl_vector *
pca_read_vector(const char *filename)
{
  FILE *fp;
  gsl_vector *v;
  size_t vsize;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "pca_read_vector: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  fread(&vsize, sizeof(size_t), 1, fp);

  v = gsl_vector_alloc(vsize);

  gsl_vector_fread(fp, v);

  fclose(fp);

  return v;
}

int
pca_write_matrix(const char *filename, const gsl_matrix * m)
{
  int s;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_matrix: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fwrite(&(m->size1), sizeof(size_t), 1, fp);
  fwrite(&(m->size2), sizeof(size_t), 1, fp);
  s = gsl_matrix_fwrite(fp, m);

  fclose(fp);

  return s;
}

gsl_matrix *
pca_read_matrix(const char *filename)
{
  FILE *fp;
  gsl_matrix *m;
  size_t size1, size2;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "pca_read_matrix: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  fread(&size1, sizeof(size_t), 1, fp);
  fread(&size2, sizeof(size_t), 1, fp);

  m = gsl_matrix_alloc(size1, size2);

  gsl_matrix_fread(fp, m);

  fclose(fp);

  return m;
}

int
pca_write_matrix_complex(const char *filename, const gsl_matrix_complex * m)
{
  int s;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_matrix: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fwrite(&(m->size1), sizeof(size_t), 1, fp);
  fwrite(&(m->size2), sizeof(size_t), 1, fp);
  s = gsl_matrix_complex_fwrite(fp, m);

  fclose(fp);

  return s;
}

gsl_matrix_complex *
pca_read_matrix_complex(const char *filename)
{
  FILE *fp;
  gsl_matrix_complex *m;
  size_t size1, size2;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "pca_read_matrix_complex: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  fread(&size1, sizeof(size_t), 1, fp);
  fread(&size2, sizeof(size_t), 1, fp);

  m = gsl_matrix_complex_alloc(size1, size2);

  gsl_matrix_complex_fread(fp, m);

  fclose(fp);

  return m;
}

int
pca_write_S(const char *filename, const size_t nmax, const size_t mmax,
            const double freq_cpd, const double window_size,
            const double window_shift, const gsl_vector *S)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_S: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fprintf(fp, "%% Singular values\n");
  fprintf(fp, "%% nmax: %zu\n", nmax);
  fprintf(fp, "%% mmax: %zu\n", mmax);
  fprintf(fp, "%% Frequency: %.6f [cpd]\n", freq_cpd);
  fprintf(fp, "%% window size:  %g [days]\n", window_size);
  fprintf(fp, "%% window slide: %g [days]\n", window_shift);
  gsl_vector_fprintf(fp, S, "%.12e");

  fclose(fp);

  return s;
}

int
pca_write_complex_U(const char *filename, const size_t nmax, const size_t mmax,
                    const double freq_cpd, const double window_size,
                    const double window_shift, const size_t nmodes, const gsl_matrix_complex *U)
{
  int s = 0;
  FILE *fp;
  gsl_matrix_complex_const_view m = gsl_matrix_complex_const_submatrix(U, 0, 0, U->size1, nmodes);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_complex_U: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fprintf(fp, "%% Left singular vectors\n");
  fprintf(fp, "%% nmax: %zu\n", nmax);
  fprintf(fp, "%% mmax: %zu\n", mmax);
  fprintf(fp, "%% number of modes (columns): %zu\n", nmodes);
  fprintf(fp, "%% Frequency: %.6f [cpd]\n", freq_cpd);
  fprintf(fp, "%% window size:  %g [days]\n", window_size);
  fprintf(fp, "%% window slide: %g [days]\n", window_shift);
  fprintf(fp, "%% Matlab read command: U = dlmread('%s',',',8,0);\n", filename);
  matlab_dlmwrite_complex(fp, &m.matrix);

  fclose(fp);

  return s;
}

int
pca_write_complex_V(const char *filename, const size_t nmax, const size_t mmax,
                    const double freq_cpd, const double window_size,
                    const double window_shift, const size_t nmodes, const gsl_matrix_complex *V)
{
  int s = 0;
  const size_t N = V->size1;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_write_complex_V: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fprintf(fp, "%% Right singular vectors\n");
  fprintf(fp, "%% nmax: %zu\n", nmax);
  fprintf(fp, "%% mmax: %zu\n", mmax);
  fprintf(fp, "%% number of time segments (rows and columns): %zu\n", N);
  fprintf(fp, "%% Frequency: %.6f [cpd]\n", freq_cpd);
  fprintf(fp, "%% window size:  %g [days]\n", window_size);
  fprintf(fp, "%% window slide: %g [days]\n", window_shift);
  fprintf(fp, "%% Matlab read command: V = dlmread('%s',',',8,0);\n", filename);
  matlab_dlmwrite_complex(fp, V);

  fclose(fp);

  return s;
}

static int
matlab_dlmwrite_complex(FILE * fp, const gsl_matrix_complex * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z = gsl_matrix_complex_get(A, i, j);
          double zr = GSL_REAL(z);
          double zi = GSL_IMAG(z);

          fprintf(fp, "%.12e%s%.12ei%s",
                  zr,
                  (zi < 0.0) ? "" : "+",
                  zi,
                  (j < N - 1) ? "," : "\n");
        }
    }

  return 0;
}
