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

int
pca_write_data(const char *filename, const size_t nmax, const size_t mmax)
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

  fclose(fp);

  return 0;
}

int
pca_read_data(const char *filename, size_t *nmax, size_t *mmax)
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