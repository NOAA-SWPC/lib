/*
 * mfield_test.h
 */

#ifndef INCLUDED_mfield_test_h
#define INCLUDED_mfield_test_h

typedef struct
{
  size_t n;
  int m;
  double gnm;   /* static */
  double dgnm;  /* SV */
  double ddgnm; /* SA */
} mfield_synth_coeff;

/*
 * Prototypes
 */

int mfield_synth_g(gsl_vector * g, mfield_workspace * w);
int mfield_synth_replace(mfield_workspace *w);

#endif /* INCLUDED_mfield_test_h */
