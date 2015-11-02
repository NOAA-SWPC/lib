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
} mfield_test_coeff;

#endif /* INCLUDED_mfield_test_h */
