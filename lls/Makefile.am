noinst_LTLIBRARIES = liblls.la

AM_CPPFLAGS =

liblls_la_SOURCES = lls.c lls_complex.c lls_lapack.c tsqr.c

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = liblls.la -lcommon -lm -llapacke -llapack ~/usr/lib/libgsl.a -lptcblas -lptf77blas -latlas -lpthread -lgfortran
