noinst_LTLIBRARIES = liblls.la

AM_CPPFLAGS = -I$(top_builddir)/common

liblls_la_SOURCES = lls.c lls_complex.c lls_lapack.c

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = liblls.la $(top_builddir)/common/libcommon.la -lm -llapacke -llapack ~/usr/lib/libgsl.a -lptcblas -lptf77blas -latlas -lpthread -lgfortran
