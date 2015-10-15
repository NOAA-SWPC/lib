lib_LTLIBRARIES = liblevmar.la

AM_CPPFLAGS = -I$(top_builddir)/common

liblevmar_la_SOURCES = levmar.c

check_PROGRAMS = test test2
test_SOURCES = test.c
test_LDADD = liblevmar.la $(top_builddir)/common/libcommon.la ~/usr/lib/libgsl.a -llevmar_lib -llapack -lf77blas -lcblas -latlas -lm -lgfortran

test2_SOURCES = test2.c
test2_LDADD = liblevmar.la $(top_builddir)/common/libcommon.la ~/usr/lib/libgsl.a -llevmar_lib -llapack -lf77blas -lcblas -latlas -lm -lgfortran
