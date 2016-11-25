lib_LTLIBRARIES = libsecs.la

libsecs_la_SOURCES = secs1d.c

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/track -I$(top_builddir)/msynth -I$(top_builddir)/green -I$(top_builddir)/lapack_wrapper -I$(top_builddir)/mageq

check_PROGRAMS = test main

test_SOURCES = test.c
test_LDADD = libsecs.la $(top_builddir)/track/libtrack.la $(top_builddir)/green/libgreen.la  $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread -lgfortran

main_SOURCES = main.c
main_LDADD = libsecs.la $(top_builddir)/track/libtrack.la $(top_builddir)/green/libgreen.la  $(top_builddir)/common/libcommon.la $(top_builddir)/mageq/libmageq.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread -lgfortran
