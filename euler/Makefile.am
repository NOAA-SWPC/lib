noinst_LTLIBRARIES = libeuler.la

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/msynth -I$(top_builddir)/magdata -I$(top_builddir)/track -I$(top_builddir)/levmar

libeuler_la_SOURCES = euler.c euler_calc.c

bin_PROGRAMS = main preproc

main_SOURCES = main.c
main_LDADD = libeuler.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/levmar/liblevmar.la $(top_builddir)/common/libcommon.la $(top_builddir)/track/libtrack.la -lm -lgsl -llevmar_lib -llapack -lf77blas -lcblas -latlas -lgfortran -lsatdata -lindices -lcdf -L/home/palken/usr/lib

preproc_SOURCES = preproc.c
preproc_LDADD = libeuler.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/common/libcommon.la $(top_builddir)/track/libtrack.la $(top_builddir)/msynth/libmsynth.la -lm -lsatdata -lindices -lcdf -L/home/palken/usr/lib -lgsl -lgslcblas

check_PROGRAMS = test test_quat
test_SOURCES = test.c
test_LDADD = libeuler.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib

test_quat_SOURCES = test_quat.c
test_quat_LDADD = libeuler.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/common/libcommon.la $(top_builddir)/track/libtrack.la -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib
