lib_LTLIBRARIES = libmag.la

libmag_la_SOURCES = \
	coord.c           \
	geod2geoc.c       \
	green.c           \
	inverteef.c       \
	lls.c             \
  log.c             \
	mag.c             \
	mag_eej.c         \
	mag_log.c         \
	mag_sqfilt.c      \
	pde.c             \
	sigma.c           \
	swarmeef.c

bin_PROGRAMS = corr main swarmeef_print compare

main_SOURCES = main.c
swarmeef_print_SOURCES = swarmeef_print.c
compare_SOURCES = compare.c
corr_SOURCES = corr.c

AM_CFLAGS = -Wall -W

AM_CPPFLAGS = -I/home/palken/usr/include -I$(top_builddir)/common -I$(top_builddir)/mageq -I$(top_builddir)/pomme -I$(top_builddir)/cond -I$(top_builddir)/apex -I$(top_builddir)/hwm -I$(top_builddir)/iri -I$(top_builddir)/msis -I$(top_builddir)/lis -I$(top_builddir)/msynth -I$(top_builddir)/julia -I$(top_builddir)/curvefit -I$(top_builddir)/efi -I$(top_builddir)/track

main_LDADD = libmag.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/common/libcommon.la $(top_builddir)/apex/libapex.la $(top_builddir)/hwm/libhwm.la $(top_builddir)/mageq/libmageq.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/cond/libcond.la $(top_builddir)/iri/libiri.la $(top_builddir)/msis/libmsis.la $(top_builddir)/lis/libmylis.la $(top_builddir)/track/libtrack.la -lm -lsatdata -lcdf -llapack -lindices -llis -lgsl -lptcblas -lptf77blas -latlas -lgfortran -lpthread

swarmeef_print_LDADD = libmag.la $(top_builddir)/common/libcommon.la -lsatdata -lcdf -lindices -lgsl -lgslcblas -lm

compare_LDADD = $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/julia/libjulia.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/efi/libefi.la -lsatdata -lcdf -lindices /home/palken/usr/lib/libgsl.a -lgslcblas

corr_LDADD = $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lsatdata -lcdf -lindices /home/palken/usr/lib/libgsl.a -lgslcblas
