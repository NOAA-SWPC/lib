lib_LTLIBRARIES = libmag.la

libmag_la_SOURCES =        \
  cfg.c                    \
	coord.c                  \
	geod2geoc.c              \
	inverteef.c              \
  log.c                    \
	mag.c                    \
	mag_grad.c               \
	mag_eej.c                \
	mag_log.c                \
	mag_log_grad.c           \
	mag_sqfilt_scalar.c      \
	mag_sqfilt_vector.c      \
	pde.c                    \
	sigma.c                  \
	swarmeef.c

check_PROGRAMS = corr main swarmeef_print compare

main_SOURCES = main.c
swarmeef_print_SOURCES = swarmeef_print.c
compare_SOURCES = compare.c
corr_SOURCES = corr.c

AM_CFLAGS = -Wall -W

AM_CPPFLAGS = -I/home/palken/usr/include -I$(top_builddir)/mageq -I$(top_builddir)/pomme -I$(top_builddir)/cond -I$(top_builddir)/hwm -I$(top_builddir)/iri -I$(top_builddir)/msis -I$(top_builddir)/lis -I$(top_builddir)/julia -I$(top_builddir)/curvefit -I$(top_builddir)/efi -I$(top_builddir)/track -I$(top_builddir)/superlu -I$(top_builddir)/estist -I$(top_builddir)/magfit -I$(top_builddir)/green

main_LDADD = libmag.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/hwm/libhwm.la $(top_builddir)/mageq/libmageq.la $(top_builddir)/cond/libcond.la $(top_builddir)/iri/libiri.la $(top_builddir)/msis/libmsis.la $(top_builddir)/lis/libmylis.la $(top_builddir)/superlu/libsuperlu.la $(top_builddir)/track/libtrack.la $(top_builddir)/estist/libestist_calc.la -lapex -lcommon -lmsynth -lm -lsatdata -lcdf -llapack -lsuperlu_mt_OPENMP -lindices -llis -lgsl -lptcblas -lptf77blas -latlas -lgfortran -lpthread -lconfig

swarmeef_print_LDADD = libmag.la -lcommon -lsatdata -lcdf -lindices -lgsl -lgslcblas -lm

compare_LDADD = $(top_builddir)/julia/libjulia.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/efi/libefi.la -lcommon -lmsynth -lm -lsatdata -lcdf -lindices /home/palken/usr/lib/libgsl.a -lgslcblas

corr_LDADD = -lcommon -lmsynth -lm -lsatdata -lcdf -lindices /home/palken/usr/lib/libgsl.a -lgslcblas
