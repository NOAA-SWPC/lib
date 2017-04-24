lib_LTLIBRARIES = libpca.la

libpca_la_SOURCES = io.c pca.c tiegcm_alloc.c tiegcm_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/apex -I$(top_builddir)/common -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/track -I$(top_builddir)/msynth -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper

common_libs = $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la $(top_builddir)/chaos/libchaos.la -L/usr/lib64 -lgsl /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lsatdata -lindices -lcdf -lmwmclmcrrt $(top_builddir)/chaos/libchaosext.so

bin_PROGRAMS = plot_maps print stage1 stage2 stage2b stage3 stage3b

plot_maps_SOURCES = plot_maps.c
plot_maps_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la -lm -lgsl ${common_libs}

print_SOURCES = print.c
print_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas -lnetcdf ${common_libs}

stage1_SOURCES = stage1.c
stage1_LDADD = libpca.la $(top_builddir)/apex/libapex.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/lls/liblls.la $(top_builddir)/green/libgreen.la  $(top_builddir)/track/libtrack.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lsatdata -lindices -lcdf -lgsl -lnetcdf ${common_libs}

stage2_SOURCES = stage2.c
stage2_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl ${common_libs}

stage2b_SOURCES = stage2b.c
stage2b_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl ${common_libs}

stage3_SOURCES = stage3.c
stage3_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl ${common_libs}

stage3b_SOURCES = stage3b.c
stage3b_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl ${common_libs}
