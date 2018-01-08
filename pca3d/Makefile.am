lib_LTLIBRARIES = libpca3d.la

libpca3d_la_SOURCES = io.c tiegcm3d_alloc.c tiegcm3d_read.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/track -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper

common_libs = $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -L/usr/lib64 -lsatdata -lindices ~/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lcdf

bin_PROGRAMS = print

print_SOURCES = print.c
print_LDADD = libpca3d.la $(top_builddir)/green/libgreen.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}
