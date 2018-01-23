lib_LTLIBRARIES = libpca3d.la

libpca3d_la_SOURCES = io.c tiegcm3d_alloc.c tiegcm3d_print.c tiegcm3d_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/track -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper

common_libs = $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -L/usr/lib64 -lsatdata -lindices ~/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lcdf

bin_PROGRAMS = print_modes print print_fft spectrogram stage2a stage3a

print_SOURCES = print.c
print_LDADD = libpca3d.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}

print_fft_SOURCES = print_fft.c
print_fft_LDADD = libpca3d.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}

print_modes_SOURCES = print_modes.c
print_modes_LDADD = libpca3d.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}

spectrogram_SOURCES = spectrogram.c
spectrogram_LDADD = libpca3d.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}

stage2a_SOURCES = stage2a.c
stage2a_CFLAGS = -fopenmp
stage2a_LDFLAGS = -fopenmp
stage2a_LDADD = libpca3d.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}

stage3a_SOURCES = stage3a.c
stage3a_LDADD = libpca3d.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm -lgsl -lgslcblas -lnetcdf ${common_libs}
