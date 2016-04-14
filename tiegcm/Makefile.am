lib_LTLIBRARIES = libtiegcm.la

libtiegcm_la_SOURCES = io.c tiegcm_alloc.c tiegcm_read.c

AM_CPPFLAGS = -I$(top_builddir)/apex -I$(top_builddir)/common -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/track -I$(top_builddir)/msynth -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper

bin_PROGRAMS = print stage1 stage2 stage2b stage3 stage3b

print_SOURCES = print.c
print_LDADD = libtiegcm.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas -lnetcdf

stage1_SOURCES = stage1.c
stage1_LDADD = libtiegcm.la $(top_builddir)/apex/libapex.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/lls/liblls.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/green/libgreen.la  $(top_builddir)/common/libcommon.la $(top_builddir)/track/libtrack.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lsatdata -lindices -lcdf -lgsl -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread -lnetcdf -lgfortran

stage2_SOURCES = stage2.c
stage2_LDADD = libtiegcm.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl -llapacke -llapack -lptcblas -latlas -lpthread -lfftw3

stage2b_SOURCES = stage2b.c
stage2b_LDADD = libtiegcm.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl -llapacke -llapack -lptcblas -latlas -lpthread -lfftw3

stage3_SOURCES = stage3.c
stage3_LDADD = libtiegcm.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl -llapacke -llapack -lptcblas -latlas -lpthread -lfftw3

stage3b_SOURCES = stage3b.c
stage3b_LDADD = libtiegcm.la $(top_builddir)/green/libgreen.la $(top_builddir)/common/libcommon.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lgsl -llapacke -llapack -lptcblas -latlas -lpthread -lfftw3
