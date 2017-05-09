bin_PROGRAMS = print stage0 stage1 stage2

commonlibs = -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lm -lcdf

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/ephemeris -I$(top_builddir)/chaos -I$(top_builddir)/msynth -I$(top_builddir)/track -I$(top_builddir)/estist -I$(top_builddir)/pomme

stage0_SOURCES = stage0.c
stage0_LDADD = $(top_builddir)/ephemeris/libephemeris.la $(top_builddir)/common/libcommon.la -lm -lcdf ${commonlibs} -lz

stage1_SOURCES = stage1.c
stage1_LDFLAGS = -fopenmp
stage1_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/chaos/libchaos.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm -lgfortran ${commonlibs} -lmwmclmcrrt $(abs_top_builddir)/chaos/libchaosext.so

stage2_SOURCES = magcal.c stage2.c
stage2_LDFLAGS = -fopenmp
stage2_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/chaos/libchaos.la $(top_builddir)/ephemeris/libephemeris.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm -lgfortran ${commonlibs} -lz

print_SOURCES = print.c
print_LDADD = $(top_builddir)/common/libcommon.la -lm -lcdf ${commonlibs}
