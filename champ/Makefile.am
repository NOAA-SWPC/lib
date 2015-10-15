bin_PROGRAMS = print print_alt print_lt print_plp stage1

commonlibs = -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lm -lcdf

print_SOURCES = print.c
print_LDADD = $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

print_alt_SOURCES = print_alt.c
print_alt_LDADD = $(top_builddir)/common/libcommon.la -lm ${commonlibs}

print_lt_SOURCES = print_lt.c
print_lt_LDADD = $(top_builddir)/common/libcommon.la -lm ${commonlibs}

print_plp_SOURCES = print_plp.c
print_plp_LDADD = $(top_builddir)/common/libcommon.la -lm ${commonlibs}

stage1_SOURCES = stage1.c
stage1_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/chaos/libchaos.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs} -lgfortran -lmwmclmcrrt $(top_builddir)/chaos/libchaosext.so

AM_CPPFLAGS = -I$(top_builddir)/pomme -I$(top_builddir)/common -I$(top_builddir)/msynth -I$(top_builddir)/apex -I$(top_builddir)/track -I$(top_builddir)/euler -I$(top_builddir)/chaos
