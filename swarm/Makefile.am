#bin_PROGRAMS = calc_crf2s calc_euler crust dynamo plot_tracks print print_all print_crust print_alt residuals stage1 swarmcdf2fday testcal
bin_PROGRAMS = compare print print_all print_efi print_lp print_lt print_crust print_alt residuals stage1 swarmcdf2fday

commonlibs = -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas

stage1_SOURCES = stage1.c
stage1_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/chaos/libchaos.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm -lgfortran ${commonlibs} -lmwmclmcrrt $(top_builddir)/chaos/libchaosext.so

#plot_tracks_SOURCES = plot_tracks.c
#plot_tracks_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/common/libcommon.la -lm ${commonlibs}

#calc_crf2s_SOURCES = calc_crf2s.c
#calc_crf2s_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/common/libcommon.la ${commonlibs}

#calc_euler_SOURCES = calc_euler.c
#calc_euler_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/common/libcommon.la ${commonlibs}

#dynamo_SOURCES = dynamo.c
#dynamo_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la -lm -lgfortran ${commonlibs}

#crust_SOURCES = crust.c
#crust_LDADD = $(top_builddir)/euler/libeuler.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs}

compare_SOURCES = compare.c
compare_LDADD = $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

print_SOURCES = print.c
print_LDADD = $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs}

print_lt_SOURCES = print_lt.c
print_lt_LDADD = $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs}

print_all_SOURCES = print_all.c
print_all_LDADD = $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

print_alt_SOURCES = print_alt.c
print_alt_LDADD = $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la -lm ${commonlibs}

print_efi_SOURCES = print_efi.c
print_efi_LDADD = $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la -lm ${commonlibs} -lgfortran

print_lp_SOURCES = print_lp.c
print_lp_LDADD = $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la -lm ${commonlibs} -lgfortran

print_crust_SOURCES = print_crust.c
print_crust_LDADD = $(top_builddir)/apex/libapex.la $(top_builddir)/common/libcommon.la -lm ${commonlibs}

residuals_SOURCES = residuals.c
residuals_LDFLAGS = -fopenmp
residuals_LDADD = $(top_builddir)/pomme/libpomme.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la ${commonlibs}

#testcal_SOURCES = testcal.c magcal.c
#testcal_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs}

swarmcdf2fday_SOURCES = swarmcdf2fday.c
swarmcdf2fday_LDADD = $(top_builddir)/common/libcommon.la -lm ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/pomme -I$(top_builddir)/common -I$(top_builddir)/msynth -I$(top_builddir)/apex -I$(top_builddir)/curvefit -I$(top_builddir)/track -I$(top_builddir)/euler -I$(top_builddir)/chaos
