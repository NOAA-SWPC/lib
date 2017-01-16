lib_LTLIBRARIES = libmsynth.la

libmsynth_la_SOURCES = msynth.c msynth_green.c msynth_arnaud.c msynth_chaos.c msynth_crust.c msynth_emm.c msynth_igrf.c msynth_pomme.c msynth_tgcm.c msynth_wmm.c

check_PROGRAMS = test calc_point extrapolate_epoch model_avg print_map print_spectrum print_test_values write_wmm

calc_point_SOURCES = calc_point.c
calc_point_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

extrapolate_epoch_SOURCES = extrapolate_epoch.c
extrapolate_epoch_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

print_map_SOURCES = print_map.c
print_map_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

print_spectrum_SOURCES = print_spectrum.c
print_spectrum_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

model_avg_SOURCES = model_avg.c
model_avg_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

print_test_values_SOURCES = print_test_values.c
print_test_values_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

write_wmm_SOURCES = write_wmm.c
write_wmm_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

test_SOURCES = test.c
test_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/common
