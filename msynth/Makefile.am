lib_LTLIBRARIES = libmsynth.la

libmsynth_la_SOURCES = msynth.c msynth_green.c msynth_grid.c msynth_arnaud.c msynth_bggm.c msynth_swarm.c msynth_crust.c msynth_emm.c msynth_igrf.c msynth_ipgp.c msynth_pomme.c msynth_tgcm.c msynth_wmm.c

check_PROGRAMS = test grid calc_point extrapolate_epoch model_avg print_map print_smatrix print_spectrum print_test_values write_wmm write_chris

calc_point_SOURCES = calc_point.c
calc_point_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

grid_SOURCES = grid.c
grid_CFLAGS = -fopenmp
grid_LDFLAGS = -fopenmp
grid_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas -lfftw3

extrapolate_epoch_SOURCES = extrapolate_epoch.c
extrapolate_epoch_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

print_map_SOURCES = print_map.c
print_map_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas -lfftw3

print_smatrix_SOURCES = print_smatrix.c
print_smatrix_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

print_spectrum_SOURCES = print_spectrum.c
print_spectrum_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

model_avg_SOURCES = model_avg.c
model_avg_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

print_test_values_SOURCES = print_test_values.c
print_test_values_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

write_wmm_SOURCES = write_wmm.c
write_wmm_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

write_chris_SOURCES = write_chris.c
write_chris_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

test_SOURCES = test.c
test_LDADD = libmsynth.la $(top_builddir)/common/libcommon.la -lm ~/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/common
