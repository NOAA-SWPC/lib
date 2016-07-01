lib_LTLIBRARIES = libtrack.la

libtrack_la_SOURCES = track.c track_filter.c track_flag.c track_offsets.c track_synth.c track_weight.c

bin_PROGRAMS = print

print_SOURCES = print.c
print_LDADD = libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/msynth -I$(top_builddir)/apex -I$(top_builddir)/pomme -I$(top_builddir)/curvefit -I$(top_builddir)/chaos -I$(top_builddir)/estist
