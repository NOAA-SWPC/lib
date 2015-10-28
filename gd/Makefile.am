bin_PROGRAMS = necorr

necorr_SOURCES = necorr.c
necorr_LDADD = $(top_builddir)/peak/libpeak.la $(top_builddir)/track/libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/euler/libeuler.la -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/msynth -I$(top_builddir)/track -I$(top_builddir)/peak
