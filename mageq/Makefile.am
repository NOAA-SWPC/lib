noinst_LTLIBRARIES = libmageq.la

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/msynth

libmageq_la_SOURCES = mageq.c

check_PROGRAMS = test

test_SOURCES = test.c

test_LDADD = libmageq.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas -lindices -L/home/palken/usr/lib
