lib_LTLIBRARIES = libcond.la

libcond_la_SOURCES = cond.c

AM_CFLAGS = -g -O2 -Wall -W

AM_CPPFLAGS = -I$(top_builddir)/iri -I$(top_builddir)/msis -I$(top_builddir)/common -I$(top_builddir)/msynth

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libcond.la $(top_builddir)/common/libcommon.la  $(top_builddir)/iri/libiri.la $(top_builddir)/msis/libmsis.la $(top_builddir)/msynth/libmsynth.la -lm /usr/lib64/libg2c.so.0 -lgsl -lgslcblas -lgfortran -lindices -L/home/palken/usr/lib
