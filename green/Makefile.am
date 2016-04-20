lib_LTLIBRARIES = libgreen.la

libgreen_la_SOURCES = green.c

AM_CPPFLAGS = -I$(top_builddir)/common

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgreen.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas
