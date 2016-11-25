lib_LTLIBRARIES = libgrobs.la

libgrobs_la_SOURCES = grobs.c iaga.c wamnet.c

AM_CPPFLAGS = -I$(top_builddir)/common

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgrobs.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas
