lib_LTLIBRARIES = libgreen.la

libgreen_la_SOURCES = current.c green.c

AM_CPPFLAGS =

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgreen.la -lcommon -lm -lgsl -lgslcblas
