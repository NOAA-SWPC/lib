lib_LTLIBRARIES = libgrobs.la

libgrobs_la_SOURCES = grobs.c iaga.c wamnet.c

AM_CPPFLAGS =

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgrobs.la -lcommon -lm -lgsl -lgslcblas
