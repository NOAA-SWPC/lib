lib_LTLIBRARIES = libjulia.la
libjulia_la_SOURCES = julia.c jicmag.c

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libjulia.la $(top_builddir)/common/libcommon.la -lm -lcdf -lgsl -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/common
