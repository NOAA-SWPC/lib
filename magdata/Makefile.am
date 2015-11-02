noinst_LTLIBRARIES = libmagdata.la

AM_CPPFLAGS = -I$(top_builddir)/common -I$(top_builddir)/track -I$(top_builddir)/msynth

libmagdata_la_SOURCES = magdata.c

bin_PROGRAMS = print

print_SOURCES = print.c
print_LDADD = libmagdata.la $(top_builddir)/common/libcommon.la $(top_builddir)/track/libtrack.la -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libmagdata.la $(top_builddir)/common/libcommon.la $(top_builddir)/track/libtrack.la -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib
