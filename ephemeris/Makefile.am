lib_LTLIBRARIES = libephemeris.la
libephemeris_la_SOURCES = eph.c eph_data.c hermite.c

AM_CPPFLAGS = -I$(top_builddir)/common

bin_PROGRAMS = print

print_SOURCES = print.c
print_LDADD = libephemeris.la $(top_builddir)/common/libcommon.la -lm -lcdf -lindices -lsatdata -lz -lgsl -lgslcblas

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libephemeris.la $(top_builddir)/common/libcommon.la -lm -lcdf -lindices -lsatdata -lz -lgsl -lgslcblas
