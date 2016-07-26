lib_LTLIBRARIES = libefi.la
libefi_la_SOURCES = efi.c

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libefi.la $(top_builddir)/common/libcommon.la -lm -lcdf -lgsl -lgslcblas -lsatdata

AM_CPPFLAGS = -I$(top_builddir)/common
