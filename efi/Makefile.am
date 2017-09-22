lib_LTLIBRARIES = libefi.la
libefi_la_SOURCES = efi.c

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libefi.la -lcommon -lm -lcdf -lgsl -lgslcblas -lsatdata -lindices

AM_CPPFLAGS =
