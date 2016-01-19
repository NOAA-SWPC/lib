lib_LTLIBRARIES = libtiegcm.la

libtiegcm_la_SOURCES = tiegcm_alloc.c tiegcm_read.c

AM_CPPFLAGS = -I$(top_builddir)/common

bin_PROGRAMS = print
print_SOURCES = print.c
print_LDADD = libtiegcm.la $(top_builddir)/common/libcommon.la -lm -lgsl -lgslcblas -lnetcdf
