lib_LTLIBRARIES = libapex.la
libapex_la_SOURCES = apexf77.f apexsh.f90 makeapexsh.f90 magfld.f apex.c

F77 = gfortran
AM_FFLAGS = -fno-second-underscore

AM_CPPFLAGS = -I$(top_builddir)/common

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libapex.la $(top_builddir)/common/libcommon.la -lm -lgfortran -lgsl -lgslcblas
