lib_LTLIBRARIES = libapex.la
libapex_la_SOURCES = apexf77.f apexsh.f90 makeapexsh.f90 magfld.f apex.c

F77 = gfortran
AM_FFLAGS = -fno-second-underscore

AM_CPPFLAGS = -I$(top_builddir)/common

check_PROGRAMS = print_qd test

test_SOURCES = test.c
test_LDADD = libapex.la $(top_builddir)/common/libcommon.la -lm -lgfortran -lgsl -lgslcblas

print_qd_SOURCES = print_qd.c
print_qd_LDADD = libapex.la $(top_builddir)/common/libcommon.la -lm -lgfortran -lgsl -lgslcblas
