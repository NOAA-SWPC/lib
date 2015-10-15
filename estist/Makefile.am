lib_LTLIBRARIES = libestist_calc.la
libestist_calc_la_SOURCES = estist_calc.c weidelt_dst.f

F77 = gfortran
AM_CPPFLAGS = -I$(top_builddir)/common

bin_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libestist_calc.la $(top_builddir)/common/libcommon.la -lm -lindices -lgfortran -lgsl -lgslcblas

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libestist_calc.la $(top_builddir)/common/libcommon.la -lm -lindices -lgfortran -lgsl -lgslcblas
