lib_LTLIBRARIES = libestist_calc.la
libestist_calc_la_SOURCES = estist_calc.c weidelt_dst.f

F77 = gfortran
AM_CPPFLAGS =

bin_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libestist_calc.la -lcommon -lm -lindices -llapack /home/palken/usr/lib/libf77blas.a /home/palken/usr/lib/libatlas.a -lgfortran -lgsl -lgslcblas

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libestist_calc.la -lcommon -lm -lindices -lgfortran -lgsl -lgslcblas
