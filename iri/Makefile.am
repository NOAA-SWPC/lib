lib_LTLIBRARIES = libiri.la

libiri_la_SOURCES = iri.c irisub.for irifun.for iritec.for iridreg.for iriflip.for igrf.for cira.for

F77 = gfortran
AM_FFLAGS = -fno-second-underscore -O2

AM_CFLAGS = -std=c99 -D_XOPEN_SOURCE -g -O2 -Wall -W

AM_CPPFLAGS =

check_PROGRAMS = test plot_map

test_SOURCES = test.c
test_LDADD = libiri.la -lm -lgsl -lgslcblas -lgfortran -lindices

plot_map_SOURCES = plot_map.c
plot_map_LDADD = libiri.la -lapex -lcommon -lm -lgsl -lgslcblas -lgfortran -lindices
