lib_LTLIBRARIES = libmsis.la

libmsis_la_SOURCES = msis.c nrlmsise00_sub.for

#F77 = g77
AM_FFLAGS = -fno-second-underscore -O2

AM_CFLAGS = -std=c99 -g -O2 -Wall -W -D_XOPEN_SOURCE

check_PROGRAMS = test nrldriver
test_SOURCES = test.c
test_LDADD = libmsis.la -lm -lindices -lgfortran

nrldriver_SOURCES = nrlmsise00_driver.for
nrldriver_LDADD = libmsis.la -lm -lgfortran
