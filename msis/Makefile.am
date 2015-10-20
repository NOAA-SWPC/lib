lib_LTLIBRARIES = libmsis.la

libmsis_la_SOURCES = msis.c gtd.f

F77 = g77
AM_FFLAGS = -fno-second-underscore -O2

AM_CFLAGS = -std=c99 -g -O2 -Wall -W -D_XOPEN_SOURCE

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libmsis.la -lm /usr/lib64/libg2c.so.0 -lindices
