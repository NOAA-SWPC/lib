lib_LTLIBRARIES = libmylis.la

libmylis_la_SOURCES = lisw.c

AM_CFLAGS = -fopenmp

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libmylis.la -lm /home/palken/usr/lib/libgsl.a -llis -lcblas -latlas
