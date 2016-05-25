lib_LTLIBRARIES = libsuperlu.la

libsuperlu_la_SOURCES = superlu_mt.c

AM_CPPFLAGS = -I/home/palken/usr/include

AM_CFLAGS = -fopenmp

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libsuperlu.la -lm /home/palken/usr/lib/libgsl.a -lsuperlu_mt_OPENMP -lcblas -lf77blas -latlas
