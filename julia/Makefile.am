lib_LTLIBRARIES = libjulia.la
libjulia_la_SOURCES = julia.c jicmag.c

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libjulia.la -lcommon -lm -lcdf -lgsl -lgslcblas

AM_CPPFLAGS =
