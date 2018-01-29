noinst_LTLIBRARIES = libmageq.la

AM_CPPFLAGS =

libmageq_la_SOURCES = mageq.c magpole.c

check_PROGRAMS = test pole

test_SOURCES = test.c
test_LDADD = libmageq.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lindices -L/home/palken/usr/lib

pole_SOURCES = pole.c
pole_LDADD = libmageq.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lindices -L/home/palken/usr/lib
