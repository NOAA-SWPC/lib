noinst_LTLIBRARIES = libmagdata.la

AM_CPPFLAGS = -I$(top_builddir)/track -I$(top_builddir)/pomme -I$(top_builddir)/estist

libmagdata_la_SOURCES = magdata.c

check_PROGRAMS = mag2swarm test print

print_SOURCES = print.c
print_LDADD = libmagdata.la $(top_builddir)/track/libtrack.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib

test_SOURCES = test.c
test_LDADD = libmagdata.la $(top_builddir)/track/libtrack.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib

mag2swarm_SOURCES = mag2swarm.c
mag2swarm_CFLAGS = -fopenmp
mag2swarm_LDFLAGS = -fopenmp
mag2swarm_LDADD = libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -lgfortran -L/home/palken/usr/lib
