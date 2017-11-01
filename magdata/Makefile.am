noinst_LTLIBRARIES = libmagdata.la

AM_CPPFLAGS = -I$(top_builddir)/track

libmagdata_la_SOURCES = magdata.c magdata_list.c preproc.c

check_PROGRAMS = test print

print_SOURCES = print.c
print_LDADD = libmagdata.la $(top_builddir)/track/libtrack.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib

test_SOURCES = test.c
test_LDADD = libmagdata.la $(top_builddir)/track/libtrack.la -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib
