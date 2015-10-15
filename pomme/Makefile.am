lib_LTLIBRARIES = libpomme.la

pmodeldir = $(top_builddir)/pomme8.1/POMME8

libpomme_la_SOURCES = $(pmodeldir)/geod2geoc.c $(pmodeldir)/gsm2geo.c $(pmodeldir)/sm2geo.c $(pmodeldir)/plmbar.c $(pmodeldir)/magfdz.c $(pmodeldir)/int-ext.c $(pmodeldir)/pomme.c pomme_wrap.c

AM_CPPFLAGS = -I$(pmodeldir) -I$(top_builddir)/common

check_PROGRAMS = test test2 test_ext
test_SOURCES = test.c
test_LDADD = libpomme.la $(top_builddir)/common/libcommon.la -lindices -lm -lgsl -lgslcblas

test_ext_SOURCES = test_ext.c
test_ext_LDADD = libpomme.la $(top_builddir)/common/libcommon.la -lindices -lm -lgsl -lgslcblas

test2_SOURCES = test2.c
test2_LDADD = libpomme.la $(top_builddir)/common/libcommon.la -lm
