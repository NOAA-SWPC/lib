lib_LTLIBRARIES = libpomme.la

libpomme_la_SOURCES = geod2geoc.c gsm2geo.c sm2geo.c plmbar.c magfdz.c int-ext.c pomme.c pomme_wrap.c

AM_CPPFLAGS = -I$(top_builddir)/estist

check_PROGRAMS = test test_ext
test_SOURCES = test.c
test_LDADD = libpomme.la $(top_builddir)/estist/libestist_calc.la -lindices -lcommon -lm -lgfortran -lgsl -lgslcblas

test_ext_SOURCES = test_ext.c
test_ext_LDADD = libpomme.la $(top_builddir)/estist/libestist_calc.la -lindices -lcommon -lm -lgsl -lgslcblas -lgfortran
