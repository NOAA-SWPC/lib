lib_LTLIBRARIES = libpomme.la

libpomme_la_SOURCES = geod2geoc.c gsm2geo.c sm2geo.c plmbar.c magfdz.c int-ext.c pomme.c pomme_wrap.c

AM_CPPFLAGS = -I$(top_builddir)/common

check_PROGRAMS = test test2 test_ext
test_SOURCES = test.c
test_LDADD = libpomme.la $(top_builddir)/common/libcommon.la -lindices -lm -lgsl -lgslcblas

test_ext_SOURCES = test_ext.c
test_ext_LDADD = libpomme.la $(top_builddir)/common/libcommon.la -lindices -lm -lgsl -lgslcblas

test2_SOURCES = test2.c
test2_LDADD = libpomme.la $(top_builddir)/common/libcommon.la -lm
