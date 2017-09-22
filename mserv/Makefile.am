AM_CPPFLAGS = -I$(top_builddir)/grobs

check_PROGRAMS = storm_rms

storm_rms_SOURCES = storm_rms.c
storm_rms_LDADD = $(top_builddir)/grobs/libgrobs.la -lcommon -lm -lgsl -lgslcblas
