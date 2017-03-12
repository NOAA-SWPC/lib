lib_LTLIBRARIES = libpoltor.la

libpoltor_la_SOURCES = poltor.c poltor_shell.c poltor_synth.c

bin_PROGRAMS = compare_plp invert plot preproc

commonlibs = -L/home/palken/usr/lib -lflow -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread -lgfortran

compare_plp_SOURCES = compare_plp.c
compare_plp_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/lls/liblls.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs}

invert_SOURCES = invert_main.c
invert_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/lls/liblls.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la -lm ${commonlibs}

plot_SOURCES = plot.c
plot_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/lls/liblls.la $(top_builddir)/msynth/libmsynth.la  $(top_builddir)/common/libcommon.la -lm -lgfortran ${commonlibs}

preproc_SOURCES = preproc.c
preproc_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/common/libcommon.la $(top_builddir)/msynth/libmsynth.la $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/pomme -I$(top_builddir)/common -I$(top_builddir)/msynth -I$(top_builddir)/apex -I$(top_builddir)/curvefit -I$(top_builddir)/lls -I$(top_builddir)/track -I$(top_builddir)/euler -I$(top_builddir)/magdata
