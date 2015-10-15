lib_LTLIBRARIES = libcurvefit.la

libcurvefit_la_SOURCES = curvefit.c gaussfit.c

check_PROGRAMS = test
test_SOURCES = test.c
test_LDADD = libcurvefit.la ~/usr/lib/libgsl.a -lgslcblas -lm
