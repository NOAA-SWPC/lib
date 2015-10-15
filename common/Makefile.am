bin_PROGRAMS = test test_zenith test_sunrise time2fday test_lunar test_lunar2 test_lunar_m2 test_doy test_interp test_omni time2str time2str2 time2str3 fday2str epoch2str test_geo test_julian test_estist2 test_ut1 doy2md

doy2md_SOURCES = doy2md.c
test_SOURCES = test.c
test_zenith_SOURCES = test_zenith.c
test_sunrise_SOURCES = test_sunrise.c
time2fday_SOURCES = time2fday.c
time2str_SOURCES = time2str.c
time2str2_SOURCES = time2str2.c
time2str3_SOURCES = time2str3.c
fday2str_SOURCES = fday2str.c
epoch2str_SOURCES = epoch2str.c
test_lunar_SOURCES = test_lunar.c
test_lunar2_SOURCES = test_lunar2.c
test_lunar_m2_SOURCES = test_lunar_m2.c
test_doy_SOURCES = test_doy.c
test_interp_SOURCES = test_interp.c
test_omni_SOURCES = test_omni.c
test_geo_SOURCES = test_geo.c
test_julian_SOURCES = test_julian.c
test_ut1_SOURCES = test_ut1.c

lib_LTLIBRARIES = libcommon.la

libcommon_la_SOURCES = \
  bin.c                \
	bin2d.c              \
	bin3d.c              \
	bsearch.c            \
	ccorr.c              \
	coord.c              \
  common.c             \
  eci.c                \
	ellipsoid.c          \
  ema.c                \
  estist2.c            \
  fit.c                \
  fit2.c               \
  geo.c                \
  histogramNd.c        \
  interp.c             \
  julian.c             \
  lse.c                \
  lapack_lls.c         \
  lunar.c              \
  oct.c                \
  omni.c               \
  plmbar.c             \
	rm.c                 \
  solarpos.c           \
  solpos.c             \
  tcalc.c              \
  tcalc2.c             \
  ut1.c

AM_CFLAGS = -std=c99 -g -O2 -Wall -W -D_XOPEN_SOURCE
AM_CPPFLAGS =

test_common_libs = -lm -lgsl -lgslcblas

doy2md_LDADD = libcommon.la ${test_common_libs}
test_LDADD = libcommon.la ${test_common_libs}
test_zenith_LDADD = libcommon.la ${test_common_libs}
test_estist2_LDADD = libcommon.la ${test_common_libs}
test_sunrise_LDADD = libcommon.la ${test_common_libs}
time2fday_LDADD = libcommon.la ${test_common_libs}
time2str_LDADD = libcommon.la ${test_common_libs}
time2str2_LDADD = libcommon.la ${test_common_libs}
time2str3_LDADD = libcommon.la ${test_common_libs}
fday2str_LDADD = libcommon.la ${test_common_libs}
epoch2str_LDADD = libcommon.la ${test_common_libs} -lsatdata -lcdf -lindices
test_lunar_LDADD = libcommon.la ${test_common_libs}
test_lunar2_LDADD = libcommon.la ${test_common_libs}
test_lunar_m2_LDADD = libcommon.la ${test_common_libs}
test_doy_LDADD = libcommon.la ${test_common_libs}
test_interp_LDADD = libcommon.la ${test_common_libs}
test_omni_LDADD = libcommon.la ${test_common_libs}
test_geo_LDADD = libcommon.la ${test_common_libs}
test_julian_LDADD = libcommon.la ${test_common_libs}
test_ut1_LDADD = libcommon.la ${test_common_libs}
