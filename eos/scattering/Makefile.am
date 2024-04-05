CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeosscattering.la
libeosscattering_la_SOURCES = \
    ee-to-ccbar.hh \
    ee-to-ccbar-hkrvd2024.cc ee-to-ccbar-hkrvd2024.hh \
    observables.cc observables.hh
libeosscattering_la_LIBADD = \
    $(top_builddir)/eos/utils/libeosutils.la \
    -lgslcblas \
    -lgsl

include_eos_scatteringdir = $(includedir)/eos/scattering
include_eos_scattering_HEADERS = \
    ee-to-ccbar.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
    export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
    ee-to-ccbar-hkrvd2024_TEST

LDADD = \
    $(top_builddir)/test/libeostest.la \
    libeosscattering.la \
    $(top_builddir)/eos/utils/libeosutils.la \
    $(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
ee_to_ccbar_hkrvd2024_TEST_SOURCES = ee-to-ccbar-hkrvd2024_TEST.cc
