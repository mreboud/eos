CLEANFILES = \
	*~ \
	references_TEST
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = constraints parameters utils optimize statistics form-factors rare-b-decays b-decays .

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeos.la
libeos_la_SOURCES = \
	constraint.cc constraint.hh \
	decays.hh \
	init.cc \
	observable.cc observable.hh observable-fwd.hh observable-impl.hh \
	reference.cc reference.hh \
	signal-pdf.cc signal-pdf.hh
libeos_la_CXXFLAGS = $(AM_CXXFLAGS) \
	-DEOS_DATADIR='"$(datadir)"' \
	$(GSL_CXXFLAGS) \
	$(YAMLCPP_CXXFLAGS)
libeos_la_LDFLAGS = $(GSL_LDFLAGS) $(YAMLCPP_LDFLAGS)
libeos_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la

include_eosdir = $(includedir)/eos
include_eos_HEADERS = \
	constraint.hh \
	decays.hh \
	observable.hh \
	reference.hh \
	signal-pdf.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
	export EOS_TESTS_REFERENCES="$(top_srcdir)/eos/";

TESTS = \
	constraint_TEST \
	reference_TEST

LDADD = \
	$(top_builddir)/test/libeostest.a \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = \
	constraint_TEST \
	reference_TEST

constraint_TEST_SOURCES = constraint_TEST.cc
constraint_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
constraint_TEST_LDADD = $(LDADD) -lyaml-cpp

reference_TEST_SOURCES = reference_TEST.cc
reference_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
reference_TEST_LDADD = $(LDADD) -lyaml-cpp

pkgdata_DATA = references.yaml
EXTRA_DIST = \
	references.yaml
