AM_CPPFLAGS = -I$(top_srcdir)/api -I$(top_srcdir)/util

noinst_PROGRAMS = testopt

testopt_SOURCES = testfuncs.c testfuncs.h testopt.cpp test_std.py
testopt_LDADD = $(top_builddir)/libnlopt@NLOPT_SUFFIX@.la
