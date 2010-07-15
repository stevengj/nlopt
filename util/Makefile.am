AM_CPPFLAGS = -I$(top_srcdir)/api

noinst_LTLIBRARIES = libutil.la
libutil_la_SOURCES = mt19937ar.c sobolseq.c soboldata.h timer.c stop.c nlopt-util.h redblack.c redblack.h qsort_r.c rescale.c

noinst_PROGRAMS = redblack_test
redblack_test_SOURCES = redblack_test.c
redblack_test_LDADD = libutil.la
