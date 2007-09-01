AM_CPPFLAGS = -I$(top_srcdir)/util -I$(top_srcdir)/api

noinst_LTLIBRARIES = libcdirect.la
libcdirect_la_SOURCES = cdirect.c hybrid.c cdirect.h redblack.c redblack.h

noinst_PROGRAMS = redblack_test
redblack_test_SOURCES = redblack_test.c
redblack_test_LDADD = libcdirect.la

EXTRA_DIST = 
