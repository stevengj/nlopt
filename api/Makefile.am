AM_CPPFLAGS = -I$(top_srcdir)/cdirect -I$(top_srcdir)/direct -I$(top_srcdir)/stogo -I$(top_srcdir)/praxis -I$(top_srcdir)/lbfgs -I$(top_srcdir)/luksan -I$(top_srcdir)/crs -I$(top_srcdir)/mlsl -I$(top_srcdir)/mma -I$(top_srcdir)/cobyla -I$(top_srcdir)/newuoa -I$(top_srcdir)/neldermead -I$(top_srcdir)/auglag -I$(top_srcdir)/util

include_HEADERS = nlopt.h nlopt.f
noinst_LTLIBRARIES = libapi.la
dist_man_MANS = nlopt_minimize.3 nlopt_minimize_constrained.3

libapi_la_SOURCES = nlopt.c nlopt.h f77api.c f77funcs.h

BUILT_SOURCES = nlopt.f

if MAINTAINER_MODE
# convert constants to F77 parameter statements

nlopt.f: nlopt.h
	rm -f $@
	(i=0; egrep 'NLOPT_[LG][DN]' $(srcdir)/nlopt.h | tr -d ' =0,' | while read n; do echo "      integer $$n"; echo "      parameter ($$n=$$i)"; i=`expr $$i + 1`; done; tail -n +`grep -n enum $(srcdir)/nlopt.h |cut -d: -f1 |tail -n 1` $(srcdir)/nlopt.h  | grep NLOPT | grep -v NLOPT_EXTERN | cut -d, -f1 | tr -d ' ' | perl -pe 's/([A-Za-z0-9_]+)=([-+0-9]+)/      integer \1\n      parameter (\1=\2)/') > $@

endif
