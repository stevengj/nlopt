AM_CPPFLAGS = -I$(top_srcdir)/cdirect -I$(top_srcdir)/direct -I$(top_srcdir)/stogo -I$(top_srcdir)/praxis -I$(top_srcdir)/lbfgs -I$(top_srcdir)/luksan -I$(top_srcdir)/crs -I$(top_srcdir)/mlsl -I$(top_srcdir)/mma -I$(top_srcdir)/cobyla -I$(top_srcdir)/newuoa -I$(top_srcdir)/neldermead -I$(top_srcdir)/auglag -I$(top_srcdir)/bobyqa -I$(top_srcdir)/isres -I$(top_srcdir)/slsqp -I$(top_srcdir)/esch -I$(top_srcdir)/util

include_HEADERS = nlopt.h nlopt.f nlopt.hpp
noinst_LTLIBRARIES = libapi.la
dist_man_MANS = nlopt.3

libapi_la_SOURCES = general.c options.c optimize.c deprecated.c	\
nlopt-internal.h nlopt.h f77api.c f77funcs.h f77funcs_.h

BUILT_SOURCES = nlopt.f nlopt.hpp
EXTRA_DIST = nlopt-in.hpp

if MAINTAINER_MODE
# convert constants to F77 parameter statements & C++

nlopt.f: nlopt.h
	rm -f $@
	(i=0; egrep 'NLOPT_[LG][DN]|NLOPT_AUGLAG|NLOPT_G_MLSL' $(srcdir)/nlopt.h | tr -d ' =0,' | while read n; do echo "      integer $$n"; echo "      parameter ($$n=$$i)"; i=`expr $$i + 1`; done; tail -n +`grep -n enum $(srcdir)/nlopt.h |cut -d: -f1 |tail -n 1` $(srcdir)/nlopt.h  | grep NLOPT | egrep -v 'EXTERN|DEPRECATED|MINF_MAX' | cut -d, -f1 | tr -d ' ' | perl -pe 's/([A-Za-z0-9_]+)=([-+0-9]+)/      integer \1\n      parameter (\1=\2)/') > $(srcdir)/$@

nlopt.hpp: nlopt.h nlopt-in.hpp
	rm -f $@
	(n=`grep -n GEN_ENUMS_HERE $(srcdir)/nlopt-in.hpp | cut -d: -f1`; head -n $$n $(srcdir)/nlopt-in.hpp;  echo "  enum algorithm {"; egrep 'NLOPT_[LG][DN]|NLOPT_AUGLAG|NLOPT_G_MLSL|NLOPT_NUM_ALGORITHMS' $(srcdir)/nlopt.h | sed 's/NLOPT_//g'; echo "  };"; echo "  enum result {"; egrep 'NLOPT_[A-Z_]* =' $(srcdir)/nlopt.h | egrep -v 'NLOPT_[LG][DN]' | sed 's/NLOPT_//g'; echo "  };"; tail -n +$$n $(srcdir)/nlopt-in.hpp) > $(srcdir)/$@

endif
