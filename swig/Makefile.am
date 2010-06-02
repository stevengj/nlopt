EXTRA_DIST = nlopt.i nlopt-guile.i

BUILT_SOURCES = nlopt-guile.cpp nlopt-renames.i

HDR = $(top_srcdir)/api/nlopt.hpp

if MAINTAINER_MODE

nlopt-guile.cpp: nlopt.i nlopt-enum-renames.i nlopt-guile.i $(top_builddir)/api/nlopt.hpp
	swig -I$(top_builddir)/api -outdir $(builddir) -c++ -guile -scmstub -o $@ nlopt.i

nlopt-enum-renames.i: $(top_srcdir)/api/nlopt.h
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; egrep 'NLOPT_[LG][DN]|NLOPT_NUM_ALGORITHMS' $(top_srcdir)/api/nlopt.h | sed 's/NLOPT_//g' |tr -d ' ' |tr '/' ',' |tr '=' ',' |cut -d, -f1 |while read name; do echo "%rename(NLOPT_$$name) nlopt::$$name;"; done; egrep 'NLOPT_[A-Z_]* =' $(top_srcdir)/api/nlopt.h | egrep -v 'NLOPT_[LG][DN]' | sed 's/NLOPT_//g' |tr -d ' ' |cut -d'=' -f1 | while read name; do echo "%rename(NLOPT_$$name) nlopt::$$name;"; done) > $@

endif
