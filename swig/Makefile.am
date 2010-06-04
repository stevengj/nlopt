SWIG_SRC = nlopt.i nlopt-exceptions.i nlopt-enum-renames.i 
EXTRA_DIST = $(SWIG_SRC) nlopt-guile.i nlopt-python.i nlopt.scm.in

BUILT_SOURCES = nlopt-guile.cpp nlopt-enum-renames.i nlopt.scm.in

HDR = $(top_builddir)/api/nlopt.hpp

##############################################################################
# Guile wrapper

libnlopt@NLOPT_SUFFIX@_guile_la_SOURCES = nlopt-guile.cpp
libnlopt@NLOPT_SUFFIX@_guile_la_LIBADD = ../libnlopt@NLOPT_SUFFIX@.la
libnlopt@NLOPT_SUFFIX@_guile_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@
libnlopt@NLOPT_SUFFIX@_guile_la_CPPFLAGS = $(GUILE_CPPFLAGS) -I$(top_builddir)/api

guiledir = $(GUILE_SITEDIR)
if WITH_GUILE
guile_DATA = nlopt.scm
guilelib = libnlopt@NLOPT_SUFFIX@_guile.la
endif



##############################################################################

lib_LTLIBRARIES = $(guilelib)

##############################################################################
# Generate wrappers via SWIG:

if MAINTAINER_MODE

nlopt-guile.cpp nlopt.scm.in: $(SWIG_SRC) nlopt-guile.i $(HDR)
	swig -I$(top_builddir)/api -outdir $(builddir) -c++ -guile -scmstub -o $@ nlopt.i
	rm -f nlopt.scm.in
	mv nlopt.scm nlopt.scm.in

nlopt-python.cpp nlopt.py: $(SWIG_SRC) nlopt-python.i $(HDR)
	swig -I$(top_builddir)/api -outdir $(builddir) -c++ -python -o $@ nlopt.i

nlopt-enum-renames.i: $(top_srcdir)/api/nlopt.h
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; egrep 'NLOPT_[LG][DN]|NLOPT_NUM_ALGORITHMS' $(top_srcdir)/api/nlopt.h | sed 's/NLOPT_//g' |tr -d ' ' |tr '/' ',' |tr '=' ',' |cut -d, -f1 |while read name; do echo "%rename(NLOPT_$$name) nlopt::$$name;"; done; egrep 'NLOPT_[A-Z_]* =' $(top_srcdir)/api/nlopt.h | egrep -v 'NLOPT_[LG][DN]' | sed 's/NLOPT_//g' |tr -d ' ' |cut -d'=' -f1 | while read name; do echo "%rename(NLOPT_$$name) nlopt::$$name;"; done) > $@

endif
