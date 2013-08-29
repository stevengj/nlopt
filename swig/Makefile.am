SWIG_SRC = nlopt.i nlopt-exceptions.i nlopt-enum-renames.i 
EXTRA_DIST = $(SWIG_SRC) nlopt-guile.i nlopt-python.i nlopt.scm.in nlopt.py numpy.i

BUILT_SOURCES = nlopt-guile.cpp nlopt-python.cpp nlopt-enum-renames.i nlopt.scm.in

HDR = $(top_srcdir)/api/nlopt.hpp

##############################################################################
# Guile wrapper

libnlopt@NLOPT_SUFFIX@_guile_la_SOURCES = nlopt-guile.cpp
libnlopt@NLOPT_SUFFIX@_guile_la_LIBADD = ../libnlopt@NLOPT_SUFFIX@.la
libnlopt@NLOPT_SUFFIX@_guile_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@
libnlopt@NLOPT_SUFFIX@_guile_la_CPPFLAGS = $(GUILE_CPPFLAGS) -I$(top_srcdir)/api

guiledir = $(GUILE_INSTALL_DIR)
if WITH_GUILE
guile_DATA = nlopt.scm
guilelib = libnlopt@NLOPT_SUFFIX@_guile.la
endif

##############################################################################
# Python wrapper

_nlopt_la_SOURCES = nlopt-python.cpp
_nlopt_la_LIBADD = ../libnlopt@NLOPT_SUFFIX@.la
_nlopt_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_nlopt_la_CPPFLAGS = $(PYTHON_INCLUDES) -I$(top_srcdir)/api

if WITH_PYTHON
python_PYTHON = nlopt.py
pyexec_LTLIBRARIES = _nlopt.la
endif

##############################################################################


lib_LTLIBRARIES = $(guilelib)

##############################################################################
# Generate wrappers via SWIG:

if MAINTAINER_MODE

nlopt-guile.cpp nlopt.scm.in: $(SWIG_SRC) nlopt-guile.i $(HDR)
	swig -I$(top_srcdir)/api -outdir $(builddir) -c++ -guile -scmstub -o $@ $(srcdir)/nlopt.i
	rm -f nlopt.scm.in
	mv nlopt.scm nlopt.scm.in

nlopt-python.cpp nlopt.py: $(SWIG_SRC) nlopt-python.i numpy.i $(HDR)
	swig -I$(top_srcdir)/api -outdir $(builddir) -c++ -python -o $@ $(srcdir)/nlopt.i

nlopt-enum-renames.i: $(top_srcdir)/api/nlopt.h
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; egrep 'NLOPT_[LG][DN]|NLOPT_AUGLAG|NLOPT_G_MLSL|NLOPT_NUM_ALGORITHMS' $(top_srcdir)/api/nlopt.h | sed 's/NLOPT_//g' |tr -d ' ' |tr '/' ',' |tr '=' ',' |cut -d, -f1 |while read name; do echo "%rename(NLOPT_$$name) nlopt::$$name;"; done; egrep 'NLOPT_[A-Z_]* =' $(top_srcdir)/api/nlopt.h | egrep -v 'NLOPT_[LG][DN]|NLOPT_AUGLAG|NLOPT_G_MLSL' | sed 's/NLOPT_//g' |tr -d ' ' |cut -d'=' -f1 | while read name; do echo "%rename(NLOPT_$$name) nlopt::$$name;"; done) > $@

endif
