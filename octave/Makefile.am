AM_CPPFLAGS = -I$(top_srcdir)/api -I$(top_srcdir)/util

MFILES = nlopt_minimize.m nlopt_minimize_constrained.m nlopt_optimize.m \
	nlopt_algorithm.m

#######################################################################
# dummy C program to fool automake into making sure that CPPFLAGS etc.
# are all defined
noinst_PROGRAMS = dummy
dummy_SOURCES = dummy.c

#######################################################################
octdir = $(OCT_INSTALL_DIR)
mdir = $(M_INSTALL_DIR)
if WITH_OCTAVE
oct_DATA = nlopt_optimize.oct
m_DATA = $(MFILES)
endif

nlopt_optimize.oct: nlopt_optimize-oct.cc nlopt_optimize_usage.h $(top_builddir)/.libs/libnlopt@NLOPT_SUFFIX@.la
	$(MKOCTFILE) -o $@ $(DEFS) \
		$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
		$(srcdir)/nlopt_optimize-oct.cc \
		$(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

BUILT_SOURCES = nlopt_optimize_usage.h
MAINTAINERCLEANFILES = $(BUILT_SOURCES)

nlopt_optimize_usage.h: $(srcdir)/nlopt_optimize.m
	echo "#define NLOPT_OPTIMIZE_USAGE \\" > $@
	cat $(srcdir)/nlopt_optimize.m \
		| sed 's/\"/\\"/g' \
		| sed 's,^% ,\",;s,^%,\",;s,$$,\\n\" \\,' >> $@
	echo "" >> $@

#######################################################################
mexdir = $(MEX_INSTALL_DIR)
if WITH_MATLAB
mex_DATA = nlopt_optimize.$(MEXSUFF) $(MFILES)
endif

nlopt_optimize.$(MEXSUFF): nlopt_optimize-mex.c $(top_builddir)/.libs/libnlopt@NLOPT_SUFFIX@.la
	$(MEX) -output nlopt_optimize -O $(DEFS) \
		$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
		$(srcdir)/nlopt_optimize-mex.c \
		$(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

#######################################################################

EXTRA_DIST = nlopt_optimize-oct.cc nlopt_optimize-mex.c \
	nlopt_optimize_usage.h $(MFILES)

CLEANFILES = nlopt_optimize.$(MEXSUFF) nlopt_optimize.oct \
	nlopt_optimize-oct.o nlopt_optimize_usage.h
