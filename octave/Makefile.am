AM_CPPFLAGS = -I$(top_srcdir)/api 

MFILES = NLOPT_GLOBAL_DIRECT_L.m NLOPT_GLOBAL_DIRECT_L_RANDOMIZED.m	\
NLOPT_GLOBAL_DIRECT.m NLOPT_GLOBAL_ORIG_DIRECT_L.m			\
NLOPT_GLOBAL_ORIG_DIRECT.m NLOPT_GLOBAL_STOGO.m				\
NLOPT_GLOBAL_STOGO_RANDOMIZED.m NLOPT_LOCAL_LBFGS.m			\
NLOPT_LOCAL_SUBPLEX.m nlopt_minimize.m

octdir = $(OCT_INSTALL_DIR)
mdir = $(M_INSTALL_DIR)

if WITH_OCTAVE
oct_DATA = nlopt_minimize.oct
m_DATA = $(MFILES)
endif

nlopt_minimize.oct: nlopt_minimize.cc nlopt_minimize_usage.h
	$(MKOCTFILE) $(DEFS) $(CPPFLAGS) $(srcdir)/nlopt_minimize.cc $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt

nlopt_minimize_usage.h: $(srcdir)/nlopt_minimize.m
	echo "#define NLOPT_MINIMIZE_USAGE \\" > $@
	sed 's/\"/\\"/g' $(srcdir)/nlopt_minimize.m | sed 's,^% ,\",;s,^%,\",;s,$$,\\n\" \\,' >> $@
	echo "" >> $@

EXTRA_DIST = nlopt_minimize.cc $(MFILES)

CLEANFILES = nlopt_minimize.oct nlopt_minimize_usage.h
