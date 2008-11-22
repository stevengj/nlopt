AM_CPPFLAGS = -I$(top_srcdir)/api -I$(top_srcdir)/util

MFILES = NLOPT_GN_DIRECT.m NLOPT_GN_DIRECT_L.m NLOPT_GN_DIRECT_L_RAND.m NLOPT_GN_DIRECT_NOSCAL.m NLOPT_GN_DIRECT_L_NOSCAL.m NLOPT_GN_DIRECT_L_RAND_NOSCAL.m NLOPT_GN_ORIG_DIRECT.m NLOPT_GN_ORIG_DIRECT_L.m NLOPT_GD_STOGO.m NLOPT_GD_STOGO_RAND.m NLOPT_LD_LBFGS_NOCEDAL.m NLOPT_LD_LBFGS.m NLOPT_LN_PRAXIS.m NLOPT_LD_VAR1.m NLOPT_LD_VAR2.m NLOPT_LD_TNEWTON.m NLOPT_LD_TNEWTON_RESTART.m NLOPT_LD_TNEWTON_PRECOND.m NLOPT_LD_TNEWTON_PRECOND_RESTART.m NLOPT_GN_CRS2_LM.m NLOPT_GN_MLSL.m NLOPT_GD_MLSL.m NLOPT_GN_MLSL_LDS.m NLOPT_GD_MLSL_LDS.m NLOPT_LD_MMA.m NLOPT_LN_COBYLA.m NLOPT_LN_NEWUOA.m NLOPT_LN_NEWUOA_BOUND.m NLOPT_LN_NELDERMEAD.m NLOPT_LN_SBPLX.m 

#######################################################################
# dummy C program to fool automake into making sure that CPPFLAGS etc.
# are all defined
noinst_PROGRAMS = dummy
dummy_SOURCES = dummy.c

#######################################################################
octdir = $(OCT_INSTALL_DIR)
mdir = $(M_INSTALL_DIR)
if WITH_OCTAVE
oct_DATA = nlopt_minimize_constrained.oct
m_DATA = $(MFILES) nlopt_minimize.m nlopt_minimize_constrained.m
endif

nlopt_minimize_constrained.oct: nlopt_minimize_constrained-oct.cc nlopt_minimize_constrained_usage.h $(top_builddir)/.libs/libnlopt@NLOPT_SUFFIX@.la
	$(MKOCTFILE) -o $@ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(srcdir)/nlopt_minimize_constrained-oct.cc $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

nlopt_minimize_constrained_usage.h: $(srcdir)/nlopt_minimize_constrained.m
	echo "#define NLOPT_MINIMIZE_CONSTRAINED_USAGE \\" > $@
	sed 's/\"/\\"/g' $(srcdir)/nlopt_minimize_constrained.m | sed 's,^% ,\",;s,^%,\",;s,$$,\\n\" \\,' >> $@
	echo "" >> $@

#######################################################################
mexdir = $(MEX_INSTALL_DIR)
if WITH_MATLAB
mex_DATA = nlopt_minimize_constrained.$(MEXSUFF) $(MFILES) nlopt_minimize.m nlopt_minimize_constrained.m
endif

nlopt_minimize_constrained.$(MEXSUFF): nlopt_minimize_constrained-mex.c $(top_builddir)/.libs/libnlopt@NLOPT_SUFFIX@.la
	$(MEX) -output nlopt_minimize_constrained -O $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(srcdir)/nlopt_minimize_constrained-mex.c $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

#######################################################################

EXTRA_DIST = nlopt_minimize_constrained-oct.cc nlopt_minimize_constrained-mex.c $(MFILES) nlopt_minimize.m nlopt_minimize_constrained.m

CLEANFILES = nlopt_minimize_constrained.oct nlopt_minimize_constrained_usage.h nlopt_minimize_constrained.$(MEXSUFF) nlopt_minimize_constrained-oct.o
