AM_CPPFLAGS = -I$(top_srcdir)/api -I$(top_srcdir)/util

MFILES = NLOPT_GN_DIRECT.m NLOPT_GN_DIRECT_L.m NLOPT_GN_DIRECT_L_RAND.m NLOPT_GN_DIRECT_NOSCAL.m NLOPT_GN_DIRECT_L_NOSCAL.m NLOPT_GN_DIRECT_L_RAND_NOSCAL.m NLOPT_GN_ORIG_DIRECT.m NLOPT_GN_ORIG_DIRECT_L.m NLOPT_GD_STOGO.m NLOPT_GD_STOGO_RAND.m NLOPT_LD_LBFGS_NOCEDAL.m NLOPT_LD_LBFGS.m NLOPT_LN_PRAXIS.m NLOPT_LD_VAR1.m NLOPT_LD_VAR2.m NLOPT_LD_TNEWTON.m NLOPT_LD_TNEWTON_RESTART.m NLOPT_LD_TNEWTON_PRECOND.m NLOPT_LD_TNEWTON_PRECOND_RESTART.m NLOPT_GN_CRS2_LM.m NLOPT_GN_MLSL.m NLOPT_GD_MLSL.m NLOPT_GN_MLSL_LDS.m NLOPT_GD_MLSL_LDS.m NLOPT_LD_MMA.m NLOPT_LN_COBYLA.m NLOPT_LN_NEWUOA.m NLOPT_LN_NEWUOA_BOUND.m NLOPT_LN_NELDERMEAD.m NLOPT_LN_SBPLX.m NLOPT_LN_AUGLAG.m NLOPT_LD_AUGLAG.m NLOPT_LN_AUGLAG_EQ.m NLOPT_LD_AUGLAG_EQ.m NLOPT_LN_BOBYQA.m NLOPT_GN_ISRES.m NLOPT_AUGLAG.m NLOPT_AUGLAG_EQ.m NLOPT_G_MLSL.m NLOPT_G_MLSL_LDS.m NLOPT_LD_SLSQP.m NLOPT_LD_CCSAQ.m NLOPT_GN_ESCH.m 

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
m_DATA = $(MFILES) nlopt_minimize.m nlopt_minimize_constrained.m nlopt_optimize.m
endif

nlopt_optimize.oct: nlopt_optimize-oct.cc nlopt_optimize_usage.h $(top_builddir)/.libs/libnlopt@NLOPT_SUFFIX@.la
	$(MKOCTFILE) -o $@ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(srcdir)/nlopt_optimize-oct.cc $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

BUILT_SOURCES = nlopt_optimize_usage.h
MAINTAINERCLEANFILES = $(BUILT_SOURCES)

nlopt_optimize_usage.h: $(srcdir)/nlopt_optimize.m
	echo "#define NLOPT_OPTIMIZE_USAGE \\" > $@
	sed 's/\"/\\"/g' $(srcdir)/nlopt_optimize.m | sed 's,^% ,\",;s,^%,\",;s,$$,\\n\" \\,' >> $@
	echo "" >> $@

#######################################################################
mexdir = $(MEX_INSTALL_DIR)
if WITH_MATLAB
mex_DATA = nlopt_optimize.$(MEXSUFF) $(MFILES) nlopt_minimize.m nlopt_minimize_constrained.m nlopt_optimize.m
endif

nlopt_optimize.$(MEXSUFF): nlopt_optimize-mex.c $(top_builddir)/.libs/libnlopt@NLOPT_SUFFIX@.la
	$(MEX) -output nlopt_optimize -O $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(srcdir)/nlopt_optimize-mex.c $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

#######################################################################

EXTRA_DIST = nlopt_optimize-oct.cc nlopt_optimize-mex.c $(MFILES) nlopt_minimize.m nlopt_minimize_constrained.m nlopt_optimize.m nlopt_optimize_usage.h

CLEANFILES = nlopt_optimize.oct nlopt_optimize.h nlopt_optimize.$(MEXSUFF) nlopt_optimize-oct.o
