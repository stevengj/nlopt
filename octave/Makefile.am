AM_CPPFLAGS = -I$(top_srcdir)/api 

MFILES = NLOPT_GN_DIRECT.m NLOPT_GN_DIRECT_L.m NLOPT_GN_DIRECT_L_RAND.m NLOPT_GN_DIRECT_NOSCAL.m NLOPT_GN_DIRECT_L_NOSCAL.m NLOPT_GN_DIRECT_L_RAND_NOSCAL.m NLOPT_GN_ORIG_DIRECT.m NLOPT_GN_ORIG_DIRECT_L.m NLOPT_LN_SUBPLEX.m NLOPT_GD_STOGO.m NLOPT_GD_STOGO_RAND.m NLOPT_LD_LBFGS_NOCEDAL.m NLOPT_LD_LBFGS.m NLOPT_LN_PRAXIS.m NLOPT_LD_VAR1.m NLOPT_LD_VAR2.m NLOPT_LD_TNEWTON.m NLOPT_LD_TNEWTON_RESTART.m NLOPT_LD_TNEWTON_PRECOND.m NLOPT_LD_TNEWTON_PRECOND_RESTART.m NLOPT_GN_CRS2_LM.m NLOPT_GN_MLSL.m NLOPT_GD_MLSL.m NLOPT_GN_MLSL_LDS.m NLOPT_GD_MLSL_LDS.m NLOPT_LD_MMA.m NLOPT_LN_COBYLA.m 

#######################################################################
octdir = $(OCT_INSTALL_DIR)
mdir = $(M_INSTALL_DIR)
if WITH_OCTAVE
oct_DATA = nlopt_minimize.oct
m_DATA = $(MFILES) nlopt_minimize.m
endif

nlopt_minimize.oct: nlopt_minimize-oct.cc nlopt_minimize_usage.h
	$(MKOCTFILE) -o $@ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(srcdir)/nlopt_minimize-oct.cc $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

nlopt_minimize_usage.h: $(srcdir)/nlopt_minimize.m
	echo "#define NLOPT_MINIMIZE_USAGE \\" > $@
	sed 's/\"/\\"/g' $(srcdir)/nlopt_minimize.m | sed 's,^% ,\",;s,^%,\",;s,$$,\\n\" \\,' >> $@
	echo "" >> $@

#######################################################################
mexdir = $(MEX_INSTALL_DIR)
if WITH_MATLAB
mex_DATA = nlopt_minimize.$(MEXSUFF) $(MFILES) nlopt_minimize.m
endif

nlopt_minimize.$(MEXSUFF): nlopt_minimize-mex.c
	$(MEX) -output nlopt_minimize -O $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(srcdir)/nlopt_minimize-mex.c $(LDFLAGS) -L$(top_builddir)/.libs -lnlopt@NLOPT_SUFFIX@

#######################################################################

EXTRA_DIST = nlopt_minimize-oct.cc nlopt_minimize-mex.c $(MFILES) nlopt_minimize.m

CLEANFILES = nlopt_minimize.oct nlopt_minimize_usage.h nlopt_minimize.$(MEXSUFF) nlopt_minimize-oct.o
