OPTIONS_AUTOMAKE=gnu
lib_LTLIBRARIES = libnlopt.la

ACLOCAL_AMFLAGS=-I ./m4

SUBDIRS=lbfgs subplex direct stogo api
EXTRA_DIST=COPYRIGHT autogen.sh nlopt.pc.in m4

libnlopt_la_SOURCES = 
libnlopt_la_LIBADD = lbfgs/liblbfgs.la subplex/libsubplex.la direct/libdirect.la stogo/libstogo.la api/libapi.la

libnlopt_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = nlopt.pc

# Somewhat hackish.  The "right" way to do this is by a dist-hook target,
# but then darcs check will fail because it doesn't run in the darcs
# repository.
darcs-dist: distdir
	darcs changes --summary > $(distdir)/ChangeLog
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	$(am__remove_distdir)
