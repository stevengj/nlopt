OPTIONS_AUTOMAKE=gnu
lib_LTLIBRARIES = libnlopt@NLOPT_SUFFIX@.la

ACLOCAL_AMFLAGS=-I m4

if WITH_CXX
CXX_DIRS = stogo
CXX_LIBS = stogo/libstogo.la
endif

SUBDIRS = util direct cdirect $(CXX_DIRS) praxis luksan crs mlsl mma cobyla newuoa lbfgs neldermead auglag bobyqa isres api . octave test
EXTRA_DIST = autogen.sh nlopt.pc.in m4

if WITH_NOCEDAL
NOCEDAL_LBFGS=lbfgs/liblbfgs.la
endif

libnlopt@NLOPT_SUFFIX@_la_SOURCES = 
libnlopt@NLOPT_SUFFIX@_la_LIBADD = direct/libdirect.la			\
cdirect/libcdirect.la $(CXX_LIBS) praxis/libpraxis.la $(NOCEDAL_LBFGS)	\
luksan/libluksan.la crs/libcrs.la mlsl/libmlsl.la mma/libmma.la		\
cobyla/libcobyla.la newuoa/libnewuoa.la neldermead/libneldermead.la	\
auglag/libauglag.la bobyqa/libbobyqa.la isres/libisres.la		\
api/libapi.la util/libutil.la

libnlopt@NLOPT_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = nlopt.pc

# Somewhat hackish.  The "right" way to do this is by a dist-hook target,
# but then darcs check will fail because it doesn't run in the darcs
# repository.
darcs-dist: distdir
	darcs changes --summary > $(distdir)/ChangeLog
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	$(am__remove_distdir)
