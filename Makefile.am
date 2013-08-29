OPTIONS_AUTOMAKE=gnu
lib_LTLIBRARIES = libnlopt@NLOPT_SUFFIX@.la

ACLOCAL_AMFLAGS=-I m4

if WITH_CXX
CXX_DIRS = stogo
CXX_LIBS = stogo/libstogo.la
endif

SUBDIRS = util direct cdirect $(CXX_DIRS) praxis luksan crs mlsl mma cobyla newuoa neldermead auglag bobyqa isres slsqp esch api . octave test swig
EXTRA_DIST = autogen.sh nlopt.pc.in m4

libnlopt@NLOPT_SUFFIX@_la_SOURCES = 
libnlopt@NLOPT_SUFFIX@_la_LIBADD = direct/libdirect.la			\
cdirect/libcdirect.la $(CXX_LIBS) praxis/libpraxis.la \
luksan/libluksan.la crs/libcrs.la mlsl/libmlsl.la mma/libmma.la		\
cobyla/libcobyla.la newuoa/libnewuoa.la neldermead/libneldermead.la	\
auglag/libauglag.la bobyqa/libbobyqa.la isres/libisres.la		\
slsqp/libslsqp.la esch/libesch.la api/libapi.la util/libutil.la

if WITH_CXX
libnlopt@NLOPT_SUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@
else
libnlopt@NLOPT_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = nlopt.pc

if MAINTAINER_MODE

README: README.md
	cp -f $(srcdir)/README.md $@

endif
