OPTIONS_AUTOMAKE=gnu
lib_LTLIBRARIES = libnlopt@NLOPT_SUFFIX@.la

ACLOCAL_AMFLAGS=-I m4

if WITH_CXX
CXX_DIRS = stogo
CXX_LIBS = stogo/libstogo.la
endif

SUBDIRS = util direct cdirect $(CXX_DIRS) praxis luksan crs mlsl mma cobyla newuoa neldermead auglag bobyqa isres slsqp esch api . octave test swig
EXTRA_DIST = autogen.sh nlopt.pc.in m4 CMakeLists.txt config.cmake.h.in cmake/FindOctave.cmake cmake/FindGuile.cmake 

libnlopt@NLOPT_SUFFIX@_la_SOURCES =
if WITH_CXX
# Dummy C++ source to cause C++ linking
nodist_EXTRA_libnlopt@NLOPT_SUFFIX@_la_SOURCES = dummy.cxx
endif
libnlopt@NLOPT_SUFFIX@_la_LIBADD = direct/libdirect.la			\
cdirect/libcdirect.la $(CXX_LIBS) praxis/libpraxis.la \
luksan/libluksan.la crs/libcrs.la mlsl/libmlsl.la mma/libmma.la		\
cobyla/libcobyla.la newuoa/libnewuoa.la neldermead/libneldermead.la	\
auglag/libauglag.la bobyqa/libbobyqa.la isres/libisres.la		\
slsqp/libslsqp.la esch/libesch.la api/libapi.la util/libutil.la

libnlopt@NLOPT_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = nlopt.pc

if MAINTAINER_MODE

README: README.md
	cp -f $(srcdir)/README.md $@

endif
