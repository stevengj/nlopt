PROJECT(NLOPT)

#==============================================================================
# NLOPT CMake file
# 
# NLopt is a free/open-source library for nonlinear optimization, providing 
# a common interface for a number of different free optimization routines 
# available online as well as original implementations of various other 
# algorithms
# WEBSITE: http://ab-initio.mit.edu/wiki/index.php/NLopt 
# AUTHOR: Steven G. Johnson
#
# This CMakeLists.txt file was created to compile NLOPT with the CMAKE utility.
# Benoit Scherrer, 2010 CRL, Harvard Medical School
# Copyright (c) 2008-2009 Children's Hospital Boston 
#==============================================================================


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckTypeSize)


#==============================================================================
# COMPILATION CHECKINGS and CONFIGURATION GENERATION
#==============================================================================
CHECK_INCLUDE_FILES( unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES( string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES( strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES( sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES( inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES( memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILES( stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES( stdint.h HAVE_STDINT_H)

CHECK_INCLUDE_FILES( sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES( sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES( sys/types.h HAVE_SYS_TYPES_H)

CHECK_FUNCTION_EXISTS( isinf HAVE_ISINF)
CHECK_FUNCTION_EXISTS( isinf HAVE_ISNAN)
CHECK_FUNCTION_EXISTS( gettimeofday HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS( qsort_r HAVE_QSORT_R)
CHECK_FUNCTION_EXISTS( time HAVE_TIME)
CHECK_FUNCTION_EXISTS( copysign HAVE_COPYSIGN)

CHECK_TYPE_SIZE(uint32_t    UINT32_T)
IF(NOT HAVE_UINT32_T)
  IF(MSVC)
    SET(uint32_t "unsigned int")
    SET(SIZEOF_UNSIGNED_INT 4)
  ENDIF(MSVC)
ENDIF(NOT HAVE_UINT32_T)

#==============================================================================
# LOOK FOR THE VERSION NUMBER IN configure.ac
#==============================================================================

FILE( READ configure.ac CONFIG_AUTOMAKE )

STRING( REGEX MATCH "AC_INIT\\(nlopt, ([0-9]+\\.[0-9]+\\.?[0-9]*)" NLOPT_AUTOMAKE_LINE_VERSION "${CONFIG_AUTOMAKE}") 

#---------------------------------
# If version matches 'X.X.X'
#---------------------------------
IF(${NLOPT_AUTOMAKE_LINE_VERSION} MATCHES ".*[0-9]+\\.[0-9]+\\.[0-9]+")
    STRING( REGEX REPLACE ".*([0-9]+)\\.[0-9]+\\.[0-9]*.*" "\\1" NLOPT_MAJOR_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
    STRING( REGEX REPLACE ".*[0-9]+\\.([0-9]+)\\.[0-9]*.*" "\\1" NLOPT_MINOR_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
    STRING( REGEX REPLACE ".*[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" NLOPT_BUGFIX_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
#---------------------------------
# Else (consider that version matches 'X.X')
#---------------------------------
ELSE(${NLOPT_AUTOMAKE_LINE_VERSION} MATCHES ".*[0-9]+\\.[0-9]+\\.[0-9]+")
    STRING( REGEX REPLACE ".*([0-9]+)\\.[0-9]+.*" "\\1" NLOPT_MAJOR_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
    STRING( REGEX REPLACE ".*[0-9]+\\.([0-9]+).*" "\\1" NLOPT_MINOR_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
    SET (NLOPT_BUGFIX_VERSION "0")
ENDIF(${NLOPT_AUTOMAKE_LINE_VERSION} MATCHES ".*[0-9]+\\.[0-9]+\\.[0-9]+")

MESSAGE(STATUS "NLOPT: Version number ${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION}.${NLOPT_BUGFIX_VERSION} found in configure.ac" )


#==============================================================================
# CREATE config.h
#==============================================================================

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h IMMEDIATE )


#==============================================================================
# INCLUDE DIRECTORIES
#==============================================================================
SET(${INCLUDE_DIRECTORIES} "")
INCLUDE_DIRECTORIES (
	${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
	stogo 
	util 
	direct 
	cdirect 
	praxis 
	luksan 
	crs 
	mlsl 
	mma 
	cobyla 
	newuoa 
	neldermead 
	auglag 
	bobyqa 
	isres 
	slsqp
        esch
	api   )


#==============================================================================
# nlopt LIBRARY TARGET (SHARED OR STATIC)
#==============================================================================

SET ( NLOPT_SOURCES   
	direct/DIRect.c direct/direct_wrap.c direct/DIRserial.c direct/DIRsubrout.c direct/direct-internal.h direct/direct.h
	cdirect/cdirect.c cdirect/hybrid.c cdirect/cdirect.h
	praxis/praxis.c praxis/praxis.h
	luksan/plis.c luksan/plip.c luksan/pnet.c luksan/mssubs.c luksan/pssubs.c luksan/luksan.h
	crs/crs.c crs/crs.h
	mlsl/mlsl.c mlsl/mlsl.h
	mma/mma.c mma/mma.h mma/ccsa_quadratic.c
	cobyla/cobyla.c cobyla/cobyla.h
	newuoa/newuoa.c newuoa/newuoa.h 
	neldermead/nldrmd.c neldermead/neldermead.h neldermead/sbplx.c	 
	auglag/auglag.c auglag/auglag.h
	bobyqa/bobyqa.c bobyqa/bobyqa.h
	isres/isres.c isres/isres.h 
	slsqp/slsqp.c slsqp/slsqp.h 
	esch/esch.c esch/esch.h
	api/general.c api/options.c api/optimize.c api/deprecated.c api/nlopt-internal.h api/nlopt.h api/f77api.c api/f77funcs.h api/f77funcs_.h api/nlopt.hpp api/nlopt-in.hpp
	util/mt19937ar.c util/sobolseq.c util/soboldata.h util/timer.c util/stop.c util/nlopt-util.h util/redblack.c util/redblack.h util/qsort_r.c util/rescale.c 
	stogo/global.cc stogo/linalg.cc stogo/local.cc stogo/stogo.cc stogo/tools.cc stogo/global.h stogo/linalg.h stogo/local.h stogo/stogo_config.h stogo/stogo.h stogo/tools.h 
        )

OPTION(NLOPT_BUILD_SHARED "Build NLOPT as a shared library" OFF )

IF(NLOPT_BUILD_SHARED)
  ADD_DEFINITIONS(-DNLOPT_DLL)
  ADD_DEFINITIONS(-DNLOPT_DLL_EXPORT )
  ADD_LIBRARY (nlopt SHARED ${NLOPT_SOURCES} )
ELSE(NLOPT_BUILD_SHARED)
  ADD_LIBRARY (nlopt STATIC ${NLOPT_SOURCES} )
ENDIF(NLOPT_BUILD_SHARED)



