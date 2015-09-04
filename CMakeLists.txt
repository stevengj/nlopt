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
cmake_minimum_required (VERSION 2.6)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project (nlopt C)

option (WITH_CXX "enable cxx routines" OFF)
option (BUILD_SHARED_LIBS "Build NLOPT as a shared library" OFF)
option (BUILD_PYTHON "build python bindings" ON)

if (WITH_CXX)
  enable_language (CXX)
endif ()

INCLUDE(CheckIncludeFiles)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckTypeSize)
include (CheckCCompilerFlag)

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

SET ( NLOPT_HEADERS 
  api/nlopt.h api/nlopt.hpp api/nlopt.f
)

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
        )

if (WITH_CXX)
  list (APPEND NLOPT_SOURCES stogo/global.cc stogo/linalg.cc stogo/local.cc stogo/stogo.cc stogo/tools.cc stogo/global.h stogo/linalg.h stogo/local.h stogo/stogo_config.h stogo/stogo.h stogo/tools.h)
endif ()


if (BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DNLOPT_DLL)
  ADD_DEFINITIONS(-DNLOPT_DLL_EXPORT )
  
  if (NOT CMAKE_INSTALL_PREFIX MATCHES "^/usr")
    set (CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
    set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

INSTALL ( FILES ${NLOPT_HEADERS} DESTINATION include )

ADD_LIBRARY (nlopt ${NLOPT_SOURCES} )

# pass -fPIC in case swig module is built with static library
if (NOT BUILD_SHARED_LIBS)
  if (CMAKE_VERSION VERSION_LESS 2.8.9)
    check_c_compiler_flag (-fPIC HAS_FPIC)
    if (HAS_FPIC)
      set (CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
    endif ()
  else ()
    set_target_properties (${nlopt_library} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  endif ()
endif ()

# parse configure.ac for lib version
STRING( REGEX MATCH "SHARED_VERSION_INFO=\"[0-9]+:[0-9]+:[0-9]+\"" _VERSION_INFO_LINE "${CONFIG_AUTOMAKE}") 
IF (_VERSION_INFO_LINE)
  STRING( REGEX REPLACE "SHARED_VERSION_INFO=\"([0-9]+):[0-9]+:[0-9]+\"" "\\1" _VERSION_INFO_CURRENT "${_VERSION_INFO_LINE}") 
  STRING( REGEX REPLACE "SHARED_VERSION_INFO=\"[0-9]+:([0-9]+):[0-9]+\"" "\\1" _VERSION_INFO_REVISION "${_VERSION_INFO_LINE}") 
  STRING( REGEX REPLACE "SHARED_VERSION_INFO=\"[0-9]+:[0-9]+:([0-9]+)\"" "\\1" _VERSION_INFO_AGE "${_VERSION_INFO_LINE}") 
  MATH(EXPR _INFO_SOVERSION "${_VERSION_INFO_CURRENT} - ${_VERSION_INFO_AGE} ")
  SET(_INFO_VERSION "${_INFO_SOVERSION}.${_VERSION_INFO_AGE}.${_VERSION_INFO_REVISION}")
  SET_TARGET_PROPERTIES(nlopt PROPERTIES SOVERSION ${_INFO_SOVERSION})
  SET_TARGET_PROPERTIES(nlopt PROPERTIES VERSION ${_INFO_VERSION} )
ENDIF ()

INSTALL ( TARGETS nlopt
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib${LIB_SUFFIX}
          ARCHIVE DESTINATION lib${LIB_SUFFIX}
        )

if (BUILD_PYTHON)
  find_package (PythonInterp)
  find_package (PythonLibs)
endif ()

if (PYTHONINTERP_FOUND)
  execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
                    OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

  get_filename_component (_ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE)
  file (RELATIVE_PATH _REL_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH})

  set (PYTHON_MODULE_PATH ${_REL_PYTHON_MODULE_PATH})

endif ()

add_subdirectory (swig)

enable_testing ()
add_subdirectory (test)
