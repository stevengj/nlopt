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

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option (WITH_CXX "enable cxx routines" OFF)
option (BUILD_SHARED_LIBS "Build NLOPT as a shared library" OFF)
option (BUILD_PYTHON "build python bindings" ON)
option (BUILD_OCTAVE "build octave bindings" ON)
option (BUILD_GUILE "build guile bindings" ON)

set (NLOPT_SUFFIX)
if (WITH_CXX)
  set (NLOPT_SUFFIX _cxx)
endif ()

if (WITH_CXX OR BUILD_PYTHON)
  enable_language (CXX)
endif ()

include (CheckIncludeFiles)
include (CheckFunctionExists)
include (CheckTypeSize)
include (CheckCCompilerFlag)

#==============================================================================
# COMPILATION CHECKINGS and CONFIGURATION GENERATION
#==============================================================================
check_include_file (dlfcn.h HAVE_DLFCN_H)
check_include_file (getopt.h HAVE_GETOPT_H)
check_include_file (unistd.h HAVE_UNISTD_H)
check_include_file (string.h HAVE_STRING_H)
check_include_file (strings.h HAVE_STRINGS_H)
check_include_file (inttypes.h HAVE_INTTYPES_H)
check_include_file (memory.h HAVE_MEMORY_H)
check_include_file (stdlib.h HAVE_STDLIB_H)
check_include_file (stdint.h HAVE_STDINT_H)
check_include_file (time.h HAVE_TIME_H)
check_include_file (sys/types.h HAVE_SYS_TYPES_H)
check_include_file (sys/stat.h HAVE_SYS_STAT_H)
check_include_file (sys/time.h HAVE_SYS_TIME_H)
if (HAVE_TIME_H AND HAVE_SYS_TIME_H)
  set (TIME_WITH_SYS_TIME TRUE)
endif ()
check_function_exists (getpid HAVE_GETPID)
check_function_exists (syscall HAVE_GETTID_SYSCALL)
check_function_exists (isinf HAVE_ISINF)
check_function_exists (isnan HAVE_ISNAN)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (qsort_r HAVE_QSORT_R)
check_function_exists (time HAVE_TIME)
check_function_exists (copysign HAVE_COPYSIGN)
check_type_size ("uint32_t" SIZEOF_UINT32_T)
set (HAVE_UINT32_T ${SIZEOF_UINT32_T})
check_type_size ("unsigned int" SIZEOF_UNSIGNED_INT)
check_type_size ("unsigned long" SIZEOF_UNSIGNED_LONG)

if (NOT DEFINED HAVE_FPCLASSIFY)
  message(STATUS "Looking for fpclassify")
  file (WRITE ${CMAKE_BINARY_DIR}/fpclassify.c "#include <math.h>\n")
  file (APPEND ${CMAKE_BINARY_DIR}/fpclassify.c "int main(void) {\n")
  file (APPEND ${CMAKE_BINARY_DIR}/fpclassify.c "if (!fpclassify(3.14159)) fpclassify(2.7183);\n")
  file (APPEND ${CMAKE_BINARY_DIR}/fpclassify.c "  return 0; }\n")
  try_compile (HAVE_FPCLASSIFY
  ${CMAKE_BINARY_DIR}/build_fpclassify
  ${CMAKE_BINARY_DIR}/fpclassify.c
  CMAKE_FLAGS -DLINK_LIBRARIES=m)
  message(STATUS "Looking for fpclassify - ${HAVE_FPCLASSIFY}")
endif ()

option (WITH_THREADLOCAL "check thread local keyword" ON)
if (WITH_THREADLOCAL AND NOT DEFINED HAVE_THREAD_LOCAL_STORAGE)
  foreach (_THREADLOCAL_KEY "__thread" "__declspec(thread)")
    unset (HAVE_THREAD_LOCAL_STORAGE CACHE)
    check_c_source_compiles("
    ${_THREADLOCAL_KEY} int tls;

    int main(void) {
        return 0;
    }" HAVE_THREAD_LOCAL_STORAGE)
    if (${HAVE_THREAD_LOCAL_STORAGE})
      set (THREADLOCAL ${_THREADLOCAL_KEY})
    endif ()
  endforeach()
endif ()

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
    set (NLOPT_VERSION_STRING ${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION}.${NLOPT_BUGFIX_VERSION})
#---------------------------------
# Else (consider that version matches 'X.X')
#---------------------------------
ELSE(${NLOPT_AUTOMAKE_LINE_VERSION} MATCHES ".*[0-9]+\\.[0-9]+\\.[0-9]+")
    STRING( REGEX REPLACE ".*([0-9]+)\\.[0-9]+.*" "\\1" NLOPT_MAJOR_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
    STRING( REGEX REPLACE ".*[0-9]+\\.([0-9]+).*" "\\1" NLOPT_MINOR_VERSION "${NLOPT_AUTOMAKE_LINE_VERSION}") 
    SET (NLOPT_BUGFIX_VERSION "0")
    set (NLOPT_VERSION_STRING ${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION})
endif()

message (STATUS "NLOPT: Version number ${NLOPT_VERSION_STRING} found in configure.ac" )

#==============================================================================
# CREATE config.h
#==============================================================================

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h IMMEDIATE )

# pkgconfig file
set (prefix "${CMAKE_INSTALL_PREFIX}")
set (exec_prefix "\${prefix}")
set (libdir "\${exec_prefix}/lib")
set (includedir "\${prefix}/include")
set (VERSION "${NLOPT_VERSION_STRING}")
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt.pc.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc @ONLY)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc DESTINATION lib/pkgconfig)

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

set (nlopt_lib nlopt${NLOPT_SUFFIX})
add_library (${nlopt_lib} ${NLOPT_SOURCES})

# pass -fPIC in case swig module is built with static library
if (NOT BUILD_SHARED_LIBS)
  check_c_compiler_flag (-fPIC HAS_FPIC)
  if (HAS_FPIC)
    set (CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
  endif ()
endif ()

# parse configure.ac for lib version
string (REGEX MATCH "SHARED_VERSION_INFO=\"[0-9]+:[0-9]+:[0-9]+\"" _VERSION_INFO_LINE "${CONFIG_AUTOMAKE}") 
if (_VERSION_INFO_LINE)
  string (REGEX REPLACE "SHARED_VERSION_INFO=\"([0-9]+):[0-9]+:[0-9]+\"" "\\1" _VERSION_INFO_CURRENT "${_VERSION_INFO_LINE}") 
  string (REGEX REPLACE "SHARED_VERSION_INFO=\"[0-9]+:([0-9]+):[0-9]+\"" "\\1" _VERSION_INFO_REVISION "${_VERSION_INFO_LINE}") 
  string (REGEX REPLACE "SHARED_VERSION_INFO=\"[0-9]+:[0-9]+:([0-9]+)\"" "\\1" _VERSION_INFO_AGE "${_VERSION_INFO_LINE}") 
  math (EXPR _INFO_SOVERSION "${_VERSION_INFO_CURRENT} - ${_VERSION_INFO_AGE} ")
  set (_INFO_VERSION "${_INFO_SOVERSION}.${_VERSION_INFO_AGE}.${_VERSION_INFO_REVISION}")
  set_target_properties (${nlopt_lib} PROPERTIES SOVERSION ${_INFO_SOVERSION})
  set_target_properties (${nlopt_lib} PROPERTIES VERSION ${_INFO_VERSION} )
endif ()

install ( TARGETS ${nlopt_lib}
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

if (BUILD_GUILE)
  find_package (Guile)
endif ()

add_subdirectory (swig)

if (BUILD_OCTAVE)
  find_package (Octave)
endif ()

if (OCTAVE_FOUND)
  add_subdirectory (octave)
endif ()

enable_testing ()
add_subdirectory (test)

set (CPACK_PACKAGE_NAME          "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_VERSION_MAJOR "${NLOPT_MAJOR_VERSION}")
set (CPACK_PACKAGE_VERSION_MINOR "${NLOPT_MINOR_VERSION}")
set (CPACK_PACKAGE_VERSION_PATCH "${NLOPT_BUGFIX_VERSION}")
set (CPACK_SOURCE_GENERATOR      "TBZ2;TGZ"    )
set (CPACK_BINARY_STGZ           "OFF"        )
set (CPACK_BINARY_TBZ2           "ON"         )
set (CPACK_BINARY_TGZ            "ON"         )
set (CPACK_BINARY_TZ             "OFF"        )
set (CPACK_SOURCE_IGNORE_FILES ".git;/build;.*~;${CPACK_SOURCE_IGNORE_FILES}")

set (CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${NLOPT_VERSION_STRING})

include (CPack)