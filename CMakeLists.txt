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

cmake_minimum_required (VERSION 2.8.5)

# default install directory locally inside build-tree
# (this must come before calling PROJECT)
if (NOT DEFINED CMAKE_INSTALL_PREFIX
    OR CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install"
    CACHE PATH "Installation directory prefix" FORCE)
endif ()

# default build type to optimized release
if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Build type: Debug|Release|RelWithDebInfo|MinSizeRel")
endif ()

# C language project
project (nlopt C)

# append path to custom cmake modules FindXXX
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# provide user options
option (BUILD_SHARED_LIBS   "Build NLopt as a shared library" ON)
option (WITH_CXX            "Enable C++-based routines (StoGO)" OFF)
option (WITH_SWIG           "Use SWIG to generate bindings" ON)
option (BUILD_PYTHON        "Build Python bindings (using SWIG)" ON)
option (BUILD_GUILE         "Build Guile bindings (using SWIG)" ON)
option (BUILD_OCTAVE        "Build Octave bindings" ON)
option (BUILD_MATLAB        "Build MATLAB bindings" ON)
option (BUILD_TESTS         "Build Tests" ON)
option (WITH_EXTRA_WARNINGS "Enable extra compiler warnings (GNU GCC)" OFF)

# validate options combination
if (BUILD_PYTHON OR BUILD_GUILE)
  if (NOT WITH_SWIG)
    message (WARNING "Python/Guile wrappers require SWIG")
  endif ()
endif ()

# enable C++ language compiler if needed
if (WITH_CXX OR BUILD_PYTHON OR BUILD_GUILE OR BUILD_OCTAVE)
  enable_language (CXX)
endif ()

set (NLOPT_SUFFIX)
if (WITH_CXX)
  set (NLOPT_SUFFIX _cxx)
endif ()

include (GNUInstallDirs)

# Offer the user the choice of overriding the installation directories
# (these are relative to CMAKE_INSTALL_PREFIX)
set (INSTALL_LIB_DIR     "${CMAKE_INSTALL_LIBDIR}"
  CACHE PATH "Installation directory for libraries")
set (INSTALL_BIN_DIR     "${CMAKE_INSTALL_BINDIR}"
  CACHE PATH "Installation directory for executables")
set (INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}"
  CACHE PATH "Installation directory for header files")
set (INSTALL_DATA_DIR    "${CMAKE_INSTALL_DATADIR}/nlopt"
  CACHE PATH "Installation directory for data files")
set (INSTALL_MAN_DIR     "${CMAKE_INSTALL_MANDIR}"
  CACHE PATH "Installation directory for man pages")
set (INSTALL_DOC_DIR       "${CMAKE_INSTALL_DOCDIR}"
  CACHE PATH "Installation directory for doc files")
set (INSTALL_CMAKE_DIR   "${INSTALL_LIB_DIR}/cmake/nlopt"
  CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA MAN DOC CMAKE)
  set (var INSTALL_${p}_DIR)
  set (RELATIVE_INSTALL_${p}_DIR "${INSTALL_${p}_DIR}")
  if (NOT IS_ABSOLUTE "${${var}}")
    set (${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif ()
endforeach ()

# rpath
set (CMAKE_INSTALL_RPATH "${INSTALL_LIB_DIR}")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if (POLICY CMP0042)
  # Set MACOSX_RPATH to ON
  cmake_policy (SET CMP0042 NEW)
endif ()

# set compiler-specific flags
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
  if (WITH_EXTRA_WARNINGS)
    # add extra warnings in maintainer mode
    foreach (_FLG -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align
        -Wno-long-long -pedantic -Wshadow -Wwrite-strings -Wredundant-decls
        -Wundef -Wconversion)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_FLG}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_FLG}")
    endforeach ()
    foreach (_FLG -Wbad-function-cast -Wstrict-prototypes -Wnested-externs
        -Wmissing-prototypes -Wmissing-declarations)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_FLG}")
    endforeach ()
    unset (_FLG)
  else ()
    # turn off deprecated warnings
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-declarations")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  endif ()
endif ()

#==============================================================================
# VERSION NUMBERS
#==============================================================================

# release version MAJOR.MINOR.BUGFIX
set (NLOPT_MAJOR_VERSION  2)
set (NLOPT_MINOR_VERSION  5)
set (NLOPT_BUGFIX_VERSION 0)
set (NLOPT_VERSION_STRING
  "${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION}.${NLOPT_BUGFIX_VERSION}")

# interface version CURRENT:REVISION:AGE (libtool's versioning)
set (NLOPT_SOVERSION_CURRENT  9)
set (NLOPT_SOVERSION_REVISION 0)
set (NLOPT_SOVERSION_AGE      9)
math (EXPR NLOPT_SOVERSION_INTERFACE
  "${NLOPT_SOVERSION_CURRENT} - ${NLOPT_SOVERSION_AGE}")
set (NLOPT_SOVERSION_STRING
  "${NLOPT_SOVERSION_INTERFACE}.${NLOPT_SOVERSION_AGE}.${NLOPT_SOVERSION_REVISION}")

#==============================================================================
# COMPILATION CHECKINGS
#==============================================================================

include (CheckIncludeFile)
include (CheckIncludeFiles)
include (CheckSymbolExists)
include (CheckTypeSize)
include (CheckCSourceCompiles)
include (CMakePushCheckState)

# AC_HEADER_STDC
check_include_files ("stdlib.h;stdarg.h;string.h;float.h" STDC_HEADERS)

# AC_HEADER_TIME
check_include_files ("sys/time.h;time.h" TIME_WITH_SYS_TIME)
check_include_file (sys/time.h HAVE_SYS_TIME_H)

# AC_CHECK_HEADERS for some POSIX headers
check_include_file (unistd.h HAVE_UNISTD_H)
check_include_file (getopt.h HAVE_GETOPT_H)
check_include_file (stdint.h HAVE_STDINT_H)

# AC_CHECK_LIB (libm mostly needed on UNIX)
cmake_push_check_state ()
set (CMAKE_REQUIRED_LIBRARIES m)
check_symbol_exists (sqrt math.h HAVE_LIBM)
cmake_pop_check_state ()

# AC_CHECK_FUNCS from math.h
cmake_push_check_state ()
if (UNIX OR HAVE_LIBM)
  set (CMAKE_REQUIRED_LIBRARIES m)
endif ()
check_symbol_exists (isnan math.h HAVE_ISNAN)
check_symbol_exists (isinf math.h HAVE_ISINF)
check_symbol_exists (isfinite math.h HAVE_ISFINITE)
check_symbol_exists (fpclassify math.h HAVE_FPCLASSIFY)
check_symbol_exists (copysign math.h HAVE_COPYSIGN)
cmake_pop_check_state ()

# AC_CHECK_FUNCS from sys/time.h;time.h
check_symbol_exists (gettimeofday sys/time.h HAVE_GETTIMEOFDAY)
check_symbol_exists (time time.h HAVE_TIME)

# AC_CHECK_FUNCS from unistd.h
check_symbol_exists (syscall unistd.h HAVE_GETTID_SYSCALL)
check_symbol_exists (getpid unistd.h HAVE_GETPID)

# AC_CHECK_FUNCS for reentrant quicksort version
check_symbol_exists (qsort_r stdlib.h HAVE_QSORT_R)

# AC_CHECK_FUNCS for formatted IO functions
check_symbol_exists (vsnprintf stdio.h HAVE_VSNPRINTF)

# check for 32-bit unsigned integer type uint32_t
# (NOTE: it autmatically checks and includes necessary standard headers)
check_type_size ("uint32_t" SIZEOF_UINT32_T)
set (HAVE_UINT32_T ${SIZEOF_UINT32_T})
check_type_size ("unsigned int" SIZEOF_UNSIGNED_INT)
check_type_size ("unsigned long" SIZEOF_UNSIGNED_LONG)

# check TLS keyword
if (NOT DEFINED HAVE_THREAD_LOCAL_STORAGE)
  foreach (_TLS_KEY "__thread" "__declspec(thread)")
    unset (HAVE_THREAD_LOCAL_STORAGE CACHE)
    check_c_source_compiles (
      "${_TLS_KEY} int tls = 0; int main(void) { tls = 1; return 0; }"
      HAVE_THREAD_LOCAL_STORAGE)
    if (HAVE_THREAD_LOCAL_STORAGE)
      set (THREADLOCAL ${_TLS_KEY} CACHE INTERNAL "Keyword TLS")
      break ()
    endif ()
  endforeach ()
  unset (_TLS_KEY)
endif ()

# check inline for C
if (NOT DEFINED HAVE_C_INLINE)
  foreach (_INLINE_KEY "inline" "__inline__" "__inline")
    unset (HAVE_C_INLINE CACHE)
    check_c_source_compiles ("
        typedef int foo_t;
        static ${_INLINE_KEY} foo_t foo() { return 0; }
        int main(void) { return foo(); }
      " HAVE_C_INLINE)
    if (HAVE_C_INLINE)
      set (C_INLINE ${_INLINE_KEY} CACHE INTERNAL "Keyword inline")
      break ()
    endif ()
  endforeach ()
  unset (_INLINE_KEY)
endif ()

# check for broken Solaris HUGE_VAL macro under gcc 3.4.x and similar
if (NOT DEFINED HAVE_HUGE_VAL)
  check_c_source_compiles (
    "#include <math.h>\nint main(void) { double x = -HUGE_VAL; return 0; }"
    HAVE_HUGE_VAL)
  if (NOT HAVE_HUGE_VAL)
    unset (HAVE_HUGE_VAL CACHE)
    check_c_source_compiles ("
        #include <math.h>
        #ifdef __GNUC__
        #undef HUGE_VAL
        #define HUGE_VAL __builtin_huge_val()
        #endif
        int main(void) { double x = -HUGE_VAL; return 0; }
      " HAVE_HUGE_VAL)
    if (HAVE_HUGE_VAL)
      set (REPLACEMENT_HUGE_VAL "__builtin_huge_val()" CACHE INTERNAL
        "Replacement for HUGE_VAL")
    endif ()
  endif ()
endif ()

#==============================================================================
# CONFIGURATION GENERATION
#==============================================================================

# create config.h file
configure_file (config.cmake.h.in config.h IMMEDIATE)

#==============================================================================
# NLOPT LIBRARY TARGET
#==============================================================================

set (${INCLUDE_DIRECTORIES} "")
include_directories (
  ${CMAKE_BINARY_DIR}/api
  ${CMAKE_BINARY_DIR}
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
  esch
  slsqp
  #cquad
  #subplex
  #tensor
  api
)

# source files: C
set (NLOPT_SOURCES
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
  esch/esch.c esch/esch.h
  slsqp/slsqp.c slsqp/slsqp.h
  #cquad/cquad.c cquad/cquad.h
  #subplex/subplex.c subplex/subplex.h
  #tensor/tensor.c
  api/general.c api/options.c api/optimize.c api/deprecated.c api/nlopt-internal.h api/nlopt.h api/f77api.c api/f77funcs.h api/f77funcs_.h api/nlopt.hpp api/nlopt-in.hpp
  util/mt19937ar.c util/sobolseq.c util/soboldata.h util/timer.c util/stop.c util/nlopt-util.h util/redblack.c util/redblack.h util/qsort_r.c util/rescale.c
)

# source files: C++
if (WITH_CXX)
  list (APPEND NLOPT_SOURCES
    stogo/global.cc stogo/linalg.cc stogo/local.cc stogo/stogo.cc stogo/tools.cc stogo/global.h stogo/linalg.h stogo/local.h stogo/stogo_config.h stogo/stogo.h stogo/tools.h)
endif ()

# dllimport/dllexport decorations for Windows
if (BUILD_SHARED_LIBS AND WIN32)
  add_definitions (-DNLOPT_DLL -DNLOPT_DLL_EXPORT)
endif ()

# library target (Shared or Static)
set (nlopt_lib nlopt${NLOPT_SUFFIX})
add_library (${nlopt_lib} ${NLOPT_SOURCES})

# link against libm if needed
if (UNIX OR HAVE_LIBM)
  target_link_libraries (${nlopt_lib} m)
endif ()

# library version
set_target_properties (${nlopt_lib} PROPERTIES
  VERSION "${NLOPT_SOVERSION_STRING}"
  SOVERSION ${NLOPT_SOVERSION_INTERFACE})

# pass -fPIC in case swig module is built with static library
if (NOT BUILD_SHARED_LIBS)
  check_c_compiler_flag (-fPIC HAS_FPIC)
  if (HAS_FPIC)
    set (CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
  endif ()
endif ()

# library installation
install (
  TARGETS ${nlopt_lib}
  EXPORT  NLoptLibraryDepends
  RUNTIME DESTINATION "${RELATIVE_INSTALL_BIN_DIR}"
  LIBRARY DESTINATION "${RELATIVE_INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${RELATIVE_INSTALL_LIB_DIR}")

#==============================================================================
# C/C++/FORTRAN HEADERS
#==============================================================================

# generate nlopt.hpp and nlopt.f
add_subdirectory (api)

# header files to install (public API)
set (NLOPT_HEADERS
  api/nlopt.h
  "${CMAKE_BINARY_DIR}/api/nlopt.hpp"
  "${CMAKE_BINARY_DIR}/api/nlopt.f"
)

install (
  FILES ${NLOPT_HEADERS}
  DESTINATION "${RELATIVE_INSTALL_INCLUDE_DIR}")

#==============================================================================
# MAN PAGES AND DOCS
#==============================================================================

if (UNIX)
  install (
    FILES api/nlopt.3 api/nlopt_minimize.3 api/nlopt_minimize_constrained.3
    DESTINATION "${RELATIVE_INSTALL_MAN_DIR}/man3")
endif ()

install (
  FILES AUTHORS COPYING NEWS
  DESTINATION "${RELATIVE_INSTALL_DOC_DIR}")
install (
  FILES README.md
  RENAME README
  DESTINATION "${RELATIVE_INSTALL_DOC_DIR}")

#==============================================================================
# DEV CMAKE GENERATION
#==============================================================================

# add library target to the build-tree export set
export (
  TARGETS ${nlopt_lib}
  FILE "${PROJECT_BINARY_DIR}/NLoptLibraryDepends.cmake")

# install the export set for use from the install-tree
install (
  EXPORT NLoptLibraryDepends
  DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}"
  COMPONENT "Development")

# create cmake config files for build-tree and install-tree
set (NLOPT_LIBRARIES ${nlopt_lib})
file (RELATIVE_PATH RELATIVE_NLOPT_INCLUDE_DIRS
  "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
file (RELATIVE_PATH RELATIVE_NLOPT_LIB_DIR
  "${INSTALL_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
configure_file (cmake/NLoptConfig.cmake.in NLoptConfig.cmake @ONLY)
configure_file (cmake/NLoptConfigVersion.cmake.in NLoptConfigVersion.cmake @ONLY)

# install package configuration file and package version file
install (
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/NLoptConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/NLoptConfigVersion.cmake"
  DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}"
  COMPONENT "Development")

#==============================================================================
# DEV PKGCONFIG GENERATION
#==============================================================================

# create and install pkgconfig file
if (UNIX OR MINGW)
  configure_file (nlopt.pc.in nlopt.pc @ONLY)

  install (
    FILES "${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc"
    DESTINATION "${RELATIVE_INSTALL_LIB_DIR}/pkgconfig"
    COMPONENT "Development")
endif ()

#==============================================================================
# PYTHON/GUILE BINDINGS (SWIG)
#==============================================================================

if (BUILD_PYTHON)
  find_package (PythonInterp)
  find_package (PythonLibs)
  find_package (NumPy)
endif ()

if (BUILD_GUILE)
  find_package (Guile)
endif ()

if (WITH_SWIG)
  find_package (SWIG)
endif ()

# python module installation directory
if (PYTHONINTERP_FOUND)
  execute_process (
    COMMAND "${PYTHON_EXECUTABLE}" -c
"from distutils.sysconfig import get_python_lib
print(get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  get_filename_component (_ABS_PYTHON_MODULE_PATH
    "${_ABS_PYTHON_MODULE_PATH}" ABSOLUTE)
  file (RELATIVE_PATH _REL_PYTHON_MODULE_PATH
    "${CMAKE_INSTALL_PREFIX}" "${_ABS_PYTHON_MODULE_PATH}")
  set (PYTHON_MODULE_PATH "${_REL_PYTHON_MODULE_PATH}")
endif ()

# generate wrappers + Python/Guile modules targets
add_subdirectory (swig)

#==============================================================================
# OCTAVE/MATLAB BINDINGS
#==============================================================================

if (BUILD_OCTAVE)
  find_package (Octave)
endif ()

if (BUILD_MATLAB)
  find_package (Matlab COMPONENTS MX_LIBRARY MAIN_PROGRAM)
endif ()

# oct/mex module targets
if (OCTAVE_FOUND OR Matlab_FOUND)
  add_subdirectory (octave)
endif ()

#==============================================================================
# CTEST
#==============================================================================

# tests: C, Python, Guile, Octave
if (BUILD_TESTS)
  enable_testing ()
  add_subdirectory (test)
endif ()

#==============================================================================
# CPACK
#==============================================================================

set (CPACK_PACKAGE_NAME                "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_CONTACT             "Steven G. Johnson <stevenj@alum.mit.edu>")
set (CPACK_PACKAGE_VENDOR              "Massachusetts Institute of Technology")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "NLopt: a nonlinear optimization library")
set (CPACK_PACKAGE_DESCRIPTION_FILE    "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_README        "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_LICENSE       "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set (CPACK_PACKAGE_VERSION_MAJOR       "${NLOPT_MAJOR_VERSION}")
set (CPACK_PACKAGE_VERSION_MINOR       "${NLOPT_MINOR_VERSION}")
set (CPACK_PACKAGE_VERSION_PATCH       "${NLOPT_BUGFIX_VERSION}")
if (NOT WIN32)
  set (CPACK_GENERATOR                 "TGZ")
  set (CPACK_SOURCE_GENERATOR          "TGZ")
  set (CPACK_BINARY_TGZ                "ON")
  set (CPACK_SOURCE_TGZ                "ON")
endif ()
set (CPACK_SOURCE_IGNORE_FILES         ".git;.*~;build/;${CPACK_SOURCE_IGNORE_FILES}")
set (CPACK_SOURCE_PACKAGE_FILE_NAME    "${CPACK_PACKAGE_NAME}-${NLOPT_VERSION_STRING}")

include (CPack)
