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

cmake_minimum_required (VERSION 2.8.9)

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
list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# provide user options
option (BUILD_SHARED_LIBS       "Build NLopt as a shared library" ON)
option (WITH_CXX                "Enable C++-based routines (StoGO)" OFF)
option (WITH_SWIG               "Use SWIG to generate bindings" ON)
option (BUILD_CPP               "Build C++ interace (nlopt.hpp)" ON)
option (BUILD_F77               "Build Fortran iterface" ON)
option (BUILD_PYTHON            "Build Python bindings (using SWIG)" ON)
option (BUILD_GUILE             "Build Guile bindings (using SWIG)" ON)
option (BUILD_OCTAVE            "Build Octave bindings" ON)
option (BUILD_MATLAB            "Build MATLAB bindings" ON)
option (BUILD_TESTS             "Build tests" ON)
option (BUILD_TESTS_UTIL        "Build util tests" OFF)
option (BUILD_TESTS_ALG         "Build algorithms tests" OFF)
option (WITH_TESTS_INSTALL_TREE "Run tests from install-tree or build-tree" OFF)
option (WITH_EXTRA_WARNINGS     "Enable extra compiler warnings (GNU GCC)" OFF)
option (WITH_FOLDERS            "Enable grouping of files/targets in folders" ON)

# validate options combination
if (BUILD_PYTHON OR BUILD_GUILE)
  if (NOT WITH_SWIG)
    message (WARNING "Python/Guile wrappers require SWIG")
  endif ()
  if (NOT BUILD_CPP)
    message (WARNING "Python/Guile wrappers require the C++ interface (nlopt.hpp)")
  endif ()
endif ()

# enable C++ language compiler if needed
if (WITH_CXX OR BUILD_PYTHON OR BUILD_GUILE OR BUILD_OCTAVE)
  enable_language (CXX)
endif ()

# default installation directories
# (provides variables: CMAKE_INSTALL_<dir> and CMAKE_INSTALL_FULL_<dir>)
include (GNUInstallDirs)

# set RPATH of installed executables and libraries
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if (POLICY CMP0042)
  # Set MACOSX_RPATH to ON
  cmake_policy (SET CMP0042 NEW)
endif ()

# organize target projects inside IDE
# (group cmake generated targets like ALL_BUILD/ZERO_CHECK/..)
if (WITH_FOLDERS)
  set_property (GLOBAL PROPERTY USE_FOLDERS ON)
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

# library name
set (NLOPT_SUFFIX)
if (WITH_CXX)
  set (NLOPT_SUFFIX _cxx)
endif ()
set (NLOPT_LIB nlopt${NLOPT_SUFFIX})

# lists of source files
set (NLOPT_SOURCES_DIRECT direct/DIRect.c direct/direct_wrap.c
  direct/DIRserial.c direct/DIRsubrout.c direct/direct-internal.h
  direct/direct.h)
set (NLOPT_SOURCES_CDIRECT cdirect/cdirect.c cdirect/hybrid.c
  cdirect/cdirect.h)
set (NLOPT_SOURCES_STOGO stogo/global.cc stogo/linalg.cc stogo/local.cc
  stogo/tools.cc stogo/stogo.cc stogo/global.h stogo/linalg.h stogo/local.h
  stogo/tools.h stogo/stogo_config.h stogo/stogo.h)
set (NLOPT_SOURCES_PRAXIS praxis/praxis.c praxis/praxis.h)
set (NLOPT_SOURCES_LUKSAN luksan/plis.c luksan/plip.c luksan/pnet.c
  luksan/mssubs.c luksan/pssubs.c luksan/luksan.h)
set (NLOPT_SOURCES_CRS crs/crs.c crs/crs.h)
set (NLOPT_SOURCES_MLSL mlsl/mlsl.c mlsl/mlsl.h)
set (NLOPT_SOURCES_MMA mma/mma.c mma/mma.h mma/ccsa_quadratic.c)
set (NLOPT_SOURCES_COBYLA cobyla/cobyla.c cobyla/cobyla.h)
set (NLOPT_SOURCES_NEWUOA newuoa/newuoa.c newuoa/newuoa.h)
set (NLOPT_SOURCES_NELDERMEAD neldermead/nldrmd.c neldermead/sbplx.c
  neldermead/neldermead.h)
set (NLOPT_SOURCES_AUGLAG auglag/auglag.c auglag/auglag.h)
set (NLOPT_SOURCES_BOBYQA bobyqa/bobyqa.c bobyqa/bobyqa.h)
set (NLOPT_SOURCES_ISRES isres/isres.c isres/isres.h)
set (NLOPT_SOURCES_ESCH esch/esch.c esch/esch.h)
set (NLOPT_SOURCES_SLSQP slsqp/slsqp.c slsqp/slsqp.h)
set (NLOPT_SOURCES_CQUAD cquad/cquad.c cquad/cquad.h)
set (NLOPT_SOURCES_SUBPLEX subplex/subplex.c subplex/subplex.h)
set (NLOPT_SOURCES_TENSOR tensor/tensor.c)
set (NLOPT_SOURCES_API api/general.c api/options.c api/optimize.c
  api/deprecated.c api/nlopt-internal.h api/nlopt.h)
set (NLOPT_SOURCES_F77 api/f77api.c api/f77funcs.h api/f77funcs_.h)
set (NLOPT_SOURCES_UTIL util/mt19937ar.c util/sobolseq.c util/timer.c
  util/stop.c util/redblack.c util/qsort_r.c util/rescale.c util/soboldata.h
  util/nlopt-util.h util/redblack.h)

# source files: C
set (NLOPT_SOURCES
  ${NLOPT_SOURCES_DIRECT}
  ${NLOPT_SOURCES_CDIRECT}
  ${NLOPT_SOURCES_PRAXIS}
  ${NLOPT_SOURCES_LUKSAN}
  ${NLOPT_SOURCES_CRS}
  ${NLOPT_SOURCES_MLSL}
  ${NLOPT_SOURCES_MMA}
  ${NLOPT_SOURCES_COBYLA}
  ${NLOPT_SOURCES_NEWUOA}
  ${NLOPT_SOURCES_NELDERMEAD}
  ${NLOPT_SOURCES_AUGLAG}
  ${NLOPT_SOURCES_BOBYQA}
  ${NLOPT_SOURCES_ISRES}
  ${NLOPT_SOURCES_ESCH}
  ${NLOPT_SOURCES_SLSQP}
  #${NLOPT_SOURCES_CQUAD}
  #${NLOPT_SOURCES_SUBPLEX}
  #${NLOPT_SOURCES_TENSOR}
  ${NLOPT_SOURCES_API}
  ${NLOPT_SOURCES_UTIL}
)

# source files: Fortran interface
if (BUILD_F77)
  list (APPEND NLOPT_SOURCES ${NLOPT_SOURCES_F77})
endif ()

# source files: C++
if (WITH_CXX)
  list (APPEND NLOPT_SOURCES ${NLOPT_SOURCES_STOGO})
endif ()

# include dirs
set (${INCLUDE_DIRECTORIES} "")
include_directories (
  ${CMAKE_BINARY_DIR}/api
  ${CMAKE_BINARY_DIR}
  api util
  direct cdirect praxis luksan crs mlsl mma cobyla newuoa neldermead
  auglag bobyqa isres esch slsqp stogo
  #cquad subplex tensor
)

# dllimport/dllexport decorations for Windows
if (BUILD_SHARED_LIBS AND WIN32)
  add_definitions (-DNLOPT_DLL -DNLOPT_DLL_EXPORT)
endif ()

# library target (Shared or Static)
add_library (${NLOPT_LIB} ${NLOPT_SOURCES})

# link against libm if needed
if (UNIX OR HAVE_LIBM)
  target_link_libraries (${NLOPT_LIB} m)
endif ()

# if built as static lib, pass -fPIC in case building extension modules
if (NOT BUILD_SHARED_LIBS AND
    (BUILD_PYTHON OR BUILD_GUILE OR BUILD_OCTAVE OR BUILD_MATLAB))
  set_target_properties (${NLOPT_LIB} PROPERTIES POSITION_INDEPENDENT_CODE True)
endif ()

# library version
set_target_properties (${NLOPT_LIB} PROPERTIES
  VERSION "${NLOPT_SOVERSION_STRING}"
  SOVERSION ${NLOPT_SOVERSION_INTERFACE})

# group files and target projects in IDE
if (WITH_FOLDERS)
  set_target_properties (${NLOPT_LIB} PROPERTIES FOLDER "Libraries")
  source_group (direct     FILES ${NLOPT_SOURCES_DIRECT})
  source_group (cdirect    FILES ${NLOPT_SOURCES_CDIRECT})
  source_group (praxis     FILES ${NLOPT_SOURCES_PRAXIS})
  source_group (luksan     FILES ${NLOPT_SOURCES_LUKSAN})
  source_group (crs        FILES ${NLOPT_SOURCES_CRS})
  source_group (mlsl       FILES ${NLOPT_SOURCES_MLSL})
  source_group (mma        FILES ${NLOPT_SOURCES_MMA})
  source_group (cobyla     FILES ${NLOPT_SOURCES_COBYLA})
  source_group (newuoa     FILES ${NLOPT_SOURCES_NEWUOA})
  source_group (neldermead FILES ${NLOPT_SOURCES_NELDERMEAD})
  source_group (auglag     FILES ${NLOPT_SOURCES_AUGLAG})
  source_group (bobyqa     FILES ${NLOPT_SOURCES_BOBYQA})
  source_group (isres      FILES ${NLOPT_SOURCES_ISRES})
  source_group (esch       FILES ${NLOPT_SOURCES_ESCH})
  source_group (slsqp      FILES ${NLOPT_SOURCES_SLSQP})
  source_group (cquad      FILES ${NLOPT_SOURCES_CQUAD})
  source_group (subplex    FILES ${NLOPT_SOURCES_SUBPLEX})
  source_group (tensor     FILES ${NLOPT_SOURCES_TENSOR})
  source_group (api        FILES ${NLOPT_SOURCES_API})
  source_group (util       FILES ${NLOPT_SOURCES_UTIL})
  source_group (f77        FILES ${NLOPT_SOURCES_F77})
  source_group (stogo      FILES ${NLOPT_SOURCES_STOGO})
endif ()

# library installation
install (
  TARGETS ${NLOPT_LIB}
  EXPORT  NLoptLibraryDepends
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT "Libraries"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT "Libraries"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT "Libraries")

#==============================================================================
# C/C++/FORTRAN HEADERS
#==============================================================================

# generate nlopt.hpp and nlopt.f
if (BUILD_CPP OR BUILD_F77)
  add_subdirectory (api)
endif ()

# header files to install (public API)
set (NLOPT_HEADERS api/nlopt.h)
if (BUILD_CPP)
  list (APPEND NLOPT_HEADERS "${CMAKE_BINARY_DIR}/api/nlopt.hpp")
endif ()
if (BUILD_F77)
  list (APPEND NLOPT_HEADERS "${CMAKE_BINARY_DIR}/api/nlopt.f")
endif ()

install (
  FILES ${NLOPT_HEADERS}
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  COMPONENT "Headers")

#==============================================================================
# DEV CMAKE GENERATION
#==============================================================================

# add library target to the build-tree export set
export (
  TARGETS ${NLOPT_LIB}
  FILE "${PROJECT_BINARY_DIR}/NLoptLibraryDepends.cmake")

# install the export set for use from the install-tree
install (
  EXPORT NLoptLibraryDepends
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nlopt"
  COMPONENT "Development")

# create cmake config files for build-tree and install-tree
file (RELATIVE_PATH RELATIVE_NLOPT_INCLUDE_DIRS
  "${CMAKE_INSTALL_FULL_LIBDIR}/cmake/nlopt" "${CMAKE_INSTALL_FULL_INCLUDEDIR}")
file (RELATIVE_PATH RELATIVE_NLOPT_LIB_DIR
  "${CMAKE_INSTALL_FULL_LIBDIR}/cmake/nlopt" "${CMAKE_INSTALL_FULL_LIBDIR}")
configure_file (cmake/NLoptConfig.cmake.in NLoptConfig.cmake @ONLY)
configure_file (cmake/NLoptConfigVersion.cmake.in NLoptConfigVersion.cmake @ONLY)

# install package configuration file and package version file
install (
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/NLoptConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/NLoptConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nlopt"
  COMPONENT "Development")

#==============================================================================
# DEV PKGCONFIG GENERATION
#==============================================================================

# create and install pkgconfig file
configure_file (nlopt.pc.in nlopt.pc @ONLY)

install (
  FILES "${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
  COMPONENT "Development")

#==============================================================================
# MAN PAGES AND DOCS
#==============================================================================

install (
  FILES api/nlopt.3 api/nlopt_minimize.3 api/nlopt_minimize_constrained.3
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man3"
  COMPONENT "Documentation")

install (
  FILES AUTHORS COPYING NEWS
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
  COMPONENT "Documentation")
install (
  FILES README.md
  RENAME README
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
  COMPONENT "Documentation")

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
    OUTPUT_VARIABLE INSTALL_PYTHON_FULL_LIB_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  get_filename_component (INSTALL_PYTHON_FULL_LIB_DIR
    "${INSTALL_PYTHON_FULL_LIB_DIR}" ABSOLUTE)
  file (RELATIVE_PATH INSTALL_PYTHON_LIB_DIR
    "${CMAKE_INSTALL_PREFIX}" "${INSTALL_PYTHON_FULL_LIB_DIR}")
endif ()

# guile extension installation directories
if (GUILE_FOUND)
  file (RELATIVE_PATH INSTALL_GUILE_EXT_DIR
    "${GUILE_ROOT_DIR}" "${GUILE_EXTENSION_DIR}")
  file (RELATIVE_PATH INSTALL_GUILE_SITE_DIR
    "${GUILE_ROOT_DIR}" "${GUILE_SITE_DIR}")
  set (INSTALL_GUILE_FULL_EXT_DIR
    "${CMAKE_INSTALL_PREFIX}/${INSTALL_GUILE_EXT_DIR}")
  set (INSTALL_GUILE_FULL_SITE_DIR
    "${CMAKE_INSTALL_PREFIX}/${INSTALL_GUILE_SITE_DIR}")
endif ()

# generate wrappers + Python/Guile modules targets
if (WITH_SWIG OR BUILD_GUILE OR BUILD_PYTHON)
  add_subdirectory (swig)
endif ()

#==============================================================================
# OCTAVE/MATLAB BINDINGS
#==============================================================================

if (BUILD_OCTAVE)
  find_package (Octave)
endif ()

if (BUILD_MATLAB)
  find_package (Matlab COMPONENTS MX_LIBRARY MAIN_PROGRAM)
endif ()

# octave package installation directories
if (OCTAVE_FOUND)
  file (RELATIVE_PATH INSTALL_OCTAVE_OCT_DIR
    "${OCTAVE_ROOT_DIR}" "${OCTAVE_OCT_SITE_DIR}")
  file (RELATIVE_PATH INSTALL_OCTAVE_M_DIR
    "${OCTAVE_ROOT_DIR}" "${OCTAVE_M_SITE_DIR}")
  set (INSTALL_OCTAVE_FULL_OCT_DIR
    "${CMAKE_INSTALL_PREFIX}/${INSTALL_OCTAVE_OCT_DIR}")
  set (INSTALL_OCTAVE_FULL_M_DIR
    "${CMAKE_INSTALL_PREFIX}/${INSTALL_OCTAVE_M_DIR}")
endif ()

# matlab package installation directories
if (Matlab_FOUND)
  set (INSTALL_MATLAB_DIR "matlab")
  set (INSTALL_MATLAB_FULL_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_MATLAB_DIR}")
endif ()

# oct/mex module targets
if (BUILD_OCTAVE OR BUILD_MATLAB)
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

# customize binary/source installers
set (CPACK_PACKAGE_NAME                "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_CONTACT             "Steven G. Johnson <stevenj@alum.mit.edu>")
set (CPACK_PACKAGE_VENDOR              "Massachusetts Institute of Technology")
set (CPACK_PACKAGE_VERSION_MAJOR       "${NLOPT_MAJOR_VERSION}")
set (CPACK_PACKAGE_VERSION_MINOR       "${NLOPT_MINOR_VERSION}")
set (CPACK_PACKAGE_VERSION_PATCH       "${NLOPT_BUGFIX_VERSION}")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "NLopt: a nonlinear optimization library")
set (CPACK_PACKAGE_DESCRIPTION_FILE    "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_README        "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_LICENSE       "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set (CPACK_SOURCE_IGNORE_FILES         ".git/;build/;.*~;${CPACK_SOURCE_IGNORE_FILES}")
set (CPACK_SOURCE_PACKAGE_FILE_NAME    "${CPACK_PACKAGE_NAME}-${NLOPT_VERSION_STRING}")
if (WIN32)
  set (CPACK_GENERATOR                 "ZIP")  # NSIS
  set (CPACK_SOURCE_GENERATOR          "ZIP")
  set (CPACK_BINARY_ZIP                "ON")
  set (CPACK_SOURCE_ZIP                "ON")
else ()
  set (CPACK_GENERATOR                 "TGZ")
  set (CPACK_SOURCE_GENERATOR          "TGZ")
  set (CPACK_BINARY_TGZ                "ON")
  set (CPACK_SOURCE_TGZ                "ON")
endif ()

# generates package/package_source targets
include (CPack)

#==============================================================================
# SUMMARY
#==============================================================================

# for debugging purposes only
function (show_variables)
    get_cmake_property (_varnames VARIABLES)
    foreach (_varname ${_varnames})
        message (STATUS "${_varname}=${${_varname}}")
    endforeach ()
endfunction ()

message (STATUS "")
message (STATUS "====== ${CMAKE_PROJECT_NAME} ${NLOPT_VERSION_STRING} ======")
message (STATUS "  Platform:")
message (STATUS "    System: "           ${CMAKE_SYSTEM})
message (STATUS "    CMake: "            ${CMAKE_VERSION})
message (STATUS "    CMake generator: "  ${CMAKE_GENERATOR})
message (STATUS "    CMake build tool: " ${CMAKE_BUILD_TOOL})
message (STATUS "  Compiler:")
if (WITH_CXX)
  message (STATUS "    C++ Compiler: ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}")
  message (STATUS "    C++ Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")
  message (STATUS "    C++ flags (Release): ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  message (STATUS "    C++ flags (Debug): ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif ()
message (STATUS "    C Compiler: ${CMAKE_C_COMPILER} ${CMAKE_C_COMPILER_ARG1}")
message (STATUS "    C Compiler Version: ${CMAKE_C_COMPILER_VERSION}")
message (STATUS "    C flags (Release): ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
message (STATUS "    C flags (Debug): ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
message (STATUS "    Linker: ${CMAKE_LINKER} ${CMAKE_LINKER_ARG1}")
if (BUILD_SHARED_LIBS)
  message (STATUS "    Linker flags (Release): ${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_RELEASE}")
  message (STATUS "    Linker flags (Debug): ${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
else ()
  message (STATUS "    Linker flags (Release): ${CMAKE_STATIC_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS_RELEASE}")
  message (STATUS "    Linker flags (Debug): ${CMAKE_STATIC_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS_DEBUG}")
endif ()
message (STATUS "  Configuration:")
if (BUILD_SHARED_LIBS)
  message (STATUS "    Library: Shared")
else ()
  message (STATUS "    Library: Static")
endif ()
message (STATUS "    Build: ${CMAKE_BUILD_TYPE}")
message (STATUS "    config.h: ${CMAKE_CURRENT_BINARY_DIR}/config.h")
message (STATUS "    Install prefix: ${CMAKE_INSTALL_PREFIX}")
message (STATUS "  Bindings:")
message (STATUS "    C++: ${BUILD_CPP}")
message (STATUS "    Fortran: ${BUILD_F77}")
if (WITH_SWIG AND SWIG_FOUND)
  message (STATUS "    SWIG: YES (${SWIG_VERSION})")
else ()
  message (STATUS "    SWIG: NO")
endif ()
if (BUILD_PYTHON AND PYTHONINTERP_FOUND AND NUMPY_FOUND)
  message (STATUS "    Python: YES (${PYTHON_VERSION_STRING}), NumPy: YES (${NUMPY_VERSION_STRING})")
else ()
  message (STATUS "    Python: NO")
endif ()
if (BUILD_GUILE AND GUILE_FOUND)
  message (STATUS "    Guile: YES (${GUILE_VERSION_STRING})")
else ()
  message (STATUS "    Guile: NO")
endif ()
if (BUILD_OCTAVE AND OCTAVE_FOUND)
  message (STATUS "    Octave: YES (${OCTAVE_VERSION_STRING})")
else ()
  message (STATUS "    Octave: NO")
endif ()
if (BUILD_MATLAB AND Matlab_FOUND)
  message (STATUS "    MATLAB: YES (${Matlab_VERSION_STRING})")
else ()
  message (STATUS "    MATLAB: NO")
endif ()
message (STATUS "  Testing:")
message (STATUS "    Tests: ${BUILD_TESTS}")
message (STATUS "    Tests Util: ${BUILD_TESTS_UTIL}")
message (STATUS "    Tests Algorithms: ${BUILD_TESTS_ALG}")
if (WITH_TESTS_INSTALL_TREE)
  message (STATUS "    Run from: install-tree")
else ()
  message (STATUS "    Run from: build-tree")
endif ()
message (STATUS "")
