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
cmake_minimum_required (VERSION 3.18)

set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (nlopt)

#==============================================================================
# version
set (NLOPT_MAJOR_VERSION "2")
set (NLOPT_MINOR_VERSION "11")
set (NLOPT_BUGFIX_VERSION "0")
set (NLOPT_VERSION_STRING ${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION}.${NLOPT_BUGFIX_VERSION})
message (STATUS "NLopt version ${NLOPT_VERSION_STRING}")

# This is the ABI version number, which differes from the API version above
# (it indicates ABI compatibility), but they are typically incremented together.
set(SO_MAJOR 1)
set(SO_MINOR 1)
set(SO_PATCH 0)
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option (NLOPT_CXX "enable cxx routines" ON)
option (NLOPT_FORTRAN "enable fortran" OFF)
option (BUILD_SHARED_LIBS "Build NLopt as a shared library" ON)
option (NLOPT_PYTHON "build python bindings" ON)
option (NLOPT_OCTAVE "build octave bindings" ON)
option (NLOPT_MATLAB "build matlab bindings" OFF)
option (NLOPT_GUILE "build guile bindings" ON)
option (NLOPT_JAVA "build java bindings" ON)
option (NLOPT_SWIG "use SWIG to build bindings" ON)
option (NLOPT_LUKSAN "enable LGPL Luksan solvers" ON)
option (NLOPT_TESTS "build unit tests" OFF)

if (NLOPT_FORTRAN)
  enable_language (Fortran)
endif ()

include (GNUInstallDirs)

# Offer the user the choice of overriding the installation directories
set (NLOPT_INSTALL_LIBDIR     ${CMAKE_INSTALL_LIBDIR}        CACHE STRING "Installation directory for libraries")
set (NLOPT_INSTALL_BINDIR     ${CMAKE_INSTALL_BINDIR}        CACHE STRING "Installation directory for executables")
set (NLOPT_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}    CACHE STRING "Installation directory for header files")
set (NLOPT_INSTALL_DATADIR    ${CMAKE_INSTALL_DATADIR}/nlopt CACHE STRING "Installation directory for data files")
set (NLOPT_INSTALL_MANDIR     ${CMAKE_INSTALL_MANDIR}        CACHE STRING "Installation directory for man documentation")
set (NLOPT_INSTALL_CMAKEDIR   ${CMAKE_INSTALL_LIBDIR}/cmake/nlopt CACHE STRING "Installation directory for cmake config files")

set (CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${NLOPT_INSTALL_LIBDIR} CACHE PATH "rpath")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE CACHE BOOL "use link path")

include (CheckIncludeFiles)
include (CheckFunctionExists)
include (CheckTypeSize)
include (CheckCSourceCompiles)

#==============================================================================
# COMPILATION CHECKINGS and CONFIGURATION GENERATION
#==============================================================================
check_include_file (getopt.h HAVE_GETOPT_H)
check_include_file (unistd.h HAVE_UNISTD_H)
check_include_file (stdint.h HAVE_STDINT_H)
check_include_file (time.h HAVE_TIME_H)
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
check_function_exists (getopt HAVE_GETOPT)
check_type_size ("uint32_t" SIZEOF_UINT32_T)
set (HAVE_UINT32_T ${SIZEOF_UINT32_T})
check_type_size ("unsigned int" SIZEOF_UNSIGNED_INT)
check_type_size ("unsigned long" SIZEOF_UNSIGNED_LONG)

find_library(M_LIBRARY m)
if (M_LIBRARY)
  set (LIBS_PRIVATE "-lm")
endif()

# BLAS/LAPACK for the algorithms that use it (src/util/nlopt-blas.h): AUTO uses
# one if found, ON requires one, OFF uses the built-in reference routines;
# anything else is the name or path of the library (or a list of them, e.g.
# "blas;lapack") to use.
set (NLOPT_BLAS AUTO CACHE STRING "BLAS/LAPACK: AUTO, ON, OFF, or library name(s)/path(s)")
set_property (CACHE NLOPT_BLAS PROPERTY STRINGS AUTO ON OFF)
# The keywords are case-insensitive.
string (TOUPPER "${NLOPT_BLAS}" _blas_mode)
if (_blas_mode MATCHES "^(1|ON|YES|TRUE|Y)$")
  set (_blas_mode ON)
elseif (_blas_mode MATCHES "^(0|OFF|NO|FALSE|N|)$")
  set (_blas_mode OFF)
elseif (NOT _blas_mode STREQUAL "AUTO")
  set (_blas_mode LIBRARY)
endif ()
set (NLOPT_BLAS_LIBRARIES)
set (_blas_candidates)
if (NOT _blas_mode STREQUAL "OFF")
  if (_blas_mode STREQUAL "AUTO" OR _blas_mode STREQUAL "ON")
    find_package (LAPACK QUIET)
    if (LAPACK_FOUND)
      set (_blas_candidates "${LAPACK_LIBRARIES}")
    endif ()
  else ()
    set (_blas_candidates "${NLOPT_BLAS}")
  endif ()
endif ()
if (_blas_candidates)
  # Find the symbol naming convention under which all routines in nlopt-blas.h
  # link, preferring the interface with integers of the platform's pointer size
  # for libraries that provide both (e.g. MKL, libblastrampoline).
  string (MD5 _blas_hash "${_blas_candidates}")
  set (CMAKE_REQUIRED_LIBRARIES ${_blas_candidates})
  if (M_LIBRARY)
    list (APPEND CMAKE_REQUIRED_LIBRARIES ${M_LIBRARY})
  endif ()
  set (CMAKE_REQUIRED_QUIET ON)
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set (_blas_suffixes _ "" _64_ _64)
  else ()
    set (_blas_suffixes _64_ _64 _ "")
  endif ()
  foreach (_suffix IN LISTS _blas_suffixes)
    set (_decls)
    set (_calls)
    foreach (_r daxpy dscal dnrm2 ddot dgemv dtpmv dtpsv dtrsm dtrsv dgelsy dgeqr2 dgerq2
                dlaic1 dlarf dlarfgp dlartgp dorm2r dormr2)
      string (APPEND _decls "char ${_r}${_suffix}(void);\n")
      string (APPEND _calls "  ${_r}${_suffix}();\n")
    endforeach ()
    string (MAKE_C_IDENTIFIER "NLOPT_BLAS_${_blas_hash}${_suffix}" _var)
    check_c_source_compiles ("${_decls}int main(void) {\n${_calls}  return 0;\n}\n" ${_var})
    if (${_var})
      set (NLOPT_BLAS_LIBRARIES ${_blas_candidates})
      set (NLOPT_BLAS_SUFFIX "${_suffix}")
      break ()
    endif ()
  endforeach ()
  unset (CMAKE_REQUIRED_LIBRARIES)
  unset (CMAKE_REQUIRED_QUIET)
  if (NOT NLOPT_BLAS_LIBRARIES AND NOT _blas_mode STREQUAL "AUTO")
    message (FATAL_ERROR "NLOPT_BLAS: cannot link the BLAS/LAPACK routines needed by NLopt from ${_blas_candidates}")
  elseif (NOT NLOPT_BLAS_LIBRARIES)
    message (STATUS "Cannot link the BLAS/LAPACK routines needed by NLopt from ${_blas_candidates}")
  endif ()
elseif (_blas_mode STREQUAL "ON")
  message (FATAL_ERROR "NLOPT_BLAS: BLAS/LAPACK not found")
endif ()
if (NLOPT_BLAS_LIBRARIES)
  # Integer size from the symbol names: "64" suffixes denote 64-bit integers.
  if (NLOPT_BLAS_SUFFIX MATCHES "64")
    set (NLOPT_BLAS_INT int64_t)
  else ()
    set (NLOPT_BLAS_INT int)
  endif ()
  message (STATUS "BLAS/LAPACK: ${NLOPT_BLAS_LIBRARIES} (symbol suffix \"${NLOPT_BLAS_SUFFIX}\", ${NLOPT_BLAS_INT} integers)")
  # For nlopt.pc: paths and flags as they are, bare names as -l<name>.
  foreach (_lib IN LISTS NLOPT_BLAS_LIBRARIES)
    if (_lib MATCHES "^-" OR _lib MATCHES "/")
      set (LIBS_PRIVATE "${LIBS_PRIVATE} ${_lib}")
    else ()
      set (LIBS_PRIVATE "${LIBS_PRIVATE} -l${_lib}")
    endif ()
  endforeach ()
else ()
  message (STATUS "BLAS/LAPACK: built-in reference routines")
endif ()

if (NOT DEFINED HAVE_FPCLASSIFY)
  message(STATUS "Looking for fpclassify")
  file (WRITE ${PROJECT_BINARY_DIR}/fpclassify.c "#include <math.h>\n")
  file (APPEND ${PROJECT_BINARY_DIR}/fpclassify.c "int main(void) {\n")
  file (APPEND ${PROJECT_BINARY_DIR}/fpclassify.c "if (!fpclassify(3.14159)) fpclassify(2.7183);\n")
  file (APPEND ${PROJECT_BINARY_DIR}/fpclassify.c "  return 0; }\n")
  try_compile (HAVE_FPCLASSIFY
  ${PROJECT_BINARY_DIR}/build_fpclassify
  ${PROJECT_BINARY_DIR}/fpclassify.c
  CMAKE_FLAGS -DLINK_LIBRARIES=m)
  message(STATUS "Looking for fpclassify - ${HAVE_FPCLASSIFY}")
endif ()

option (WITH_THREADLOCAL "check thread local keyword" ON)
if (WITH_THREADLOCAL AND NOT DEFINED THREADLOCAL)
  foreach (_THREADLOCAL_KEY "__thread" "__declspec(thread)")
    unset (HAVE_THREAD_LOCAL_STORAGE CACHE)
    check_c_source_compiles("
    ${_THREADLOCAL_KEY} int tls;

    int main(void) {
        return 0;
    }" HAVE_THREAD_LOCAL_STORAGE)
    if (HAVE_THREAD_LOCAL_STORAGE)
      set (THREADLOCAL ${_THREADLOCAL_KEY} CACHE STRING "Thread local keyword")
    endif ()
  endforeach()
endif ()

#==============================================================================
# CREATE nlopt_config.h
#==============================================================================

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt_config.h IMMEDIATE)

# pkgconfig file
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt.pc.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc @ONLY)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc DESTINATION ${NLOPT_INSTALL_LIBDIR}/pkgconfig)

#==============================================================================
# nlopt LIBRARY TARGET (SHARED OR STATIC)
#==============================================================================

set (NLOPT_HEADERS
  src/api/nlopt.h ${PROJECT_BINARY_DIR}/nlopt.hpp
)

add_custom_command (OUTPUT nlopt.hpp
  COMMAND ${CMAKE_COMMAND} -DAPI_SOURCE_DIR=${PROJECT_SOURCE_DIR}/src/api -P ${PROJECT_SOURCE_DIR}/cmake/generate-cpp.cmake
  MAIN_DEPENDENCY src/api/nlopt-in.hpp)
add_custom_target (generate-cpp ALL DEPENDS nlopt.hpp)

if (NLOPT_FORTRAN)
  add_custom_command (OUTPUT nlopt.f
    COMMAND ${CMAKE_COMMAND} -DAPI_SOURCE_DIR=${PROJECT_SOURCE_DIR}/src/api -P ${PROJECT_SOURCE_DIR}/cmake/generate-fortran.cmake
    MAIN_DEPENDENCY src/api/nlopt.h)
  add_custom_target (generate-fortran ALL DEPENDS nlopt.f)
  list (APPEND NLOPT_HEADERS ${PROJECT_BINARY_DIR}/nlopt.f)
endif()

set (NLOPT_SOURCES
  src/algs/direct/DIRect.c src/algs/direct/direct_wrap.c src/algs/direct/DIRserial.c src/algs/direct/DIRsubrout.c src/algs/direct/direct-internal.h src/algs/direct/direct.h
  src/algs/cdirect/cdirect.c src/algs/cdirect/hybrid.c src/algs/cdirect/cdirect.h
  src/algs/praxis/praxis.c src/algs/praxis/praxis.h
  src/algs/crs/crs.c src/algs/crs/crs.h
  src/algs/mlsl/mlsl.c src/algs/mlsl/mlsl.h
  src/algs/mma/mma.c src/algs/mma/mma.h src/algs/mma/ccsa_quadratic.c
  src/algs/cobyla/cobyla.c src/algs/cobyla/cobyla.h
  src/algs/newuoa/newuoa.c src/algs/newuoa/newuoa.h
  src/algs/neldermead/nldrmd.c src/algs/neldermead/neldermead.h src/algs/neldermead/sbplx.c
  src/algs/auglag/auglag.c src/algs/auglag/auglag.h
  src/algs/bobyqa/bobyqa.c src/algs/bobyqa/bobyqa.h
  src/algs/isres/isres.c src/algs/isres/isres.h
  src/algs/slsqp/slsqp.c src/algs/slsqp/slsqp.h
  src/algs/esch/esch.c src/algs/esch/esch.h
  src/api/general.c src/api/options.c src/api/optimize.c src/api/deprecated.c src/api/nlopt-internal.h src/api/nlopt.h src/api/f77api.c src/api/f77funcs.h src/api/f77funcs_.h ${PROJECT_BINARY_DIR}/nlopt.hpp
  src/util/nlopt-blas.h src/util/nlopt-blas-builtin.h src/util/mt19937ar.c src/util/sobolseq.c src/util/soboldata.h src/util/timer.c src/util/stop.c src/util/nlopt-util.h src/util/redblack.c src/util/redblack.h src/util/qsort_r.c src/util/rescale.c
)

if(NLOPT_LUKSAN)
  list(APPEND NLOPT_SOURCES
    src/algs/luksan/plis.c src/algs/luksan/plip.c src/algs/luksan/pnet.c src/algs/luksan/mssubs.c src/algs/luksan/pssubs.c src/algs/luksan/luksan.h)
endif()


set_property(SOURCE src/algs/bobyqa/bobyqa.c src/algs/cdirect/hybrid.c src/algs/mma/ccsa_quadratic.c src/algs/cobyla/cobyla.c
                    src/util/redblack.c src/algs/neldermead/nldrmd.c src/algs/newuoa/newuoa.c src/util/qsort_r.c src/util/mt19937ar.c
             PROPERTY SKIP_UNITY_BUILD_INCLUSION ON)

# workaround for false positive in newuoa
if (CMAKE_C_COMPILER_ID MATCHES "GNU" AND CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
  set_source_files_properties(src/algs/newuoa/newuoa.c PROPERTIES COMPILE_FLAGS -U_FORTIFY_SOURCE)
endif ()

if (NLOPT_CXX)
  list (APPEND NLOPT_SOURCES
    src/algs/stogo/global.cc
    src/algs/stogo/linalg.cc
    src/algs/stogo/local.cc
    src/algs/stogo/stogo.cc
    src/algs/stogo/tools.cc
    src/algs/stogo/global.h
    src/algs/stogo/linalg.h
    src/algs/stogo/local.h
    src/algs/stogo/stogo_config.h
    src/algs/stogo/stogo.h
    src/algs/stogo/tools.h
    src/algs/ags/data_types.hpp
    src/algs/ags/evolvent.hpp
    src/algs/ags/evolvent.cc
    src/algs/ags/solver.hpp
    src/algs/ags/solver.cc
    src/algs/ags/local_optimizer.hpp
    src/algs/ags/local_optimizer.cc
    src/algs/ags/ags.h
    src/algs/ags/ags.cc)

set_property(SOURCE src/algs/ags/solver.cc src/algs/ags/local_optimizer.cc src/algs/ags/ags.cc src/algs/slsqp/slsqp.c
             PROPERTY SKIP_UNITY_BUILD_INCLUSION ON)
endif ()

install (FILES ${NLOPT_HEADERS} DESTINATION ${NLOPT_INSTALL_INCLUDEDIR})

set (nlopt_lib nlopt)
add_library (${nlopt_lib} ${NLOPT_SOURCES})
add_dependencies(${nlopt_lib} generate-cpp)
if (NLOPT_LUKSAN)
  target_include_directories(${nlopt_lib} PRIVATE src/algs/luksan)
  target_compile_definitions (${nlopt_lib} PRIVATE NLOPT_LUKSAN)
endif ()
if (M_LIBRARY)
  target_link_libraries (${nlopt_lib} PRIVATE ${M_LIBRARY})
endif ()
if (NLOPT_BLAS_LIBRARIES)
  target_link_libraries (${nlopt_lib} PRIVATE ${NLOPT_BLAS_LIBRARIES})
  target_compile_definitions (${nlopt_lib} PRIVATE NLOPT_HAVE_BLAS
                              NLOPT_BLAS_SUFFIX=${NLOPT_BLAS_SUFFIX} NLOPT_BLAS_INT=${NLOPT_BLAS_INT})
endif ()
set_target_properties (${nlopt_lib} PROPERTIES SOVERSION ${SO_MAJOR})
set_target_properties (${nlopt_lib} PROPERTIES VERSION "${SO_MAJOR}.${SO_MINOR}.${SO_PATCH}")
if (NLOPT_CXX)
  target_compile_features (${nlopt_lib} PUBLIC cxx_lambdas)
endif ()

#==============================================================================
# INCLUDE DIRECTORIES
#==============================================================================
target_include_directories (${nlopt_lib} PRIVATE
  ${PROJECT_BINARY_DIR}
  src/algs/stogo
  src/algs/ags
  src/util
  src/algs/direct
  src/algs/cdirect
  src/algs/praxis
  src/algs/crs
  src/algs/mlsl
  src/algs/mma
  src/algs/cobyla
  src/algs/newuoa
  src/algs/neldermead
  src/algs/auglag
  src/algs/bobyqa
  src/algs/isres
  src/algs/slsqp
  src/algs/esch
  src/api)

get_target_property (NLOPT_PRIVATE_INCLUDE_DIRS ${nlopt_lib} INCLUDE_DIRECTORIES)
target_include_directories (${nlopt_lib} INTERFACE "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/api;${PROJECT_BINARY_DIR}>" "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

if (BUILD_SHARED_LIBS)
  target_compile_definitions (${nlopt_lib} PUBLIC NLOPT_DLL)

  if (WIN32)
    # needed to switch between dllexport/dllimport
    target_compile_definitions (${nlopt_lib} PRIVATE NLOPT_DLL_EXPORT)
  endif ()
endif ()

# needed for the language plugins
set_target_properties (${nlopt_lib} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# for consistent floating-point operations across architectures
option (DISABLE_FP_CONTRACT "Disable floating-point contractions" ON)
if (DISABLE_FP_CONTRACT AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffp-contract=off")
endif ()

install ( TARGETS ${nlopt_lib}
          EXPORT  NLoptLibraryDepends
          RUNTIME DESTINATION ${NLOPT_INSTALL_BINDIR}
          LIBRARY DESTINATION ${NLOPT_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${NLOPT_INSTALL_LIBDIR}
        )

if (MSVC AND BUILD_SHARED_LIBS)
  install (FILES $<TARGET_PDB_FILE:${nlopt_lib}> DESTINATION ${NLOPT_INSTALL_BINDIR} CONFIGURATIONS Debug RelWithDebInfo COMPONENT Debug)
endif ()

# install man
if (UNIX)
  install (FILES src/api/nlopt.3 src/api/nlopt_minimize.3 src/api/nlopt_minimize_constrained.3 DESTINATION ${NLOPT_INSTALL_MANDIR}/man3)
endif ()

if (NLOPT_PYTHON)
  option (NLOPT_PYTHON_SABI "Use Python stable ABI" OFF)
  if (NLOPT_PYTHON_SABI AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.26)
    find_package (Python 3.6 COMPONENTS Interpreter Development.SABIModule NumPy)
    add_library (Python::Module ALIAS Python::SABIModule)
  else ()
    find_package (Python 3.6 COMPONENTS Interpreter Development.Module NumPy)
  endif ()
endif ()

if (WIN32)
  set (INSTALL_PYTHON_DIR Lib/site-packages CACHE STRING "site packages dir")
else ()
  set (INSTALL_PYTHON_DIR ${NLOPT_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages CACHE STRING "site packages dir")
endif ()

if (NLOPT_GUILE)
  find_package (Guile)
endif ()

if (NLOPT_JAVA)
  # we do not really need any component, only the main JNI target, but if the
  # list of components is left empty, FindJNI defaults to "JVM AWT", and we
  # specifically do not want to check for the AWT library that is not available
  # on headless installations
  find_package (JNI COMPONENTS JVM)
  # FindJNI.cmake in CMake versions prior to 3.24 does not export any targets
  if(JNI_FOUND AND NOT TARGET JNI::JNI)
    add_library(JNI::JNI IMPORTED INTERFACE)
    set_property(TARGET JNI::JNI PROPERTY INTERFACE_INCLUDE_DIRECTORIES
      ${JAVA_INCLUDE_PATH})
    if(NOT JNI_INCLUDE_PATH2_OPTIONAL AND JAVA_INCLUDE_PATH2)
      set_property(TARGET JNI::JNI APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        ${JAVA_INCLUDE_PATH2})
    endif()
  endif()
  find_package (Java 1.5)
endif ()

if (NLOPT_SWIG)
  find_package (SWIG 3)
  if (SWIG_FOUND)
    add_subdirectory (src/swig)
  endif ()
endif ()

if (NLOPT_OCTAVE)
  find_package (Octave)
endif ()

if (NLOPT_MATLAB)
  find_package (Matlab COMPONENTS MX_LIBRARY MAIN_PROGRAM)
endif ()

if (OCTAVE_FOUND OR Matlab_FOUND)
  add_subdirectory (src/octave)
endif ()

if (NLOPT_TESTS)
  enable_testing ()
  add_subdirectory (test)
endif ()

set (CPACK_PACKAGE_NAME          "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_VERSION_MAJOR "${NLOPT_MAJOR_VERSION}")
set (CPACK_PACKAGE_VERSION_MINOR "${NLOPT_MINOR_VERSION}")
set (CPACK_PACKAGE_VERSION_PATCH "${NLOPT_BUGFIX_VERSION}")
set (CPACK_SOURCE_GENERATOR      "TBZ2;TGZ"    )
set (CPACK_BINARY_STGZ           OFF CACHE BOOL "STGZ")
set (CPACK_BINARY_TBZ2           ON  CACHE BOOL "TBZ2")
set (CPACK_BINARY_TGZ            ON  CACHE BOOL "TGZ")
set (CPACK_BINARY_TZ             OFF CACHE BOOL "TZ")
set (CPACK_SOURCE_IGNORE_FILES ".git;/build;.*~;${CPACK_SOURCE_IGNORE_FILES}")

set (CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${NLOPT_VERSION_STRING})

include (CPack)

# configuration files
export (TARGETS ${nlopt_lib} NAMESPACE NLopt:: FILE ${PROJECT_BINARY_DIR}/NLoptLibraryDepends.cmake)

# Install the export set for use with the install-tree
install(EXPORT NLoptLibraryDepends
        NAMESPACE NLopt::
        DESTINATION ${NLOPT_INSTALL_CMAKEDIR}
        COMPONENT Development)

# Create a NLOPTConfig.cmake file for the use from the install tree
# and install it
set (NLOPT_LIBRARIES "NLopt::${nlopt_lib}")

configure_file (cmake/NLoptConfig.cmake.in NLoptConfig.cmake @ONLY)
configure_file (cmake/NLoptConfigVersion.cmake.in NLoptConfigVersion.cmake @ONLY)
install (FILES
          ${CMAKE_CURRENT_BINARY_DIR}/NLoptConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/NLoptConfigVersion.cmake
         DESTINATION ${NLOPT_INSTALL_CMAKEDIR}
         COMPONENT Development)
