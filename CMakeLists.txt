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
cmake_minimum_required (VERSION 3.0)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project (nlopt)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option (WITH_CXX "enable cxx routines" OFF)
option (BUILD_SHARED_LIBS "Build NLopt as a shared library" ON)
option (BUILD_PYTHON "build python bindings" ON)
option (BUILD_OCTAVE "build octave bindings" ON)
option (BUILD_MATLAB "build matlab bindings" ON)
option (BUILD_GUILE "build guile bindings" ON)
option (USE_SWIG "use SWIG to build bindings" ON)

set (NLOPT_SUFFIX)
if (WITH_CXX)
  set (NLOPT_SUFFIX _cxx)
endif ()

include (GNUInstallDirs)

# Offer the user the choice of overriding the installation directories
set (INSTALL_LIB_DIR     ${CMAKE_INSTALL_LIBDIR}        CACHE PATH "Installation directory for libraries")
set (INSTALL_BIN_DIR     ${CMAKE_INSTALL_BINDIR}        CACHE PATH "Installation directory for executables")
set (INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}    CACHE PATH "Installation directory for header files")
set (INSTALL_DATA_DIR    ${CMAKE_INSTALL_DATADIR}/nlopt CACHE PATH "Installation directory for data files")
set (INSTALL_MAN_DIR     ${CMAKE_INSTALL_MANDIR}        CACHE PATH "Installation directory for man documentation")
set (INSTALL_CMAKE_DIR   ${INSTALL_LIB_DIR}/cmake/nlopt CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
  set (var INSTALL_${p}_DIR)
  set (RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
  if (NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif ()
endforeach ()


set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_DIR})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(POLICY CMP0042)
  # Set MACOSX_RPATH to ON
  cmake_policy(SET CMP0042 NEW)
endif()

include (CheckIncludeFiles)
include (CheckFunctionExists)
include (CheckTypeSize)
include (CheckCCompilerFlag)
include (CheckCXXSymbolExists)
include (CheckCXXCompilerFlag)
include (CheckLibraryExists)

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

check_library_exists ("m" sqrt "" HAVE_LIBM)
if (HAVE_LIBM)
  set (M_LIBRARY m)
endif()

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

if (WITH_CXX OR BUILD_PYTHON OR BUILD_GUILE OR BUILD_OCTAVE)
  check_cxx_symbol_exists (_LIBCPP_VERSION string SYSTEM_HAS_LIBCPP)
  if (SYSTEM_HAS_LIBCPP)
    check_cxx_compiler_flag ("-std=c++11" SUPPORTS_STDCXX11)
    if (SUPPORTS_STDCXX11)
      set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    endif ()
  endif ()
endif ()

#==============================================================================
# version
#==============================================================================
set (NLOPT_MAJOR_VERSION "2")
set (NLOPT_MINOR_VERSION "5")
set (NLOPT_BUGFIX_VERSION "0")
set (NLOPT_VERSION_STRING ${NLOPT_MAJOR_VERSION}.${NLOPT_MINOR_VERSION}.${NLOPT_BUGFIX_VERSION})
message (STATUS "NLopt version ${NLOPT_VERSION_STRING}")

#==============================================================================
# CREATE nlopt_config.h
#==============================================================================

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt_config.h IMMEDIATE)

# pkgconfig file
if (UNIX OR MINGW)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt.pc.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc @ONLY)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.pc DESTINATION ${RELATIVE_INSTALL_LIB_DIR}/pkgconfig)
endif ()

#==============================================================================
# nlopt LIBRARY TARGET (SHARED OR STATIC)
#==============================================================================

configure_file (api/nlopt.h ${PROJECT_BINARY_DIR}/api/nlopt.h COPYONLY)

set (NLOPT_HEADERS
  ${PROJECT_BINARY_DIR}/api/nlopt.h ${PROJECT_BINARY_DIR}/api/nlopt.hpp ${PROJECT_BINARY_DIR}/api/nlopt.f
)

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
  slsqp/slsqp.c slsqp/slsqp.h
  esch/esch.c esch/esch.h
  api/general.c api/options.c api/optimize.c api/deprecated.c api/nlopt-internal.h api/nlopt.h api/f77api.c api/f77funcs.h api/f77funcs_.h api/nlopt.hpp api/nlopt-in.hpp
  util/mt19937ar.c util/sobolseq.c util/soboldata.h util/timer.c util/stop.c util/nlopt-util.h util/redblack.c util/redblack.h util/qsort_r.c util/rescale.c
)

if (WITH_CXX)
  list (APPEND NLOPT_SOURCES stogo/global.cc stogo/linalg.cc stogo/local.cc stogo/stogo.cc stogo/tools.cc stogo/global.h stogo/linalg.h stogo/local.h stogo/stogo_config.h stogo/stogo.h stogo/tools.h)
endif ()

install (FILES ${NLOPT_HEADERS} DESTINATION ${INSTALL_INCLUDE_DIR})

set (nlopt_lib nlopt${NLOPT_SUFFIX})
add_library (${nlopt_lib} ${NLOPT_SOURCES})
target_link_libraries (${nlopt_lib} ${M_LIBRARY})

set_target_properties (${nlopt_lib} PROPERTIES SOVERSION 0)
set_target_properties (${nlopt_lib} PROPERTIES VERSION 0.9.0)

#==============================================================================
# INCLUDE DIRECTORIES
#==============================================================================
target_include_directories (${nlopt_lib} PRIVATE
  ${PROJECT_BINARY_DIR}/api
  ${PROJECT_BINARY_DIR}
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
  api)

get_target_property (NLOPT_PRIVATE_INCLUDE_DIRS ${nlopt_lib} INCLUDE_DIRECTORIES)
target_include_directories (${nlopt_lib} INTERFACE "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/api>")

if (BUILD_SHARED_LIBS)
  target_compile_definitions (${nlopt_lib} INTERFACE -DNLOPT_DLL)
  target_compile_definitions (${nlopt_lib} INTERFACE -DNLOPT_DLL_EXPORT)
endif ()

# pass -fPIC in case swig module is built with static library
if (NOT BUILD_SHARED_LIBS)
  check_c_compiler_flag (-fPIC HAS_FPIC)
  if (HAS_FPIC)
    set (CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
  endif ()
endif ()

install ( TARGETS ${nlopt_lib}
          EXPORT  NLoptLibraryDepends
          RUNTIME DESTINATION ${RELATIVE_INSTALL_BIN_DIR}
          LIBRARY DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
          ARCHIVE DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
        )

add_subdirectory (api)

if (BUILD_PYTHON)
  find_package (PythonInterp)
  find_package (PythonLibs)
  find_package (NumPy)
endif ()

if (NOT DEFINED INSTALL_PYTHON_DIR AND PYTHONINTERP_FOUND)
  execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
                    OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

  get_filename_component (_ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE)
  file (RELATIVE_PATH _REL_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH})

  set (INSTALL_PYTHON_DIR ${_REL_PYTHON_MODULE_PATH})

endif ()

if (BUILD_GUILE)
  find_package (Guile)
endif ()

if (USE_SWIG)
  find_package (SWIG)
endif ()

add_subdirectory (swig)

if (BUILD_OCTAVE)
  find_package (Octave)
endif ()

if (BUILD_MATLAB)
  find_package (Matlab)
endif ()

if (OCTAVE_FOUND OR Matlab_FOUND)
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

# configuration files
export (TARGETS ${nlopt_lib} FILE ${PROJECT_BINARY_DIR}/NLoptLibraryDepends.cmake)

# Install the export set for use with the install-tree
install(EXPORT NLoptLibraryDepends
        DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR}
        COMPONENT Development)

# Create a NLOPTConfig.cmake file for the use from the install tree
# and install it
set (NLOPT_LIBRARIES ${nlopt_lib})

set (NLOPT_CMAKE_DIR "${INSTALL_CMAKE_DIR}")
file (RELATIVE_PATH rel_include_dir "${NLOPT_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
list (APPEND RELATIVE_NLOPT_INCLUDE_DIRS ${rel_include_dir})

file (RELATIVE_PATH rel_lib_dir "${NLOPT_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
list (APPEND RELATIVE_NLOPT_LIB_DIR ${rel_lib_dir})

configure_file (cmake/NLoptConfig.cmake.in NLoptConfig.cmake @ONLY)
configure_file (cmake/NLoptConfigVersion.cmake.in NLoptConfigVersion.cmake @ONLY)
install (FILES
          ${CMAKE_CURRENT_BINARY_DIR}/NLoptConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/NLoptConfigVersion.cmake
         DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR}
         COMPONENT Development)

