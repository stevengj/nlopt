
#==============================================================================
# Helpers
#==============================================================================

set (BUILD_TESTS_UTIL ON)

# add test executable linked against nlopt
macro (add_executable_nlopt _NAME)
  add_executable (${_NAME} EXCLUDE_FROM_ALL ${ARGN})
  target_link_libraries (${_NAME} ${NLOPT_LIB})
  if (WITH_FOLDERS)
    set_target_properties (${_NAME} PROPERTIES FOLDER "Testing")
  endif ()
endmacro ()

#==============================================================================
# Test Executables
#==============================================================================

# nlopt tests
# (note: have to add some util/*.c as symbols are declared extern)
set (TESTOPT_SOURCES
  testopt.c testfuncs.c testfuncs.h
  ../util/timer.c ../util/mt19937ar.c)
if (NOT HAVE_GETOPT_H)
  list (APPEND TESTOPT_SOURCES ../util/nlopt-getopt.c)
endif ()
add_executable_nlopt (testopt ${TESTOPT_SOURCES})
add_executable_nlopt (lorentzfit lorentzfit.c ../util/mt19937ar.c)
add_executable_nlopt (box box.c)

# util tests
if (BUILD_TESTS_UTIL)
  add_executable (mt19937ar_test EXCLUDE_FROM_ALL
    ../util/mt19937ar_test.c ../util/mt19937ar.c)
  add_executable (redblack_test EXCLUDE_FROM_ALL
    ../util/redblack_test.c ../util/redblack.c)
  add_executable (sobolseq_test EXCLUDE_FROM_ALL
    ../util/sobolseq_test.c ../util/sobolseq.c ../util/mt19937ar.c)

  if (UNIX OR HAVE_LIBM)
    target_link_libraries (mt19937ar_test m)
    target_link_libraries (redblack_test m)
    target_link_libraries (sobolseq_test m)
  endif ()
endif ()

#==============================================================================
# Convenience targets
#==============================================================================

# group all tests executables under one target
add_custom_target (tests)
add_dependencies (tests testopt lorentzfit box)
if (BUILD_TESTS_UTIL)
  add_dependencies (tests mt19937ar_test redblack_test sobolseq_test)
endif ()

# emulate Autotools "make check" to build and run tests
add_custom_target (check
  COMMAND ${CMAKE_CTEST_COMMAND} -C ${CMAKE_BUILD_TYPE}
  DEPENDS tests)

if (WITH_FOLDERS)
  set_target_properties (tests check PROPERTIES FOLDER "Testing")
endif ()

#==============================================================================
# C/C++ tests
#==============================================================================

# testopt (sweep over various algorithms and objective functions)
foreach (_algo_index RANGE 29)  # 42
  foreach (_obj_index RANGE 1)  # 21
    set (_enable TRUE)
    # cxx stogo
    if (NOT WITH_CXX)
      if (_algo_index STREQUAL 8 OR _algo_index STREQUAL 9)
        set (_enable FALSE)
      endif ()
    endif ()
    # L-BFGS
    if (_algo_index STREQUAL 10)
      set (_enable FALSE)
    endif ()
    # LN_NEWUOA only works for dimensionality >= 2
    if (_algo_index STREQUAL 26 AND _obj_index STREQUAL 19)
      set (_enable FALSE)
    endif ()
    # LN_NEWUOA_BOUND timeouts on msvc
    if (_algo_index STREQUAL 27)
      set (_enable FALSE)
    endif ()
    if (_enable)
      add_test (NAME "testopt_algo${_algo_index}_obj${_obj_index}"
        COMMAND testopt -a ${_algo_index} -o ${_obj_index})
    endif ()
  endforeach ()
endforeach ()

# more nlopt tests
add_test (NAME "lorentzfit" COMMAND lorentzfit)
add_test (NAME "box" COMMAND box)

# util tests
if (BUILD_TESTS_UTIL)
  add_test (NAME "mt19937ar_test" COMMAND mt19937ar_test)
  add_test (NAME "redblack_test" COMMAND redblack_test 20)
  add_test (NAME "sobolseq_test" COMMAND sobolseq_test 2 100)
endif ()

#==============================================================================
# Python tests
#==============================================================================

if (PYTHONINTERP_FOUND AND NUMPY_FOUND AND SWIG_FOUND)
  add_test (NAME test_python
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_python.py")

  # setup test environment
  set (PYINSTALLCHECK_ENVIRONMENT
    "LD_LIBRARY_PATH=${INSTALL_DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}:$ENV{LD_LIBRARY_PATH}"
    "PYTHONPATH=${INSTALL_DESTDIR}${PYTHON_MODULE_PATH}${PATH_SEP}$ENV{PYTHONPATH}"
  )
  set_tests_properties (test_python
    PROPERTIES ENVIRONMENT "${PYINSTALLCHECK_ENVIRONMENT}")
endif ()

#==============================================================================
# Guile tests
#==============================================================================

if (GUILE_FOUND AND SWIG_FOUND AND SWIG_VERSION VERSION_GREATER 2.0.9)
  add_test (NAME test_guile
    COMMAND "${GUILE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_guile.scm")

  # setup test environment
  set (GUILECHECK_ENVIRONMENT
    "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/swig"
    "GUILE_LOAD_PATH=${CMAKE_BINARY_DIR}/swig"
  )
  set_tests_properties (test_guile
    PROPERTIES ENVIRONMENT "${GUILECHECK_ENVIRONMENT}")
endif ()

#==============================================================================
# Octave tests
#==============================================================================

if (OCTAVE_FOUND)
  add_test (NAME test_octave
    COMMAND "${OCTAVE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_octave.m"
      ${CMAKE_SOURCE_DIR}/octave ${CMAKE_BINARY_DIR}/octave)
endif ()
