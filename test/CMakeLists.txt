
#==============================================================================
# Test Executables
#==============================================================================

# have to add timer.c and mt19937ar.c as symbols are declared extern
set (TESTOPT_SOURCES
  testopt.c testfuncs.c testfuncs.h
  "${CMAKE_SOURCE_DIR}/util/timer.c"
  "${CMAKE_SOURCE_DIR}/util/mt19937ar.c"
)
if (NOT HAVE_GETOPT_H)
  list (APPEND TESTOPT_SOURCES "${CMAKE_SOURCE_DIR}/util/nlopt-getopt.c")
endif ()
add_executable (testopt EXCLUDE_FROM_ALL ${TESTOPT_SOURCES})
target_link_libraries (testopt ${NLOPT_LIB})

if (WITH_FOLDERS)
  set_target_properties (testopt PROPERTIES FOLDER "Testing")
endif ()

#==============================================================================
# Convenience targets
#==============================================================================

# group all tests executables under one target
add_custom_target (tests)
add_dependencies (tests testopt)

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
foreach (_algo_index RANGE 29)# 42
  foreach (_obj_index RANGE 1)# 21
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

#==============================================================================
# Python tests
#==============================================================================

if (PYTHONINTERP_FOUND AND NUMPY_FOUND AND SWIG_FOUND)
  set (PYINSTALLCHECK_ENVIRONMENT
    "LD_LIBRARY_PATH=${INSTALL_DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}:$ENV{LD_LIBRARY_PATH}"
    "PYTHONPATH=${INSTALL_DESTDIR}${PYTHON_MODULE_PATH}${PATH_SEP}$ENV{PYTHONPATH}"
  )
  add_test (NAME test_python
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_python.py")
  set_tests_properties (test_python
    PROPERTIES ENVIRONMENT "${PYINSTALLCHECK_ENVIRONMENT}")
endif ()

#==============================================================================
# Guile tests
#==============================================================================

if (GUILE_FOUND AND SWIG_FOUND AND SWIG_VERSION VERSION_GREATER 2.0.9)
  set (GUILECHECK_ENVIRONMENT
    "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/swig"
    "GUILE_LOAD_PATH=${CMAKE_BINARY_DIR}/swig"
  )
  add_test (NAME test_guile
    COMMAND "${GUILE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_guile.scm")
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
