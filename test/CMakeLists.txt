
#==============================================================================
# Helper functions and macros
#==============================================================================

# get path separator, LD path environmnent variable name, and LD install dir
# in a cross-platform way
if (WIN32)
  set (PATH_SEP "\;")
  set (PATH_VAR "PATH")
  set (PATH_VAL "${CMAKE_INSTALL_FULL_BINDIR}")
elseif (APPLE)
  set (PATH_SEP ":")
  set (PATH_VAR "DYLD_LIBRARY_PATH")
  set (PATH_VAL "${CMAKE_INSTALL_FULL_LIBDIR}")
else ()
  set (PATH_SEP ":")
  set (PATH_VAR "LD_LIBRARY_PATH")
  set (PATH_VAL "${CMAKE_INSTALL_FULL_LIBDIR}")
endif ()

# get value of an environment variable properly escaped if needed
macro (get_env _VARIABLE _NAME)
  set (${_VARIABLE} $ENV{${_NAME}})
  if (WIN32)
    string (REPLACE ";" "\\;" ${_VARIABLE} "${${_VARIABLE}}")
  endif ()
endmacro ()

# prepend a common path to a list of paths
function (prepend VAR PREFIX)
   set (LISTVAR "")
   foreach (P ${ARGN})
      list (APPEND LISTVAR "${PREFIX}/${P}")
   endforeach (P)
   set (${VAR} "${LISTVAR}" PARENT_SCOPE)
endfunction ()

# add test executable linked against nlopt
function (add_executable_nlopt NAME)
  add_executable (${NAME} EXCLUDE_FROM_ALL ${ARGN})
  target_link_libraries (${NAME} ${NLOPT_LIB})
  if (WITH_FOLDERS)
    set_target_properties (${NAME} PROPERTIES FOLDER "Testing")
  endif ()
endfunction ()

# add test executable for util functions
function (add_executable_util NAME)
  add_executable (${NAME} EXCLUDE_FROM_ALL ${ARGN})
  if (UNIX OR HAVE_LIBM)
    target_link_libraries (${NAME} m)
  endif ()
  if (WITH_FOLDERS)
    set_target_properties (${NAME} PROPERTIES FOLDER "Testing")
  endif ()
endfunction ()

# add test executable for an internal algorithm
function (add_executable_alg NAME ALG)
  prepend (NLOPT_SOURCES_${ALG}2 .. ${NLOPT_SOURCES_${ALG}})
  prepend (NLOPT_SOURCES_UTIL2 .. ${NLOPT_SOURCES_UTIL})
  add_executable (${NAME} EXCLUDE_FROM_ALL ${ARGN}
    ${NLOPT_SOURCES_${ALG}2} ${NLOPT_SOURCES_UTIL2})
  if (UNIX OR HAVE_LIBM)
    target_link_libraries (${NAME} m)
  endif ()
  if (WITH_FOLDERS)
    set_target_properties (${NAME} PROPERTIES FOLDER "Testing")
  endif ()
endfunction ()

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
  add_executable_util (mt19937ar_test ../util/mt19937ar_test.c
    ../util/mt19937ar.c)
  add_executable_util (redblack_test ../util/redblack_test.c
    ../util/redblack.c)
  add_executable_util (sobolseq_test ../util/sobolseq_test.c
    ../util/sobolseq.c ../util/mt19937ar.c)
endif ()

# algorithm tests
if (BUILD_TESTS_ALG)
  add_executable_alg (direct_tstc DIRECT ../direct/tstc.c)
  if (WITH_CXX)
    add_executable_alg (stogo_prog STOGO ../stogo/prog.cc)
    add_executable_alg (stogo_testros STOGO ../stogo/testros.cc)
    add_executable_alg (stogo_tst STOGO ../stogo/tst.cc)
    add_executable_alg (stogo_tstc STOGO ../stogo/tstc.c)
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
if (BUILD_TESTS_ALG)
  add_dependencies (tests direct_tstc)
  if (WITH_CXX)
    add_dependencies (tests stogo_prog stogo_testros stogo_tst stogo_tstc)
  endif ()
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

# collect all C tests
set (C_TESTS "")

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
      list (APPEND C_TESTS "testopt_algo${_algo_index}_obj${_obj_index}")
    endif ()
  endforeach ()
endforeach ()

# more nlopt tests
add_test (NAME "lorentzfit" COMMAND lorentzfit)
add_test (NAME "box" COMMAND box)
list (APPEND C_TESTS lorentzfit box)

# util tests
if (BUILD_TESTS_UTIL)
  add_test (NAME "mt19937ar_test" COMMAND mt19937ar_test)
  add_test (NAME "redblack_test" COMMAND redblack_test 20)
  add_test (NAME "sobolseq_test" COMMAND sobolseq_test 2 100)
  list (APPEND C_TESTS mt19937ar_test redblack_test sobolseq_test)
endif ()

# internal tests
if (BUILD_TESTS_ALG)
  add_test (NAME "direct_tstc" COMMAND direct_tstc)
  list (APPEND C_TESTS direct_tstc)
  if (WITH_CXX)
    if (UNIX)
      # conditional, otherwise we get "echo.exe not found" error on Windows
      # (this is an interactive program that prompts for input from stdin)
      add_test (NAME "stogo_prog" COMMAND echo 1 1 y | stogo_prog)
      list (APPEND C_TESTS stogo_prog)
    endif ()
    add_test (NAME "stogo_testros" COMMAND stogo_testros)
    add_test (NAME "stogo_tst" COMMAND stogo_tst)
    add_test (NAME "stogo_tstc" COMMAND stogo_tstc)
    list (APPEND C_TESTS stogo_testros stogo_tst stogo_tstc)
  endif ()
endif ()

# setup test environment
get_env (PATH_ENV ${PATH_VAR})
if (WITH_TESTS_INSTALL_TREE)
  set (C_TEST_ENVIRONMENT
    "${PATH_VAR}=${PATH_VAL}${PATH_SEP}${PATH_ENV}")
else ()
  # NOTE: inside build tree, Visual Studio places binaries in a subdirectory
  # named after the config type (Release, Debug, ...)
  set (C_TEST_ENVIRONMENT
    "${PATH_VAR}=${CMAKE_BINARY_DIR}${PATH_SEP}${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}${PATH_SEP}${PATH_ENV}")
endif ()
set_tests_properties (${C_TESTS}
  PROPERTIES ENVIRONMENT "${C_TEST_ENVIRONMENT}")

#TODO: when setting tests env, install dirs are still in cmake format
# which can be problematic on windows, ie: PATH=C:/build/install;C:\Windows;...

#==============================================================================
# Python tests
#==============================================================================

if (PYTHONINTERP_FOUND AND NUMPY_FOUND AND SWIG_FOUND)
  add_test (NAME test_python
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_python.py")

  # setup test environment
  if (WITH_TESTS_INSTALL_TREE)
    set (PYTHON_TEST_ENVIRONMENT
      "${PATH_VAR}=${PATH_VAL}${PATH_SEP}${PATH_ENV}"
      "PYTHONPATH=${INSTALL_PYTHON_FULL_LIB_DIR}${PATH_SEP}$ENV{PYTHONPATH}")
  else ()
    set (PYTHON_TEST_ENVIRONMENT
      "${PATH_VAR}=${CMAKE_BINARY_DIR}${PATH_SEP}${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}${PATH_SEP}${PATH_ENV}"
      "PYTHONPATH=${CMAKE_BINARY_DIR}/swig${PATH_SEP}${CMAKE_BINARY_DIR}/swig/${CMAKE_BUILD_TYPE}${PATH_SEP}$ENV{PYTHONPATH}")
  endif ()
  set_tests_properties (test_python
    PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")
endif ()

#==============================================================================
# Guile tests
#==============================================================================

if (GUILE_FOUND AND SWIG_FOUND AND SWIG_VERSION VERSION_GREATER 2.0.9)
  add_test (NAME test_guile
    COMMAND "${GUILE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/t_guile.scm")

  # setup test environment
  if (WITH_TESTS_INSTALL_TREE)
    set (GUILE_TEST_ENVIRONMENT
      "${PATH_VAR}=${INSTALL_GUILE_FULL_EXT_DIR}${PATH_SEP}${PATH_VAL}${PATH_SEP}${PATH_ENV}"
      "GUILE_LOAD_PATH=${INSTALL_GUILE_FULL_SITE_DIR}${PATH_SEP}...")
  else ()
    set (GUILE_TEST_ENVIRONMENT
      "${PATH_VAR}=${CMAKE_BINARY_DIR}/swig${PATH_SEP}${CMAKE_BINARY_DIR}/swig/${CMAKE_BUILD_TYPE}${PATH_SEP}${CMAKE_BINARY_DIR}${PATH_SEP}${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}${PATH_SEP}${PATH_ENV}"
      "GUILE_LOAD_PATH=${CMAKE_BINARY_DIR}/swig${PATH_SEP}${CMAKE_BINARY_DIR}/swig/${CMAKE_BUILD_TYPE}${PATH_SEP}...")
  endif ()
  set_tests_properties (test_guile
    PROPERTIES ENVIRONMENT "${GUILE_TEST_ENVIRONMENT}")
endif ()

#==============================================================================
# Octave tests
#==============================================================================

if (OCTAVE_FOUND)
  # octave executable options
  set (OCTAVE_OPTS --quiet --norc)
  if (NOT OCTAVE_VERSION_STRING VERSION_LESS 4.0.0)
    list (APPEND OCTAVE_OPTS --no-gui --no-window-system)
  endif ()

  add_test (NAME test_octave
    COMMAND "${OCTAVE_EXECUTABLE}" ${OCTAVE_OPTS}
      --eval "cd('${CMAKE_CURRENT_SOURCE_DIR}'), t_octave, exit")

  # setup test environment
  if (WITH_TESTS_INSTALL_TREE)
    set (OCTAVE_TEST_ENVIRONMENT
      "${PATH_VAR}=${PATH_VAL}${PATH_SEP}${PATH_ENV}"
      "OCTAVE_PATH=${INSTALL_OCTAVE_FULL_OCT_DIR}${PATH_SEP}${INSTALL_OCTAVE_FULL_M_DIR}${PATH_SEP}")
  else ()
    set (OCTAVE_TEST_ENVIRONMENT
      "${PATH_VAR}=${CMAKE_BINARY_DIR}${PATH_SEP}${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}${PATH_SEP}${PATH_ENV}"
      "OCTAVE_PATH=${CMAKE_BINARY_DIR}/octave${PATH_SEP}${CMAKE_BINARY_DIR}/octave/${CMAKE_BUILD_TYPE}${PATH_SEP}${CMAKE_SOURCE_DIR}/octave${PATH_SEP}")
  endif ()
  set_tests_properties (test_octave
    PROPERTIES ENVIRONMENT "${OCTAVE_TEST_ENVIRONMENT}")
endif ()

#==============================================================================
# MATLAB tests
#==============================================================================

if (Matlab_FOUND)
  # MATLAB executable options
  set (MATLAB_OPTS -nodesktop -nodisplay -nosplash)
  if (WIN32)
    list (APPEND MATLAB_OPTS -wait)  #TODO: unix alternative ?
  endif ()

  # unfortunately MATLABPATH env var only works in Linux/Mac
  # so we will manually call addpath() instead
  if (WITH_TESTS_INSTALL_TREE)
    set (ADDPATH_CALL "addpath('${INSTALL_MATLAB_FULL_DIR}')")
  else ()
    set (ADDPATH_CALL "addpath('${CMAKE_SOURCE_DIR}/octave'), addpath('${CMAKE_BINARY_DIR}/octave'), addpath('${CMAKE_BINARY_DIR}/octave/${CMAKE_BUILD_TYPE}')")
  endif ()

  add_test (NAME test_matlab
    COMMAND "${Matlab_MAIN_PROGRAM}" ${MATLAB_OPTS}
      -r "try, ${ADDPATH_CALL}, end, cd('${CMAKE_CURRENT_SOURCE_DIR}'), t_octave, exit")

  # setup test environment
  if (WITH_TESTS_INSTALL_TREE)
    set (MATLAB_TEST_ENVIRONMENT
      "${PATH_VAR}=${PATH_VAL}${PATH_SEP}${PATH_ENV}")
  else ()
    set (MATLAB_TEST_ENVIRONMENT
      "${PATH_VAR}=${CMAKE_BINARY_DIR}${PATH_SEP}${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}${PATH_SEP}${PATH_ENV}")
  endif ()
  set_tests_properties (test_matlab
    PROPERTIES ENVIRONMENT "${MATLAB_TEST_ENVIRONMENT}")
endif ()
