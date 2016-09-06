# - Find the NumPy libraries
# This module finds if NumPy is installed, and sets the following variables
# indicating where it is.
#
# TODO: Update to provide the libraries and paths for linking npymath lib.
#
#  NUMPY_FOUND               - was NumPy found
#  NUMPY_VERSION_STRING      - the version of NumPy found as a string
#  NUMPY_VERSION_MAJOR       - the major version number of NumPy
#  NUMPY_VERSION_MINOR       - the minor version number of NumPy
#  NUMPY_VERSION_PATCH       - the patch version number of NumPy
#  NUMPY_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
#  NUMPY_INCLUDE_DIRS        - path to the NumPy include files

#============================================================================
# Copyright 2012 Continuum Analytics, Inc.
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
#============================================================================

# finding NumPy involves calling the Python interpreter
if (NOT DEFINED PYTHONINTERP_FOUND)
  if (NumPy_FIND_REQUIRED)
    find_package (PythonInterp REQUIRED)
  else ()
    find_package (PythonInterp)
  endif ()
endif ()

# find NumPy
if (PYTHONINTERP_FOUND)
  # get numpy version
  execute_process (
    COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "import numpy; print(numpy.__version__);"
    RESULT_VARIABLE _NUMPY_RESULT_CODE
    OUTPUT_VARIABLE _NUMPY_OUTPUT
    ERROR_VARIABLE _NUMPY_OUTPUT_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # check version
  unset (NUMPY_VERSION_STRING)
  if (_NUMPY_OUTPUT AND _NUMPY_RESULT_CODE MATCHES 0)
    # parse version
    string (REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)"
      _MATCH "${_NUMPY_OUTPUT}")
    if (_MATCH)
      set (NUMPY_VERSION_MAJOR ${CMAKE_MATCH_1})
      set (NUMPY_VERSION_MINOR ${CMAKE_MATCH_2})
      set (NUMPY_VERSION_PATCH ${CMAKE_MATCH_3})
      set (NUMPY_VERSION_STRING
        "${NUMPY_VERSION_MAJOR}.${NUMPY_VERSION_MINOR}.${NUMPY_VERSION_PATCH}")
      math (EXPR NUMPY_VERSION_DECIMAL "(${NUMPY_VERSION_MAJOR} * 10000) +
        (${NUMPY_VERSION_MINOR} * 100) + ${NUMPY_VERSION_PATCH}")
    endif ()
  elseif (NumPy_FIND_REQUIRED)
    message (SEND_ERROR "Could NOT get NumPy version: ${_NUMPY_OUTPUT_ERROR}")
  endif ()

  # get numpy include dirs
  execute_process (
    COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "import numpy; print(numpy.get_include());"
    RESULT_VARIABLE _NUMPY_RESULT_CODE
    OUTPUT_VARIABLE _NUMPY_OUTPUT
    ERROR_VARIABLE _NUMPY_OUTPUT_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # check headers
  unset (NUMPY_INCLUDE_DIR)
  if (_NUMPY_OUTPUT AND _NUMPY_RESULT_CODE MATCHES 0)
    # normalize path to cmake convention with '/' as separator
    file (TO_CMAKE_PATH "${_NUMPY_OUTPUT}" _NUMPY_OUTPUT)
    find_path (NUMPY_INCLUDE_DIR numpy/arrayobject.h
      HINTS
        "${_NUMPY_OUTPUT}"
        ${PYTHON_INCLUDE_DIRS}
    )
    set (NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIR})
  elseif (NumPy_FIND_REQUIRED)
    message (SEND_ERROR "Could NOT get NumPy includes: ${_NUMPY_OUTPUT_ERROR}")
  endif ()

  # cleanup
  unset (_NUMPY_RESULT_CODE)
  unset (_NUMPY_OUTPUT)
  unset (_NUMPY_OUTPUT_ERROR)
  unset (_MATCH)
endif ()

# handle REQUIRED and QUIET options
include (FindPackageHandleStandardArgs)
if (CMAKE_VERSION VERSION_LESS 2.8.3)
  find_package_handle_standard_args (NumPy
    DEFAULT_MSG
      NUMPY_INCLUDE_DIR
      NUMPY_VERSION_STRING
  )
else ()
  find_package_handle_standard_args (NumPy
    REQUIRED_VARS
      NUMPY_INCLUDE_DIR
      NUMPY_VERSION_STRING
    VERSION_VAR
      NUMPY_VERSION_STRING
  )
endif ()

mark_as_advanced (NUMPY_INCLUDE_DIR)
