# - Find Guile
# https://www.gnu.org/software/guile/
#
# This module defines:
#
#  GUILE_FOUND             - system has the GUILE library
#  GUILE_INCLUDE_DIRS      - the GUILE include directory
#  GUILE_LIBRARIES         - The libraries needed to use GUILE
#  GUILE_VERSION_STRING    - Version
#  GUILE_ROOT_DIR          - prefix dir
#  GUILE_SITE_DIR          - site dir
#  GUILE_EXTENSION_DIR     - extension dir

#=============================================================================
# Copyright (c) 2008, 2014 OpenCog.org (http://opencog.org)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# find guile executables
find_program (GUILE_EXECUTABLE
  NAMES guile
)
find_program (GUILE_CONFIG_EXECUTABLE
  NAMES guile-config
)

if (GUILE_CONFIG_EXECUTABLE)
  # query guile directories
  execute_process (
    COMMAND ${GUILE_CONFIG_EXECUTABLE} info prefix
    OUTPUT_VARIABLE GUILE_ROOT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process (
    COMMAND ${GUILE_CONFIG_EXECUTABLE} info sitedir
    OUTPUT_VARIABLE GUILE_SITE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process (
    COMMAND ${GUILE_CONFIG_EXECUTABLE} info extensiondir
    OUTPUT_VARIABLE GUILE_EXTENSION_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif ()

# Look for the header file
# Look for guile-2.2 first, then 2.0, then 1.8
# Macports for OSX puts things in /opt/local
find_path (GUILE_INCLUDE_DIR libguile.h
  PATH_SUFFIXES
    guile/2.2
    guile/2.0
    guile/1.8
    libguile
    guile
  HINTS
    /opt/local/include
)
set (GUILE_INCLUDE_DIRS ${GUILE_INCLUDE_DIR})

# Look for the library
find_library (GUILE_LIBRARY
  NAMES guile-2.2 guile-2.0 guile
  HINTS
    /opt/local/lib
)
set (GUILE_LIBRARIES ${GUILE_LIBRARY})

# check guile's version if we're using cmake >= 2.6
if (GUILE_INCLUDE_DIR)
  SET(GUILE_VERSION_MAJOR 0)
  SET(GUILE_VERSION_MINOR 0)
  SET(GUILE_VERSION_PATCH 0)

  IF(NOT EXISTS "${GUILE_INCLUDE_DIR}/libguile/version.h")
          MESSAGE(FATAL_ERROR "Found ${GUILE_INCLUDE_DIR}/libguile.h but not version.h; check your guile installation!")
  ENDIF(NOT EXISTS "${GUILE_INCLUDE_DIR}/libguile/version.h")

  # Extract the libguile version from the 'version.h' file
  SET(GUILE_MAJOR_VERSION 0)
  FILE(READ "${GUILE_INCLUDE_DIR}/libguile/version.h" _GUILE_VERSION_H_CONTENTS)

  STRING(REGEX MATCH "#define SCM_MAJOR_VERSION[	 ]+([0-9])" _MATCH "${_GUILE_VERSION_H_CONTENTS}")
  SET(GUILE_VERSION_MAJOR ${CMAKE_MATCH_1})
  STRING(REGEX MATCH "#define SCM_MINOR_VERSION[	 ]+([0-9]+)" _MATCH "${_GUILE_VERSION_H_CONTENTS}")
  SET(GUILE_VERSION_MINOR ${CMAKE_MATCH_1})
  STRING(REGEX MATCH "#define SCM_MICRO_VERSION[	 ]+([0-9]+)" _MATCH "${_GUILE_VERSION_H_CONTENTS}")
  SET(GUILE_VERSION_PATCH ${CMAKE_MATCH_1})

  SET(GUILE_VERSION_STRING "${GUILE_VERSION_MAJOR}.${GUILE_VERSION_MINOR}.${GUILE_VERSION_PATCH}")

endif ()

#if (GUILE_FOUND AND GUILE_VERSION_MAJOR EQUAL 2)
#  add_definitions (-DHAVE_GUILE2)
#endif ()

# handle REQUIRED and QUIET options
include (FindPackageHandleStandardArgs)
if (CMAKE_VERSION VERSION_LESS 2.8.3)
  find_package_handle_standard_args (Guile
    DEFAULT_MSG
      GUILE_ROOT_DIR
      GUILE_INCLUDE_DIRS
      GUILE_LIBRARIES
      GUILE_VERSION_STRING)
else ()
  find_package_handle_standard_args (Guile
    REQUIRED_VARS
      GUILE_ROOT_DIR
      GUILE_INCLUDE_DIRS
      GUILE_LIBRARIES
    VERSION_VAR
      GUILE_VERSION_STRING)
endif ()

mark_as_advanced (
  GUILE_INCLUDE_DIR
  GUILE_LIBRARY
)
