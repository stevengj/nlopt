# Copyright (c) 2008, 2014 OpenCog.org (http://opencog.org)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# - Try to find Guile; Once done this will define
#
# GUILE_FOUND            - system has the GUILE library
# GUILE_INCLUDE_DIRS     - the GUILE include directory
# GUILE_LIBRARIES        - The libraries needed to use GUILE
# GUILE_VERSION_STRING   - Version
# GUILE_SITE_DIR         - site dir
# GUILE_EXTENSION_DIR    - extension dir
# GUILE_ROOT_DIR         - prefix dir
find_package(PkgConfig REQUIRED)

pkg_check_modules(GUILE IMPORTED_TARGET guile)
if (GUILE_FOUND)
  pkg_get_variable(GUILE_ROOT_DIR guile prefix)
  pkg_get_variable(GUILE_SITE_DIR guile sitedir)
  pkg_get_variable(GUILE_EXTENSION_DIR guile extensiondir)
else()
  pkg_check_modules(GUILE IMPORTED_TARGET guile-2.0>=2.0)
endif()
if (GUILE_FOUND)
  pkg_get_variable(GUILE_ROOT_DIR guile-2.0 prefix)
  pkg_get_variable(GUILE_SITE_DIR guile-2.0 sitedir)
  pkg_get_variable(GUILE_EXTENSION_DIR guile-2.0 extensiondir)
else()
  pkg_check_modules(GUILE REQUIRED IMPORTED_TARGET guile-3.0>=3.0)
  pkg_get_variable(GUILE_ROOT_DIR guile-3.0 prefix)
  pkg_get_variable(GUILE_SITE_DIR guile-3.0 sitedir)
  pkg_get_variable(GUILE_EXTENSION_DIR guile-3.0 extensiondir)
endif()
message(STATUS "GUILE_VERSION is set to ${GUILE_VERSION}")
message(STATUS "GUILE_ROOT_DIR is set to ${GUILE_ROOT_DIR}")
message(STATUS "GUILE_SITE_DIR is set to ${GUILE_SITE_DIR}")
message(STATUS "GUILE_EXTENSION_DIR is set to ${GUILE_EXTENSION_DIR}")

# handle REQUIRED and QUIET options
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Guile REQUIRED_VARS GUILE_SITE_DIR GUILE_EXTENSION_DIR GUILE_ROOT_DIR GUILE_INCLUDE_DIRS GUILE_LIBRARIES GUILE_CFLAGS GUILE_LDFLAGS GUILE_VERSION)

mark_as_advanced (GUILE_INCLUDE_DIR GUILE_LIBRARY)

