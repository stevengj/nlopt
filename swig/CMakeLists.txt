
#==============================================================================
# Python extension
#==============================================================================

if (PYTHONLIBS_FOUND AND NUMPY_FOUND AND SWIG_FOUND)
  # UseSWIG reads INCLUDE_DIRECTORIES to retrieve include dirs
  include_directories (${PYTHON_INCLUDE_DIRS})
  include_directories (${NUMPY_INCLUDE_DIRS})

  # use swig to generate python wrapper and Python lib target
  # (NOTE: no conflict in target name, real name in swig macros is _nlopt)
  include (UseSWIG)
  set_source_files_properties (nlopt.i PROPERTIES
    SWIG_MODULE_NAME nlopt
    CPLUSPLUS ON)
  set (SWIG_MODULE_nlopt_EXTRA_DEPS
    nlopt-exceptions.i nlopt-python.i numpy.i
    "${CMAKE_BINARY_DIR}/api/nlopt.hpp"
    "${CMAKE_SOURCE_DIR}/api/nlopt.h")
  swig_add_module (nlopt python nlopt.i)
  swig_link_libraries (nlopt ${NLOPT_LIB} ${PYTHON_LIBRARIES})
  set_target_properties (${SWIG_MODULE_nlopt_REAL_NAME} PROPERTIES
    NO_SONAME ON)

  if (WITH_FOLDERS)
    set_target_properties (${SWIG_MODULE_nlopt_REAL_NAME} PROPERTIES
      FOLDER "Python")
  endif ()

  install (
    TARGETS ${SWIG_MODULE_nlopt_REAL_NAME}
    LIBRARY DESTINATION "${RELATIVE_PYTHON_MODULE_PATH}"
    COMPONENT "Bindings")
  install (
    FILES "${CMAKE_CURRENT_BINARY_DIR}/nlopt.py"
    DESTINATION "${RELATIVE_PYTHON_MODULE_PATH}"
    COMPONENT "Bindings")
endif ()

#==============================================================================
# Guile extension
#==============================================================================

# generate nlopt-enum-renames.i from nlopt.h
if (BUILD_GUILE AND "${CMAKE_SOURCE_DIR}/api/nlopt.h" IS_NEWER_THAN
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i")
  file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i"
    "// AUTOMATICALLY GENERATED -- DO NOT EDIT\n")
  file (STRINGS "${CMAKE_SOURCE_DIR}/api/nlopt.h" NLOPT_H_LINES
    REGEX "    NLOPT_[A-Z0-9_]+")
  foreach (NLOPT_H_LINE ${NLOPT_H_LINES})
    string (REGEX REPLACE ".*NLOPT_([A-Z0-9_]+).*"
      "%rename(NLOPT_\\1) nlopt::\\1;\n" ENUM_LINE ${NLOPT_H_LINE})
    file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i"
      "${ENUM_LINE}")
  endforeach ()

  # cleanup
  unset (ENUM_LINE)
  unset (NLOPT_H_LINE)
  unset (NLOPT_H_LINES)

  # mark file as generated
  set_source_files_properties (
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i"
    PROPERTIES GENERATED TRUE)
endif ()

# guile bindings with gcc only ok with swig >= 2.0.10
if (GUILE_FOUND AND SWIG_FOUND AND SWIG_VERSION VERSION_GREATER 2.0.9)
  include_directories (${GUILE_INCLUDE_DIRS})

  # generate guile wrappers using swig
  set_source_files_properties (nlopt.i PROPERTIES CPLUSPLUS ON)
  add_custom_command (
    OUTPUT
      nloptGUILE_wrap.cxx
      nlopt.scm
    COMMAND "${SWIG_EXECUTABLE}"
    ARGS
      -guile -scmstub -c++
      -I${CMAKE_BINARY_DIR}/api
      -module nlopt
      -outdir "${CMAKE_CURRENT_BINARY_DIR}"
      -o nloptGUILE_wrap.cxx
      "${CMAKE_CURRENT_SOURCE_DIR}/nlopt.i"
    MAIN_DEPENDENCY nlopt.i
    DEPENDS
      nlopt-exceptions.i nlopt-guile.i
      "${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i"
      "${CMAKE_BINARY_DIR}/api/nlopt.hpp"
      "${CMAKE_SOURCE_DIR}/api/nlopt.h"
    COMMENT "Swig source")
  set_source_files_properties (
    "${CMAKE_CURRENT_BINARY_DIR}/nloptGUILE_wrap.cxx"
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm"
    PROPERTIES GENERATED TRUE)

  # Guile lib target
  add_library (nlopt_guile MODULE
    "${CMAKE_CURRENT_BINARY_DIR}/nloptGUILE_wrap.cxx"
    "${CMAKE_BINARY_DIR}/api/nlopt.hpp"
    "${CMAKE_SOURCE_DIR}/api/nlopt.h")
  target_link_libraries (nlopt_guile ${NLOPT_LIB} ${GUILE_LIBRARIES})
  set_target_properties (nlopt_guile PROPERTIES
    NO_SONAME ON)

  if (WITH_FOLDERS)
    set_target_properties (nlopt_guile PROPERTIES FOLDER "Guile")
  endif ()

  install (
    TARGETS nlopt_guile
    LIBRARY DESTINATION "${RELATIVE_GUILE_EXTENSION_DIR}"
    COMPONENT "Bindings")
  install (
    FILES "${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm"
    DESTINATION "${RELATIVE_GUILE_SITE_DIR}"
    COMPONENT "Bindings")
endif ()
