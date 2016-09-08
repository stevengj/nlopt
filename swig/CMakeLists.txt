
#==============================================================================
# Python extension
#==============================================================================

if (NUMPY_FOUND AND PYTHONLIBS_FOUND AND (SWIG_FOUND OR
    (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt-python.cpp")))
  include_directories (${PYTHON_INCLUDE_DIRS})
  include_directories (${NUMPY_INCLUDE_DIRS})

  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt-python.cpp")
    set (SWIG_MODULE_nlopt_REAL_NAME _nlopt)
    add_library (${SWIG_MODULE_nlopt_REAL_NAME} MODULE nlopt-python.cpp)
    target_link_libraries (${SWIG_MODULE_nlopt_REAL_NAME} ${nlopt_lib})
    target_link_libraries (${SWIG_MODULE_nlopt_REAL_NAME} ${PYTHON_LIBRARIES})
    set_target_properties (${SWIG_MODULE_nlopt_REAL_NAME} PROPERTIES PREFIX "")
    if (NOT CMAKE_SOURCE_DIR MATCHES CMAKE_BINARY_DIR)
      execute_process (COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/nlopt.py"
        "${CMAKE_CURRENT_BINARY_DIR}/nlopt.py")
    endif ()
  else ()
    include (UseSWIG)
    set_source_files_properties (nlopt.i PROPERTIES CPLUSPLUS ON)
    set (SWIG_MODULE_nlopt_EXTRA_DEPS nlopt-python.i numpy.i)
    swig_add_module (nlopt python nlopt.i)
    swig_link_libraries (nlopt ${nlopt_lib})
    swig_link_libraries (nlopt ${PYTHON_LIBRARIES})
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
    "${CMAKE_CURRENT_BINARY_DIR}/nloptGUILE_wrap.cxx")
  target_link_libraries(nlopt_guile ${nlopt_lib} ${GUILE_LIBRARIES})
  set_target_properties (nlopt_guile PROPERTIES NO_SONAME ON)

  install (
    TARGETS nlopt_guile
    LIBRARY DESTINATION "${RELATIVE_GUILE_EXTENSION_DIR}"
    COMPONENT "Bindings")
  install (
    FILES "${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm"
    DESTINATION "${RELATIVE_GUILE_SITE_DIR}"
    COMPONENT "Bindings")
endif ()
