
include (UseSWIG)

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-enum-renames.i)
  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i "// AUTOMATICALLY GENERATED -- DO NOT EDIT\n")
  file (STRINGS ${PROJECT_SOURCE_DIR}/src/api/nlopt.h NLOPT_H_LINES REGEX "    NLOPT_[A-Z0-9_]+")
  foreach (NLOPT_H_LINE ${NLOPT_H_LINES})
    string (REGEX REPLACE ".*NLOPT_([A-Z0-9_]+).*" "%rename(NLOPT_\\1) nlopt::\\1;\n" ENUM_LINE ${NLOPT_H_LINE})
    file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i "${ENUM_LINE}")
  endforeach ()
endif ()

include_directories (${NLOPT_PRIVATE_INCLUDE_DIRS})
set_source_files_properties (nlopt.i PROPERTIES CPLUSPLUS ON)

if (NUMPY_FOUND AND PYTHONLIBS_FOUND)

  set (SWIG_MODULE_nlopt_EXTRA_DEPS nlopt-python.i numpy.i)

  if (CMAKE_VERSION VERSION_LESS 3.8)
    swig_add_module (nlopt python nlopt.i)
  else ()
    swig_add_library (nlopt LANGUAGE python SOURCES nlopt.i)
  endif ()

  target_include_directories (${SWIG_MODULE_nlopt_REAL_NAME} PUBLIC ${PYTHON_INCLUDE_DIRS})
  target_include_directories (${SWIG_MODULE_nlopt_REAL_NAME} PUBLIC ${NUMPY_INCLUDE_DIRS})

  swig_link_libraries (nlopt ${nlopt_lib})
  if (NLOPT_LINK_PYTHON)
    swig_link_libraries (nlopt ${PYTHON_LIBRARIES})
  endif ()

  set_target_properties (${SWIG_MODULE_nlopt_REAL_NAME} PROPERTIES OUTPUT_NAME _nlopt)

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.py DESTINATION ${INSTALL_PYTHON_DIR})
  install (TARGETS ${SWIG_MODULE_nlopt_REAL_NAME} LIBRARY DESTINATION ${INSTALL_PYTHON_DIR})
endif ()


# guile bindings with gcc only ok with swig >= 2.0.10
if (GUILE_FOUND AND (SWIG_VERSION VERSION_GREATER 2.0.9))

  set_source_files_properties (nlopt.i PROPERTIES SWIG_FLAGS "-scmstub")
  set (SWIG_MODULE_nlopt_guile_EXTRA_DEPS nlopt-guile.i)

  if (CMAKE_VERSION VERSION_LESS 3.8)
    swig_add_module (nlopt_guile guile nlopt.i)
  else ()
    swig_add_library (nlopt_guile LANGUAGE guile SOURCES nlopt.i)
  endif ()

  target_include_directories (${SWIG_MODULE_nlopt_guile_REAL_NAME} PUBLIC ${GUILE_INCLUDE_DIRS})

  swig_link_libraries (nlopt_guile ${nlopt_lib})
  swig_link_libraries (nlopt_guile ${GUILE_LIBRARIES})

  file (RELATIVE_PATH _REL_GUILE_SITE_PATH ${GUILE_ROOT_DIR} ${GUILE_SITE_DIR})
  set (GUILE_SITE_PATH ${_REL_GUILE_SITE_PATH})
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm DESTINATION ${GUILE_SITE_PATH})

  file (RELATIVE_PATH _REL_GUILE_EXTENSION_PATH ${GUILE_ROOT_DIR} ${GUILE_EXTENSION_DIR})
  set (GUILE_EXTENSION_PATH ${_REL_GUILE_EXTENSION_PATH})
  install (TARGETS ${SWIG_MODULE_nlopt_guile_REAL_NAME} LIBRARY DESTINATION ${GUILE_EXTENSION_PATH})
endif ()
