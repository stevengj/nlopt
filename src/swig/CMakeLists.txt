
if (POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif ()
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

if (NUMPY_FOUND AND Python_FOUND)

  set (SWIG_MODULE_nlopt_python_EXTRA_DEPS nlopt-python.i numpy.i generate-cpp)

  # swig_add_module is deprecated
  swig_add_library (nlopt_python LANGUAGE python SOURCES nlopt.i)

  target_include_directories (nlopt_python PRIVATE ${Python_INCLUDE_DIRS})
  target_include_directories (nlopt_python PRIVATE ${NUMPY_INCLUDE_DIRS})

  swig_link_libraries (nlopt_python ${nlopt_lib})
  if (TARGET Python::Module)
    target_link_libraries (nlopt_python Python::Module)
  else ()
    include (TargetLinkLibrariesWithDynamicLookup)
    target_link_libraries_with_dynamic_lookup (${nlopt_python} ${Python_LIBRARIES})
  endif ()

  set_target_properties (nlopt_python PROPERTIES OUTPUT_NAME nlopt)

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.py DESTINATION ${INSTALL_PYTHON_DIR})
  install (TARGETS nlopt_python DESTINATION ${INSTALL_PYTHON_DIR})

  configure_file (METADATA.in METADATA @ONLY)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/METADATA
           DESTINATION ${INSTALL_PYTHON_DIR}/nlopt-${NLOPT_VERSION_STRING}.dist-info)
endif ()


if (GUILE_FOUND)

  set_source_files_properties (nlopt.i PROPERTIES SWIG_FLAGS "-scmstub")
  set (SWIG_MODULE_nlopt_guile_EXTRA_DEPS nlopt-guile.i generate-cpp)

  # swig_add_module is deprecated
  swig_add_library (nlopt_guile LANGUAGE guile SOURCES nlopt.i)

  target_include_directories (nlopt_guile PRIVATE ${GUILE_INCLUDE_DIRS})
  target_compile_options(nlopt_guile PRIVATE ${GUILE_CFLAGS})

  swig_link_libraries (nlopt_guile ${nlopt_lib})
  target_link_libraries (nlopt_guile ${GUILE_LIBRARIES})
  target_link_options(nlopt_guile PRIVATE ${GUILE_LDFLAGS})

  file (RELATIVE_PATH _REL_GUILE_SITE_PATH ${GUILE_ROOT_DIR} ${GUILE_SITE_DIR})
  set (GUILE_SITE_PATH ${_REL_GUILE_SITE_PATH})
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm DESTINATION ${GUILE_SITE_PATH})

  file (RELATIVE_PATH _REL_GUILE_EXTENSION_PATH ${GUILE_ROOT_DIR} ${GUILE_EXTENSION_DIR})
  set (GUILE_EXTENSION_PATH ${_REL_GUILE_EXTENSION_PATH})
  install (TARGETS nlopt_guile LIBRARY DESTINATION ${GUILE_EXTENSION_PATH})
endif ()
