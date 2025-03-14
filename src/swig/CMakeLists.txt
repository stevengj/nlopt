
if (POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif ()
# clean up old generated source files before running SWIG, useful for Java
set(UseSWIG_MODULE_VERSION 2)
include (UseSWIG)

# allows one set C++ flags for swig wrappers
set (SWIG_COMPILE_FLAGS "" CACHE STRING "C++ compiler flags used for wrapper code")

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

if (Python_NumPy_FOUND AND NLOPT_CXX)

  set (SWIG_MODULE_nlopt_python_EXTRA_DEPS nlopt-python.i numpy.i generate-cpp)

  swig_add_library (nlopt_python LANGUAGE python SOURCES nlopt.i
                    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/python
                    OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR})

  target_link_libraries (nlopt_python ${nlopt_lib})
  target_link_libraries (nlopt_python Python::NumPy)
  target_link_libraries (nlopt_python Python::Module)

  set_target_properties (nlopt_python PROPERTIES OUTPUT_NAME nlopt)
  set_target_properties (nlopt_python PROPERTIES COMPILE_FLAGS "${SWIG_COMPILE_FLAGS}")

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/python/nlopt.py DESTINATION ${INSTALL_PYTHON_DIR})
  install (TARGETS nlopt_python DESTINATION ${INSTALL_PYTHON_DIR})

  configure_file (METADATA.in METADATA @ONLY)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/METADATA
           DESTINATION ${INSTALL_PYTHON_DIR}/nlopt-${NLOPT_VERSION_STRING}.dist-info)
endif ()


if (GUILE_FOUND AND NLOPT_CXX)

  set (SWIG_MODULE_nlopt_guile_EXTRA_DEPS nlopt-guile.i generate-cpp)
  set (CMAKE_SWIG_FLAGS -scmstub)

  swig_add_library (nlopt_guile LANGUAGE guile SOURCES nlopt.i
                    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/guile
                    OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR})

  set (CMAKE_SWIG_FLAGS)

  target_include_directories (nlopt_guile PRIVATE ${GUILE_INCLUDE_DIRS})

  target_link_libraries (nlopt_guile ${nlopt_lib})
  target_link_libraries (nlopt_guile ${GUILE_LIBRARIES})
  set_target_properties (nlopt_guile PROPERTIES COMPILE_FLAGS "${SWIG_COMPILE_FLAGS}")

  file (RELATIVE_PATH _REL_GUILE_SITE_PATH ${GUILE_ROOT_DIR} ${GUILE_SITE_DIR})
  set (GUILE_SITE_PATH ${_REL_GUILE_SITE_PATH})
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/guile/nlopt.scm DESTINATION ${GUILE_SITE_PATH})

  file (RELATIVE_PATH _REL_GUILE_EXTENSION_PATH ${GUILE_ROOT_DIR} ${GUILE_EXTENSION_DIR})
  set (GUILE_EXTENSION_PATH ${_REL_GUILE_EXTENSION_PATH})
  install (TARGETS nlopt_guile LIBRARY DESTINATION ${GUILE_EXTENSION_PATH})
endif ()


if (JNI_FOUND AND Java_FOUND AND SWIG_FOUND AND NLOPT_CXX)

  include (UseJava)

  set (SWIG_MODULE_nlopt_java_EXTRA_DEPS nlopt-java.i generate-cpp)
  set (CMAKE_SWIG_FLAGS -package nlopt)

  # OUTPUT_DIR is ${CMAKE_CURRENT_BINARY_DIR}/java/ + the -package above (with
  # any '.' replaced by '/'). It must also match the GLOB in glob_java.cmake.
  swig_add_library (nlopt_java LANGUAGE java SOURCES nlopt.i
                               OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/java/nlopt
                               OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR})

  set (CMAKE_SWIG_FLAGS)

  target_link_libraries (nlopt_java ${nlopt_lib})
  target_link_libraries (nlopt_java JNI::JNI)

  set_target_properties (nlopt_java PROPERTIES OUTPUT_NAME nloptjni)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set_target_properties (nlopt_java PROPERTIES
                           COMPILE_OPTIONS "-fno-strict-aliasing")
  endif ()

  # MODULE libraries are always installed to the LIBRARY DESTINATION.
  if (CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
    install (TARGETS nlopt_java LIBRARY DESTINATION ${NLOPT_INSTALL_BINDIR})
  else()
    install (TARGETS nlopt_java LIBRARY DESTINATION ${NLOPT_INSTALL_LIBDIR})
  endif()

  # unfortunately, SWIG will not tell us which .java files it generated, so we
  # have to find out ourselves - this is the only portable way to do so
  # (The nlopt*.i dependencies are there to force updating the list of sources
  # on any changes to the SWIG interface code, they are not direct inputs.)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/java_sources.txt
                     COMMAND ${CMAKE_COMMAND}
                             -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
                             -P ${CMAKE_CURRENT_SOURCE_DIR}/glob_java.cmake
                     DEPENDS "${swig_generated_file_fullname}"
                             nlopt.i nlopt-exceptions.i nlopt-java.i
                             nlopt_java_swig_compilation glob_java.cmake)

  add_jar (nlopt_jar SOURCES @${CMAKE_CURRENT_BINARY_DIR}/java_sources.txt
                     OUTPUT_NAME nlopt)
  install_jar (nlopt_jar ${CMAKE_INSTALL_DATADIR}/java)
endif ()
