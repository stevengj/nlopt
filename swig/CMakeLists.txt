
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-enum-renames.i)
  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i "// AUTOMATICALLY GENERATED -- DO NOT EDIT\n")
  file (STRINGS ${CMAKE_SOURCE_DIR}/api/nlopt.h NLOPT_H_LINES REGEX "    NLOPT_[A-Z0-9_]+")
  foreach (NLOPT_H_LINE ${NLOPT_H_LINES})
    string (REGEX REPLACE ".*NLOPT_([A-Z0-9_]+).*" "%rename(NLOPT_\\1) nlopt::\\1;\n" ENUM_LINE ${NLOPT_H_LINE})
    file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/nlopt-enum-renames.i "${ENUM_LINE}")
  endforeach ()
endif ()

if (BUILD_OCTAVE)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-python.cpp)
    find_package (SWIG)
  endif ()
endif ()

if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND AND (SWIG_FOUND OR (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-python.cpp)))

  include_directories (${PYTHON_INCLUDE_DIRS})

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-python.cpp)
    set (SWIG_MODULE_nlopt_REAL_NAME _nlopt)
    add_library (${SWIG_MODULE_nlopt_REAL_NAME} MODULE nlopt-python.cpp)
    target_link_libraries(${SWIG_MODULE_nlopt_REAL_NAME} ${nlopt_lib})
    target_link_libraries(${SWIG_MODULE_nlopt_REAL_NAME} ${PYTHON_LIBRARIES})
    set_target_properties (${SWIG_MODULE_nlopt_REAL_NAME} PROPERTIES PREFIX "")
    if (NOT CMAKE_SOURCE_DIR MATCHES CMAKE_BINARY_DIR)
      execute_process (COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/nlopt.py ${CMAKE_CURRENT_BINARY_DIR}/nlopt.py)
    endif ()
  else ()
    include (UseSWIG)
    set_source_files_properties (nlopt.i PROPERTIES CPLUSPLUS ON)
    set (SWIG_MODULE_nlopt_EXTRA_DEPS nlopt-python.i numpy.i)

    swig_add_module (nlopt python nlopt.i)
    swig_link_libraries (nlopt ${nlopt_lib})
    swig_link_libraries (nlopt ${PYTHON_LIBRARIES})
  endif ()

  install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.py DESTINATION ${PYTHON_MODULE_PATH} )
  install ( TARGETS ${SWIG_MODULE_nlopt_REAL_NAME} LIBRARY DESTINATION ${PYTHON_MODULE_PATH} )
endif ()

if (BUILD_GUILE)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-guile.cpp)
    find_package (SWIG)
  endif ()
endif ()

# guile bindings with gcc only ok with swig >= 2.0.10
if (GUILE_FOUND AND ((SWIG_FOUND AND SWIG_VERSION VERSION_GREATER 2.0.9) OR (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-guile.cpp)))
  include_directories (${GUILE_INCLUDE_DIRS})

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-guile.cpp)
    set (guile_cpp_source ${CMAKE_CURRENT_SOURCE_DIR}/nlopt-guile.cpp)
  else ()
    file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/GenericConfigureFile.cmake "configure_file (\${INFILE} \${OUTFILE} @ONLY)")

    set (guile_cpp_source ${CMAKE_CURRENT_BINARY_DIR}/nlopt-guile.cpp)
    add_custom_command (OUTPUT ${guile_cpp_source} nlopt.scm.in
                        COMMAND ${SWIG_EXECUTABLE} -I${CMAKE_BINARY_DIR}/api -outdir ${CMAKE_CURRENT_BINARY_DIR} -c++ -guile -scmstub -o ${guile_cpp_source} ${CMAKE_CURRENT_SOURCE_DIR}/nlopt.i
                        COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm ${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm.in
                        COMMAND ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm -DNLOPT_SUFFIX=${NLOPT_SUFFIX} -P ${CMAKE_CURRENT_BINARY_DIR}/GenericConfigureFile.cmake)
  endif ()

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nlopt.scm.in)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/nlopt.scm.in ${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm @ONLY)
  endif ()

  set (SWIG_MODULE_nlopt_guile_REAL_NAME nlopt_guile)
  add_library (${SWIG_MODULE_nlopt_guile_REAL_NAME} MODULE ${guile_cpp_source})
  target_link_libraries(${SWIG_MODULE_nlopt_guile_REAL_NAME} ${nlopt_lib})
  target_link_libraries(${SWIG_MODULE_nlopt_guile_REAL_NAME} ${GUILE_LIBRARIES})

  set_target_properties (${SWIG_MODULE_nlopt_guile_REAL_NAME} PROPERTIES SOVERSION ${_INFO_SOVERSION})
  set_target_properties (${SWIG_MODULE_nlopt_guile_REAL_NAME} PROPERTIES VERSION ${_INFO_VERSION} )

  file (RELATIVE_PATH _REL_GUILE_SITE_PATH ${GUILE_ROOT_DIR} ${GUILE_SITE_DIR})
  set (GUILE_SITE_PATH ${_REL_GUILE_SITE_PATH})
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.scm DESTINATION ${GUILE_SITE_PATH})

  file (RELATIVE_PATH _REL_GUILE_EXTENSION_PATH ${GUILE_ROOT_DIR} ${GUILE_EXTENSION_DIR})
  set (GUILE_EXTENSION_PATH ${_REL_GUILE_EXTENSION_PATH})
  install (TARGETS ${SWIG_MODULE_nlopt_guile_REAL_NAME} LIBRARY DESTINATION ${GUILE_EXTENSION_PATH})
endif ()
