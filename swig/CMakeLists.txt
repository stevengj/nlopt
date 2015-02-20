option (BUILD_PYTHON "build python bindings" ON)

find_package (PythonInterp)
find_package (PythonLibs)
find_package (SWIG)

if (BUILD_PYTHON AND PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND AND SWIG_FOUND AND BUILD_SHARED_LIBS)
  include (UseSWIG)
  include_directories (${CMAKE_SOURCE_DIR}/api)
  include_directories (${PYTHON_INCLUDE_DIRS})

  set_source_files_properties (nlopt.i PROPERTIES CPLUSPLUS ON)
  set (SWIG_MODULE_nlopt_EXTRA_DEPS nlopt-python.i numpy.i ${CMAKE_SOURCE_DIR}/api/nlopt.hpp)

  swig_add_module (nlopt python nlopt.i)
  swig_link_libraries (nlopt nlopt)

  execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
                    OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

  get_filename_component (_ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE)
  file (RELATIVE_PATH _REL_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH})

  set (PYTHON_MODULE_PATH ${_REL_PYTHON_MODULE_PATH})

  install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/nlopt.py DESTINATION ${PYTHON_MODULE_PATH} )
  install ( TARGETS ${SWIG_MODULE_nlopt_REAL_NAME} LIBRARY DESTINATION ${PYTHON_MODULE_PATH} )
endif ()
