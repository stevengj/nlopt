add_custom_target (tests)

add_executable (testopt EXCLUDE_FROM_ALL testfuncs.c testfuncs.h testopt.cpp)
target_link_libraries (testopt nlopt)
add_dependencies (tests testopt)

foreach (algo_index RANGE 29)# 42
  foreach (obj_index RANGE 1)# 21
    set (enable_ TRUE)
    # cxx stogo
    if (NOT WITH_CXX)
      if (algo_index STREQUAL 8 OR algo_index STREQUAL 9)
        set (enable_ FALSE)
      endif ()
    endif ()
    # L-BFGS
    if (algo_index STREQUAL 10)
      set (enable_ FALSE)
    endif ()
    if (enable_)
      add_test (NAME testopt_algo${algo_index}_obj${obj_index} COMMAND testopt -a ${algo_index} -o ${obj_index})
    endif ()
  endforeach ()
endforeach ()

if (PYTHONINTERP_FOUND)
  set (PYINSTALLCHECK_ENVIRONMENT "LD_LIBRARY_PATH=${INSTALL_DESTDIR}${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}"
                                  "PYTHONPATH=${INSTALL_DESTDIR}${CMAKE_INSTALL_PREFIX}/${PYTHON_MODULE_PATH}${PATH_SEP}$ENV{PYTHONPATH}"
    )
  add_test (NAME test_std_python COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_std.py)
  set_tests_properties (test_std_python PROPERTIES ENVIRONMENT "${PYINSTALLCHECK_ENVIRONMENT}")
endif ()
