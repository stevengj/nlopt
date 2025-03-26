add_custom_target (tests)

macro(NLOPT_add_cpp_test test_name)
  if (NLOPT_CXX)
    add_executable (${test_name} ${test_name}.cxx)
    target_link_libraries (${test_name} ${nlopt_lib})
    add_dependencies (tests ${test_name})
    target_include_directories (${test_name} PRIVATE ${NLOPT_PRIVATE_INCLUDE_DIRS})
    foreach(arg IN ITEMS ${ARGN})
      add_test (NAME check_${test_name}_${arg} COMMAND ${test_name} ${arg})
      if (CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
        set_tests_properties (check_${test_name}_${arg}
          PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}\\${CMAKE_BUILD_TYPE};$ENV{PATH}")  # to load dll
      endif ()
    endforeach()
  endif ()
endmacro()

NLOPT_add_cpp_test(t_tutorial 24 25 31 40)
NLOPT_add_cpp_test(cpp_functor 0)
NLOPT_add_cpp_test(t_fbound 0)
NLOPT_add_cpp_test(t_except 1 0)
NLOPT_add_cpp_test(t_bounded 0 1 2 3 4 5 6 7 8 19 35 42 43)

# have to add timer.c and mt19937ar.c as symbols are declared extern
set (testopt_sources testfuncs.c testfuncs.h testopt.c ${PROJECT_SOURCE_DIR}/src/util/timer.c ${PROJECT_SOURCE_DIR}/src/util/mt19937ar.c)
if (NOT HAVE_GETOPT OR NOT HAVE_GETOPT_H)
  list (APPEND testopt_sources ${PROJECT_SOURCE_DIR}/src/util/nlopt-getopt.c)
endif ()
add_executable (testopt ${testopt_sources})
target_link_libraries (testopt ${nlopt_lib})
target_include_directories (testopt PRIVATE ${NLOPT_PRIVATE_INCLUDE_DIRS})
add_dependencies (tests testopt)

if (NLOPT_CXX)
  set_target_properties(testopt PROPERTIES LINKER_LANGUAGE CXX)
endif ()

foreach (algo_index RANGE 28) # 42
  foreach (obj_index RANGE 1) # 21
    set (enable_ TRUE)
    # cxx ags/stogo
    if (NOT NLOPT_CXX)
      if (algo_index STREQUAL 8 OR algo_index STREQUAL 9 OR algo_index STREQUAL 42)
        set (enable_ FALSE)
      endif ()
    endif ()
    # nocedal
    if (algo_index STREQUAL 10)
     set (enable_ FALSE)
    endif ()
    if (enable_)
      add_test (NAME testopt_algo${algo_index}_obj${obj_index} COMMAND testopt -r 0 -a ${algo_index} -o ${obj_index})
      if (CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
        set_tests_properties (testopt_algo${algo_index}_obj${obj_index} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}\\${CMAKE_BUILD_TYPE};$ENV{PATH}")  # to load dll
      endif ()
      # Check if LUKSAN targets are available.
      if (NOT NLOPT_LUKSAN)
        set (list_of_algorithms_requiring_luksan 10 11 12 13 14 15 16 17 18)
        if (algo_index IN_LIST list_of_algorithms_requiring_luksan)
          set_tests_properties (testopt_algo${algo_index}_obj${obj_index} PROPERTIES DISABLED TRUE)
        endif()
      endif()
    endif ()
  endforeach ()
endforeach ()

if (TARGET nlopt_python)
  set (PYINSTALLCHECK_ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/src/swig"
                                  "PYTHONPATH=$<SHELL_PATH:${PROJECT_BINARY_DIR}/src/swig/python;${PROJECT_BINARY_DIR}/src/swig>"
      )

  foreach (algo_index 24 25 31 40)
    add_test (NAME test_python${algo_index} COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/t_python.py ${algo_index})
    set_tests_properties (test_python${algo_index} PROPERTIES ENVIRONMENT "${PYINSTALLCHECK_ENVIRONMENT}")
  endforeach()

  add_test (NAME test_memoize COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/t_memoize.py ${algo_index})
  set_tests_properties (test_memoize PROPERTIES ENVIRONMENT "${PYINSTALLCHECK_ENVIRONMENT}")
endif ()

# The test uses lambdas and the :: operator, which were introduced in Java 1.8.
# The binding itself should also compile with more ancient Java versions that
# have already reached their end of life, but it is not worth uglifying the test
# code for them, because people will then cargo-cult the legacy boilerplate.
if (TARGET nlopt_java AND NOT Java_VERSION VERSION_LESS 1.8)
  include (UseJava)
  add_jar (t_java SOURCES t_java.java INCLUDE_JARS nlopt_jar ENTRY_POINT t_java)
  get_property (t_java_jar TARGET t_java PROPERTY JAR_FILE)
  get_property (nlopt_jar_jar TARGET nlopt_jar PROPERTY JAR_FILE)
  set (nlopt_java_dir $<TARGET_FILE_DIR:nlopt_java>)
  foreach (algo_index 24 25 31 40)
    add_test (NAME test_java${algo_index}
              COMMAND ${Java_JAVA_EXECUTABLE}
                      -cp "$<SHELL_PATH:${t_java_jar};${nlopt_jar_jar}>"
                      -Djava.library.path=${nlopt_java_dir} t_java
                      ${algo_index})
  endforeach()
endif ()

if (TARGET nlopt_optimize)
  add_test (NAME test_octave COMMAND ${OCTAVE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/t_octave.m ${PROJECT_SOURCE_DIR}/src/octave ${PROJECT_BINARY_DIR}/src/octave)
endif ()

if (TARGET nlopt_optimize-mex)
  add_test (NAME test_matlab COMMAND ${Matlab_MAIN_PROGRAM} -nodesktop -nosplash  -r "addpath('${PROJECT_SOURCE_DIR}/src/octave'); addpath('${PROJECT_BINARY_DIR}/src/octave'); try; run('${CMAKE_CURRENT_SOURCE_DIR}/t_matlab.m'); catch; exit(1); end; quit")
endif ()

if (TARGET nlopt_guile)
  set (GUILECHECK_ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/src/swig"
                              "GUILE_LOAD_PATH=$<SHELL_PATH:${PROJECT_BINARY_DIR}/src/swig/guile;${PROJECT_BINARY_DIR}/src/swig>"
                              "GUILE_AUTO_COMPILE=0")
  add_test (NAME test_guile
            COMMAND ${GUILE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/t_guile.scm
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/src/swig)
  set_tests_properties (test_guile PROPERTIES ENVIRONMENT "${GUILECHECK_ENVIRONMENT}")
endif ()

if (NLOPT_FORTRAN)
  add_executable (t_fortran t_fortran.f90)
  target_link_libraries (t_fortran ${nlopt_lib})
  target_include_directories (t_fortran PRIVATE ${NLOPT_PRIVATE_INCLUDE_DIRS})
  add_test (NAME test_fortran COMMAND t_fortran)
endif ()
