cmake_minimum_required (VERSION 2.6)

# define project
project (nlopt_example_cxx CXX)

# find NLopt
# You may need to set NLOPT_DIR variable to the absolute path of the
# directory containing NLoptConfig.cmake configuration file,
# or set CMAKE_PREFIX_PATH to the installation prefix
find_package (NLopt REQUIRED NO_MODULE)
if (NLopt_FOUND)
  message (STATUS "NLopt: ${NLopt_VERSION}")
endif ()

# add NLopt headers location to include path
add_definitions (${NLOPT_DEFINITIONS})
include_directories (${NLOPT_INCLUDE_DIRS})

# declare C++ example target
add_executable (example_cpp example.cpp)
target_link_libraries (example_cpp ${NLOPT_LIBRARIES})
if (UNIX)
  target_link_libraries (example_cpp m)
endif ()
