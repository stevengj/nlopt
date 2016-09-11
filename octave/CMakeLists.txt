
#==============================================================================
# MATLAB/Octave common source files
#==============================================================================

set (NLOPT_SOURCES_M
  nlopt_algorithm.m
  nlopt_minimize.m nlopt_minimize_constrained.m nlopt_optimize.m)

#==============================================================================
# Octave OCT extension
#==============================================================================

# generate nlopt_optimize_usage.h (docstring) from nlopt_optimize.m
if (BUILD_OCTAVE AND "${CMAKE_CURRENT_SOURCE_DIR}/nlopt_optimize.m" IS_NEWER_THAN
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt_optimize_usage.h")
  file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/nlopt_optimize_usage.h"
    "#define NLOPT_OPTIMIZE_USAGE \\\n")
  file (STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt_optimize.m" INPUT_LINES)
  foreach (INPUT_LINE ${INPUT_LINES})
    string (REGEX REPLACE "^% " "" INPUT_LINE ${INPUT_LINE})
    string (REGEX REPLACE "^%" "" INPUT_LINE ${INPUT_LINE})
    foreach (repl_expr "tolerance" "help NLOPT_LN_SBPLX" "population")
      string (REGEX REPLACE "\"${repl_expr}\"" "${repl_expr}"
        INPUT_LINE "${INPUT_LINE}")
    endforeach ()
    file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/nlopt_optimize_usage.h"
      "  \"${INPUT_LINE}\\n\" \\\n")
  endforeach ()
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/nlopt_optimize_usage.h" "  \"\"\n")

  # cleanup
  unset (repl_expr)
  unset (INPUT_LINE)
  unset (INPUT_LINES)

  # set file as generated
  set_source_files_properties (
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt_optimize_usage.h"
    PROPERTIES GENERATED TRUE)
endif ()

if (OCTAVE_FOUND)
  if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories (BEFORE ${OCTAVE_INCLUDE_DIRS})
  endif ()
  include_directories (${CMAKE_CURRENT_BINARY_DIR})

  # Octave lib target (OCT-file)
  octave_add_oct (nlopt_octave
    SOURCES
      nlopt_optimize-oct.cc
      "${CMAKE_CURRENT_BINARY_DIR}/nlopt_optimize_usage.h"
      "${CMAKE_SOURCE_DIR}/api/nlopt.h"
    LINK_LIBRARIES ${NLOPT_LIB})
  set_target_properties (nlopt_octave PROPERTIES OUTPUT_NAME nlopt_optimize)

  if (WITH_FOLDERS)
    set_target_properties (nlopt_octave PROPERTIES FOLDER "Octave")
  endif ()

  install (
    TARGETS nlopt_octave
    DESTINATION "${INSTALL_OCTAVE_OCT_DIR}"
    COMPONENT "Bindings")
  install (
    FILES ${NLOPT_SOURCES_M}
    DESTINATION "${INSTALL_OCTAVE_M_DIR}"
    COMPONENT "Bindings")
endif ()

#==============================================================================
# MATLAB MEX extension
#==============================================================================

if (Matlab_FOUND)
  # required for the matlab_add_mex macro
  cmake_minimum_required (VERSION 3.3)

  # MATLAB lib target (MEX-file)
  matlab_add_mex (NAME nlopt_matlab
    SRC
      nlopt_optimize-mex.c
      "${CMAKE_SOURCE_DIR}/api/nlopt.h"
    OUTPUT_NAME nlopt_optimize
    LINK_TO ${NLOPT_LIB})
  set_target_properties (nlopt_matlab PROPERTIES NO_SONAME ON)

  if (WITH_FOLDERS)
    set_target_properties (nlopt_matlab PROPERTIES FOLDER "MATLAB")
  endif ()

  install (
    TARGETS nlopt_matlab
    DESTINATION "${INSTALL_MATLAB_DIR}"
    COMPONENT "Bindings")
  install (
    FILES ${NLOPT_SOURCES_M}
    DESTINATION "${INSTALL_MATLAB_DIR}"
    COMPONENT "Bindings")
endif ()
