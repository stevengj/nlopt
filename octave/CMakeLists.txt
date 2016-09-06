
#==============================================================================
# MATLAB/Octave common source files
#==============================================================================

set (NLOPT_SOURCES_M
  NLOPT_GN_DIRECT.m NLOPT_GN_DIRECT_L.m NLOPT_GN_DIRECT_L_RAND.m
  NLOPT_GN_DIRECT_NOSCAL.m NLOPT_GN_DIRECT_L_NOSCAL.m
  NLOPT_GN_DIRECT_L_RAND_NOSCAL.m NLOPT_GN_ORIG_DIRECT.m
  NLOPT_GN_ORIG_DIRECT_L.m NLOPT_GD_STOGO.m NLOPT_GD_STOGO_RAND.m
  NLOPT_LD_LBFGS_NOCEDAL.m NLOPT_LD_LBFGS.m NLOPT_LN_PRAXIS.m
  NLOPT_LD_VAR1.m NLOPT_LD_VAR2.m NLOPT_LD_TNEWTON.m
  NLOPT_LD_TNEWTON_RESTART.m NLOPT_LD_TNEWTON_PRECOND.m
  NLOPT_LD_TNEWTON_PRECOND_RESTART.m NLOPT_GN_CRS2_LM.m NLOPT_GN_MLSL.m
  NLOPT_GD_MLSL.m NLOPT_GN_MLSL_LDS.m NLOPT_GD_MLSL_LDS.m NLOPT_LD_MMA.m
  NLOPT_LN_COBYLA.m NLOPT_LN_NEWUOA.m NLOPT_LN_NEWUOA_BOUND.m
  NLOPT_LN_NELDERMEAD.m NLOPT_LN_SBPLX.m NLOPT_LN_AUGLAG.m NLOPT_LD_AUGLAG.m
  NLOPT_LN_AUGLAG_EQ.m NLOPT_LD_AUGLAG_EQ.m NLOPT_LN_BOBYQA.m
  NLOPT_GN_ISRES.m NLOPT_AUGLAG.m NLOPT_AUGLAG_EQ.m NLOPT_G_MLSL.m
  NLOPT_G_MLSL_LDS.m NLOPT_LD_SLSQP.m NLOPT_LD_CCSAQ.m NLOPT_GN_ESCH.m
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
    LINK_LIBRARIES ${nlopt_lib})
  set_target_properties (nlopt_octave PROPERTIES OUTPUT_NAME nlopt_optimize)

  install (
    TARGETS nlopt_octave
    DESTINATION "${RELATIVE_OCTAVE_OCT_SITE_DIR}"
    COMPONENT "Bindings")
  install (
    FILES ${NLOPT_SOURCES_M}
    DESTINATION "${RELATIVE_OCTAVE_M_SITE_DIR}"
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
    SRC nlopt_optimize-mex.c
    OUTPUT_NAME nlopt_optimize
    LINK_TO ${nlopt_lib})
  set_target_properties (nlopt_matlab PROPERTIES NO_SONAME ON)
endif ()
