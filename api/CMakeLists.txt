
#==============================================================================
# Fortran interface
#==============================================================================

# generate nlopt.f from nlopt.h enums
if (BUILD_F77 AND "${CMAKE_CURRENT_SOURCE_DIR}/nlopt.h" IS_NEWER_THAN
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt.f")
  file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/nlopt.f" "")
  file (STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt.h" NLOPT_H_LINES
    REGEX "    NLOPT_[A-Z0-9_]+")
  set (i 0)
  foreach (NLOPT_H_LINE ${NLOPT_H_LINES})
    if (NOT NLOPT_H_LINE MATCHES "NLOPT_NUM_")
      string (REGEX REPLACE ".*NLOPT_([A-Z0-9_]+).*" "\\1"
        ENUM_STRING ${NLOPT_H_LINE})
      string (REGEX REPLACE ".*NLOPT_[A-Z0-9_]+ = (-?[0-9]+).*" "\\1"
        ENUM_VAL ${NLOPT_H_LINE})
      if (ENUM_VAL MATCHES "^-?[0-9]+$")
        set (i ${ENUM_VAL})
      endif ()
      set (ENUM_LINE "      integer NLOPT_${ENUM_STRING}\n      parameter (NLOPT_${ENUM_STRING}=${i})\n")
      file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/nlopt.f" "${ENUM_LINE}")

      # https://public.kitware.com/Bug/print_bug_page.php?bug_id=8996
      if (i MATCHES "^-")
        math (EXPR i "1 ${i}")
      else ()
        math (EXPR i "${i} + 1")
      endif ()
    endif ()
  endforeach ()

  # cleanup
  unset (ENUM_LINE)
  unset (ENUM_VAL)
  unset (ENUM_STRING)
  unset (i)
  unset (NLOPT_H_LINE)
  unset (NLOPT_H_LINES)

  # mark file as generated
  set_source_files_properties ("${CMAKE_CURRENT_BINARY_DIR}/nlopt.f"
    PROPERTIES GENERATED TRUE)
endif ()

#==============================================================================
# C++ interface
#==============================================================================

# generate nlopt.hpp from nlopt-in.hpp
if (BUILD_CPP AND "${CMAKE_CURRENT_SOURCE_DIR}/nlopt-in.hpp" IS_NEWER_THAN
    "${CMAKE_CURRENT_BINARY_DIR}/nlopt.hpp")
  set (OUTPUT_HPP "${CMAKE_CURRENT_BINARY_DIR}/nlopt.hpp")
  file (WRITE "${OUTPUT_HPP}" "")
  file (STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt-in.hpp" NLOPT_HPP_LINES)
  foreach (NLOPT_HPP_LINE ${NLOPT_HPP_LINES})
    list (LENGTH NLOPT_HPP_LINE line_len)
    # handling trailing backlashes in "file (STRINGS" is a little tricky
    if (line_len VERSION_LESS 8)
      file (APPEND "${OUTPUT_HPP}" "${NLOPT_HPP_LINE}\n")
    else ()
      set (prev_inst FALSE)
      foreach (NLOPT_HPP_SUBLINE ${NLOPT_HPP_LINE})
        # test is we need to add the eaten semicolon
        if (NLOPT_HPP_SUBLINE MATCHES "\\)$" OR NLOPT_HPP_SUBLINE MATCHES "return")
          set (new_inst TRUE)
        else ()
          set (new_inst FALSE)
        endif ()
        if (NOT prev_inst)
          file (APPEND "${OUTPUT_HPP}" "${NLOPT_HPP_SUBLINE}")
          if (new_inst)
            file (APPEND "${OUTPUT_HPP}" ";")
          endif ()
          list (FIND NLOPT_HPP_LINE "${NLOPT_HPP_SUBLINE}" index)
          math (EXPR index "${index} + 1")
          list (LENGTH NLOPT_HPP_LINE total)
          if (NOT index STREQUAL total)
            file (APPEND "${OUTPUT_HPP}" " \\")
          endif ()
          file (APPEND "${OUTPUT_HPP}" "\n")
        endif ()
        set (prev_inst ${new_inst})
      endforeach ()
    endif ()
    if (NLOPT_HPP_LINE MATCHES "GEN_ENUMS_HERE")
      file (APPEND "${OUTPUT_HPP}" "  enum algorithm {\n")
      file (STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt.h" NLOPT_H_LINES
        REGEX "    NLOPT_[A-Z0-9_]+")
      foreach (NLOPT_H_LINE ${NLOPT_H_LINES})
        string (REGEX REPLACE "NLOPT_" "" ENUM_LINE ${NLOPT_H_LINE})
        file (APPEND "${OUTPUT_HPP}" "${ENUM_LINE}\n")
        if (NLOPT_H_LINE MATCHES "NLOPT_NUM_ALGORITHMS")
          file (APPEND "${OUTPUT_HPP}" "  };\n  enum result {\n")
        elseif (NLOPT_H_LINE MATCHES "NLOPT_MAXTIME_REACHED")
          file (APPEND "${OUTPUT_HPP}" "  };\n")
        endif ()
      endforeach ()
    endif ()
  endforeach ()

  # cleanup
  unset (ENUM_LINE)
  unset (NLOPT_H_LINE)
  unset (NLOPT_H_LINES)
  unset (NLOPT_HPP_LINE)
  unset (NLOPT_HPP_LINES)
  unset (NLOPT_HPP_SUBLINE)
  unset (total)
  unset (index)
  unset (line_len)
  unset (new_inst)
  unset (prev_inst)
  unset (OUTPUT_HPP)

  # mark file as generated
  set_source_files_properties ("${CMAKE_CURRENT_BINARY_DIR}/nlopt.hpp"
    PROPERTIES GENERATED TRUE)
endif ()
