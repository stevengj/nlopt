set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    
string(ASCII 27 Esc)
set(Green    "${Esc}[92m"     CACHE STRING "green color"    FORCE)
set(Blue     "${Esc}[1;94m"   CACHE STRING "Blue color"     FORCE)
set(BoldRed  "${Esc}[1;5;91m" CACHE STRING "bold-red color" FORCE)
set(NoColor  "${Esc}[0m"      CACHE STRING "reset color"    FORCE)

if(NOT WIN32)

    execute_process(COMMAND   ${CMAKE_C_COMPILER}  -dumpversion   OUTPUT_VARIABLE   GCC_VERSION)
    if(GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
        set(STDLIBC  "-std=c++1y")
    elseif(GCC_VERSION VERSION_EQUAL 4.8)
        set(STDLIBC  "-std=c++11")
    else()
        set(STDLIBC  " ")
    endif()
    set(STDLIBCXX_FLAG  ${STDLIBC}  CACHE STRING "stdlibc++ flag")
    message(STATUS "${Green}GCC_VERSION= ${GCC_VERSION} ")
    message(STATUS " -std= flag: ${STDLIBCXX_FLAG}  ${NoColor} " )

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   ${STDLIBCXX_FLAG} " )

    execute_process( COMMAND  mkdir -p ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

    enable_testing()
    
    if(BUILD_TESTS)
        message(STATUS "${Blue} Testing enabled: run 'ctest -VV -N' to see list of tests available and their arguments. ${NoColor}")
    endif()
   
endif()
