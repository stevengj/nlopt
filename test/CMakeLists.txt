
add_executable ( testopt EXCLUDE_FROM_ALL testfuncs.c testfuncs.h testopt.cpp )
target_link_libraries ( testopt nlopt )
