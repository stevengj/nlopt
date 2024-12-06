# This file(GLOB ...) must run at build (make) time, after the SWIG run. So it
# cannot be invoked directly from CMakeLists.txt, but must be invoked through
# cmake -P at the correct spot of the build, using add_custom_command.
file(GLOB JAVA_SOURCES ${BINARY_DIR}/java/nlopt/*.java)
list(JOIN JAVA_SOURCES "\n" JAVA_SOURCES_LINES)
file(WRITE ${BINARY_DIR}/java_sources.txt ${JAVA_SOURCES_LINES})

# SWIG hardcodes non-vararg initial elements for std::vector wrappers,
# probably to support Java versions older than 1.5. We do not really care
# about supporting a Java that old, so fix the generated code.
# See: https://github.com/swig/swig/issues/3085
file(READ ${BINARY_DIR}/java/nlopt/DoubleVector.java FILE_CONTENTS)
string(REPLACE "double[] initialElements" "double... initialElements" FILE_CONTENTS "${FILE_CONTENTS}")
file(WRITE ${BINARY_DIR}/java/nlopt/DoubleVector.java "${FILE_CONTENTS}")
