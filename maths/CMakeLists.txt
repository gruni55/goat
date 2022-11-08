#CMAKE configuration file for vector math library

set (MATHS_HDRS fresnel.h   goodies.h   matrix.h    vector.h)
set (MATHS_SRCS fresnel.cpp   goodies.cpp   matrix.cpp    vector.cpp)
add_library (goat_maths ${MATHS_SRCS} ${MATHS_HDRS})
target_include_directories (goat_maths PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
