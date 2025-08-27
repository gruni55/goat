# build examples for optical tweezers

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE "Release")
#  message(STATUS "Build type not specified: Use Release by default")
#  set(LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Release) 
#else() 
#  message(STATUS "Build type is set to Debug")
#set(LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Debug)
#endif(NOT CMAKE_BUILD_TYPE)


message (STATUS "LIB Library: ${LIBDIR}")
#counterpropagating trap
add_executable (cp_trap counterpropagating_trap.cpp)
# target_link_directories(cp_trap PUBLIC ${LIBDIR})
find_library (NAMES goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)

target_include_directories (cp_trap PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (cp_trap goat_raytracing goat_maths)

#one particle in a optical trap
add_executable (otexample example_one_sphere.cpp)
#target_link_directories(otexample PUBLIC ${LIBDIR})
target_include_directories (otexample PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (otexample goat_raytracing goat_maths)

goat_install (TARGETS otexample)
