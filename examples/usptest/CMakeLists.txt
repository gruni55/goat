# build examples for optical tweezers
cmake_minimum_required (VERSION 2.8.11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type not specified: Use Release by default")
  set(LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Release) 
else() 
  message(STATUS "Build type is set to Debug")
set(LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Debug)
endif(NOT CMAKE_BUILD_TYPE)


message (STATUS "LIB Library: ${LIBDIR}")

find_library (NAMES goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)

add_executable (usptest usptest.cpp )
target_include_directories (usptest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (usptest goat_raytracing goat_maths)

install(TARGETS usptest
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
add_executable (pnexttest pnexttest.cpp "pnexttest.cpp")
target_include_directories (pnexttest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (pnexttest goat_raytracing goat_maths)

install(TARGETS pnexttest
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
