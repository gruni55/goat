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
#counterpropagating trap
add_executable (layers layers.cpp)
add_executable (testlayer testlayer.cpp)
# target_link_directories(cp_trap PUBLIC ${LIBDIR})
find_library (NAMES goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)

target_include_directories (layers PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (layers goat_raytracing goat_maths)



target_include_directories (testlayer PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (testlayer goat_raytracing goat_maths)


install(TARGETS layers testlayer
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
