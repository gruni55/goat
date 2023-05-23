# build examples for sphere calculation
cmake_minimum_required (VERSION 2.8.11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type not specified: Use Release by default")
  set(LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}) 
else() 
  message(STATUS "Build type is set to Debug")
set(LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif(NOT CMAKE_BUILD_TYPE)


message (STATUS "LIB Library: ${LIBDIR}")
add_executable (readXML "readXML.cpp")

find_library (NAMES goat_maths PATHS ../../lib ${LIBDIR})
find_library (NAMES goat_raytracing PATHS ../../lib ${LIBDIR})
find_library (NAMES goat_xml PATHS ../../lib ${LIBDIR})
find_library (NAMES tinyxml2 PATHS ../../lib  ${LIBDIR})

target_include_directories (readXML PUBLIC ${CMAKE_CURRENT_SOURCE} ../maths ../raytracing ../tinyxml2)
target_link_libraries (readXML goat_raytracing goat_maths goat_xml tinyxml2)
