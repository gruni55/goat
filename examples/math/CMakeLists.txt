# build examples for mathematic library

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

add_executable (testmath testmath.cpp)
target_include_directories (testmath PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths)
target_link_libraries (testmath goat_maths)

install(TARGETS testmath
                RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/${CMAKE_BUILD_TYPE}
                LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                )

