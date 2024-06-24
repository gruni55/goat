# build examples for sphere calculation

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
add_executable (readXML ../tinyxml2/tinyxml2.cpp "readXML.cpp")


find_package (OpenMP)
find_library (NAMES goat_maths PATHS ../lib ${LIBDIR})
find_library (NAMES goat_raytracing PATHS ../lib ${LIBDIR})
find_library (NAMES goat_xml PATHS ../lib ${LIBDIR})


target_compile_options (readXML PUBLIC ${OpenMP_CXX_FLAGS})
target_compile_definitions (readXML PUBLIC WITH_OCTREE WITH_SUPERGITTER)
target_include_directories (readXML PUBLIC ${CMAKE_CURRENT_SOURCE} ../maths ../raytracing ../tinyxml2 ../XML)
target_link_libraries (readXML OpenMP::OpenMP_CXX goat_xml goat_raytracing_mp goat_maths)
install(TARGETS readXML
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
                )

