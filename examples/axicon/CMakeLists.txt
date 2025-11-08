# build examples for axicon calculation

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
add_executable (axicon "axicon.cpp")
add_executable (axicon_kirchhoff "axicon_kirchhoff.cpp")
add_executable (axicon2 "axicon2.cpp")
add_executable (axicon_with_box "axicon_with_box.cpp")


find_library (NAMES goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)


target_include_directories (axicon PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_include_directories (axicon_kirchhoff PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_include_directories (axicon2 PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_include_directories (axicon_with_box PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)

find_package(OpenMP)
if ((NOT OPENMP_FOUND) OR (WITHOUT_OPENMP))
   target_link_libraries (axicon goat_raytracing goat_maths)
   target_link_libraries (axicon_kirchhoff goat_raytracing goat_maths)
   target_link_libraries (axicon2 goat_raytracing goat_maths)
   target_link_libraries (axicon_with_box goat_raytracing goat_maths)

else()
   target_link_libraries (axicon goat_raytracing_mp goat_maths OpenMP::OpenMP_CXX)
   target_link_libraries (axicon_kirchhoff goat_raytracing_mp goat_maths OpenMP::OpenMP_CXX)
   target_link_libraries (axicon2 goat_raytracing_mp goat_maths OpenMP::OpenMP_CXX)
   target_link_libraries (axicon_with_box goat_raytracing_mp goat_maths OpenMP::OpenMP_CXX)
endif()

goat_install (TARGETS axicon2 axicon_with_box axicon_kirchhoff)
