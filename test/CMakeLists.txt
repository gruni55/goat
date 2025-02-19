# build examples for axicon calculation

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


find_library (goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)


# settings for fieldtest
add_executable (fieldtest "fieldtest.cpp")
target_include_directories (fieldtest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)

find_package(OpenMP)
if ((NOT OPENMP_FOUND) OR (WITHOUT_OPENMP))
   target_link_libraries (fieldtest goat_raytracing goat_maths)	

else()
   target_link_libraries (fieldtest goat_raytracing_mp goat_maths OpenMP::OpenMP_CXX)	
endif()

# settings for polarisationtest
add_executable(polarisationtest polarisationtest.cpp)
target_include_directories (polarisationtest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (polarisationtest goat_raytracing goat_maths)	



install(TARGETS fieldtest polarisationtest
                RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/${CMAKE_BUILD_TYPE}
                LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                )
