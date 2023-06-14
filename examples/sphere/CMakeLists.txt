# build examples for sphere calculation

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
add_executable (sphere "sphere.cpp")
if (WITHOUT_OPENMP) 
	message (STATUS "OpenMP support skipped due to parameter \"WITHOUT_OPENMP\"")
        target_link_libraries (sphere goat_raytracing goat_maths)
else()
	find_package(OpenMP)
	if (OPENMP_FOUND)
 		message (STATUS "OpenMP found")
	        target_compile_options (goat_raytracing PUBLIC ${OpenMP_CXX_FLAGS})
	        target_link_libraries (sphere goat_raytracing goat_maths OpenMP::OpenMP_CXX)
	else()
		message (WARNING "OpenMP not found - Calculation without parallelization")
	        target_link_libraries (sphere goat_raytracing goat_maths)
	endif()
endif()

find_library (NAMES goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)

target_include_directories (sphere PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)

