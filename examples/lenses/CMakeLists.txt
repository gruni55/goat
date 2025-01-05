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
add_executable (lenses "lenses.cpp")

find_library (NAMES goat_maths PATHS ../../lib)
find_library (NAMES goat_raytracing PATHS ../../lib)
find_package (OpenMP)
if (OpenMP_FOUND) 
	   message (STATUS "OpenMP support works")
	   target_compile_options (boxtest PUBLIC ${OpenMP_CXX_FLAGS})
	   target_include_directories (lenses PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	   target_link_libraries (lenses OpenMP::OpenMP_CXX goat_raytracing_mp goat_maths)	
else()
	   message (STATUS "OpenMP support skipped in boxtest")
	   target_include_directories (lenses PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	   target_link_libraries (lenses goat_raytracing goat_maths)	
endif()
install(TARGETS lenses
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    ) 

message (STATUS "LIB Library: ${LIBDIR}")
#add_executable (testbox "testbox.cpp")
#target_include_directories (testbox PUBLIC ${CMAKE_CURRENT_SOURCE} #../../maths ../../raytracing)
#target_link_libraries (testbox goat_raytracing goat_maths)

