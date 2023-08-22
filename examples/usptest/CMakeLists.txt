# build examples for optical tweezers

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

add_executable (usptest_spatial usptest_spatial.cpp)
target_include_directories (usptest_spatial PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (usptest_spatial goat_raytracing goat_maths)

install(TARGETS usptest_spatial
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )


add_executable (pnexttest pnexttest.cpp )
target_include_directories (pnexttest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (pnexttest goat_raytracing goat_maths)

install(TARGETS pnexttest
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

find_package (OpenMP)
if (WITHOUT_MP) 
	message (STATUS "OpenMP support skipped in boxtest")
else()
	message (STATUS "OpenMP support for boxtest ==> boxtest will be build")
	add_executable (boxtest boxtest.cpp)
	target_compile_options (boxtest PUBLIC ${OpenMP_CXX_FLAGS})
        target_compile_definitions (boxtest PUBLIC WITH_OCTREE)	
	target_include_directories (boxtest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	target_link_libraries (boxtest goat_raytracing_mp goat_maths)
endif()
install(TARGETS boxtest
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

