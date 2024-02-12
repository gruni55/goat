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
	message (STATUS "OpenMP support works")
	add_executable (boxtest boxtest.cpp)
	target_compile_options (boxtest PUBLIC ${OpenMP_CXX_FLAGS})
        target_compile_definitions (boxtest PUBLIC WITH_OCTREE)	
	target_include_directories (boxtest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	target_link_libraries (boxtest OpenMP::OpenMP_CXX goat_raytracing_mp goat_maths)

        add_executable (usp_ellipsoid usp_ellipsoid.cpp)
	target_compile_options (usp_ellipsoid PUBLIC ${OpenMP_CXX_FLAGS})
        target_compile_definitions (usp_ellipsoid PUBLIC WITH_OCTREE)	
	target_include_directories (usp_ellipsoid PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	target_link_libraries (usp_ellipsoid OpenMP::OpenMP_CXX goat_raytracing_mp goat_maths)

  	add_executable (usp_compress usp_compress.cpp)
	target_compile_options (usp_compress PUBLIC ${OpenMP_CXX_FLAGS})
        target_compile_definitions (usp_compress PUBLIC WITH_OCTREE)	
	target_include_directories (usp_compress PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	target_link_libraries (usp_compress OpenMP::OpenMP_CXX goat_raytracing_mp goat_maths)

  	add_executable (usp_lstest usp_lstest.cpp)
	target_compile_options (usp_lstest PUBLIC ${OpenMP_CXX_FLAGS})
        target_compile_definitions (usp_lstest PUBLIC WITH_OCTREE)	
	target_include_directories (usp_lstest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	target_link_libraries (usp_lstest OpenMP::OpenMP_CXX goat_raytracing_mp goat_maths)

	add_executable (usp_line usp_line.cpp)
	target_compile_options (usp_line PUBLIC ${OpenMP_CXX_FLAGS})
        target_compile_definitions (usp_line PUBLIC WITH_OCTREE)	
	target_include_directories (usp_line PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
	target_link_libraries (usp_line OpenMP::OpenMP_CXX goat_raytracing_mp goat_maths)

 endif()

add_executable (usp_raytest usp_raytest.cpp)
target_include_directories (usp_raytest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (usp_raytest goat_raytracing goat_maths)




