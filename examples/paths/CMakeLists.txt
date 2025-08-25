# here all entries for the test program raytest

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE "Release")
#  message(STATUS "Build type not specified: Use Release by default")
#else() 
#  message(STATUS "Build type is set to Debug")
#endif(NOT CMAKE_BUILD_TYPE)

add_executable (raytest raytest.cpp)
target_compile_definitions (raytest PUBLIC CALC_PATH)
target_include_directories (raytest PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (raytest goat_raytracing goat_maths)

add_executable (stlexample stlexample.cpp)
target_compile_definitions (stlexample PUBLIC CALC_PATH)
target_include_directories (stlexample PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (stlexample goat_raytracing goat_maths)

add_executable (conepath conepath.cpp)
target_compile_definitions (conepath PUBLIC CALC_PATH)
target_include_directories (conepath PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (conepath goat_raytracing goat_maths)

add_executable (mirror mirror.cpp)
target_compile_definitions (mirror PUBLIC CALC_PATH)
target_include_directories (mirror PUBLIC ${CMAKE_CURRENT_SOURCE} ../../maths ../../raytracing)
target_link_libraries (mirror goat_raytracing goat_maths)

install(TARGETS mirror conepath stlexample raytest 
                RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/$<CONFIG>
                LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                )

