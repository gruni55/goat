cmake_minimum_required (VERSION 3.2)

project (goat)
# include (CTest)
if(CMAKE_CXX_COMPILER_LOADED)
     message(STATUS "Compiler path: ${CMAKE_CXX_COMPILER}")
     message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
     message(STATUS "Compiler version:    
             ${CMAKE_CXX_COMPILER_VERSION}")
     message(STATUS "Compiler is part of GCC: 
             ${CMAKE_COMPILER_IS_GNUCXX}")
  #   add_compile_options(-std=c++11)
   endif()


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message (STATUS "Library path is set to: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

add_subdirectory (maths) 
add_subdirectory (raytracing)
add_subdirectory (tinyxml2)
add_subdirectory (examples/ot)
add_subdirectory (examples/paths)
add_subdirectory (examples/layers)
add_subdirectory (examples/axicon)
add_subdirectory (examples/mctest)
add_subdirectory (examples/usptest)
add_subdirectory (examples/sphere)
# add_subdirectory (examples/math)
add_subdirectory (examplesXML)
add_subdirectory (examples/lenses)
add_subdirectory (XML)

