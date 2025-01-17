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
    add_compile_options(-std=c++17)
   endif()

# Set a default build type if none was specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
# Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

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
add_subdirectory (test)

