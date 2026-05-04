cmake_minimum_required (VERSION 3.5)
project (goat
   VERSION 1.0.0
    LANGUAGES CXX
)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


include(GNUInstallDirs)

# Subdir-Genex für Debug/Release
set(GOAT_CFG_SUBDIR
  "$<IF:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>,Release,Debug>"
)

option(WITH_FFTW "Enable FFTW support" ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")



function(goat_install)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs TARGETS)
  cmake_parse_arguments(GOATINST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT GOATINST_TARGETS)
    message(FATAL_ERROR "goat_install() benötigt TARGETS")
  endif()

  install(TARGETS ${GOATINST_TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${GOAT_CFG_SUBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${GOAT_CFG_SUBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${GOAT_CFG_SUBDIR}
  )
endfunction()


set(GOAT_CFG_SUBDIR
  "$<IF:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>,Release,Debug>"
)

option(WITH_OPENMP      "Enable OpenMP support"         ON)
option(WITH_SUPERGITTER "Enable Supergitter features"   ON)
option(WITH_OCTREE      "Enable Octree features"        ON)

if(CMAKE_CXX_COMPILER_LOADED)
     message(STATUS "Compiler path: ${CMAKE_CXX_COMPILER}")
     message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
     message(STATUS "Compiler version:    
             ${CMAKE_CXX_COMPILER_VERSION}")
     message(STATUS "Compiler is part of GCC: 
             ${CMAKE_COMPILER_IS_GNUCXX}")
  #  add_compile_options(-std=c++17)
   endif()

if(UNIX)
      if (NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
          set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release")
      endif()
endif()
# include (CTest)



#if (NOT CMAKE_CONFIGURATION_TYPES)
#    # Set a default build type if none was specified
#     if (CMAKE_BUILD_TYPE STREQUAL "")
#        set(default_build_type "Release")
#        #if(NOT DEFINED CMAKE_BUILD_TYPE AND NOT DEFINED CMAKE_CONFIGURATION_TYPES)
#          message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
#          set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
#              STRING "Choose the type of build." )
#        # Set the possible values of build type for cmake-gui
#          set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
#            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
#        #endif()
#    endif()
#endif()

if(WITH_FFTW)
    include(cmake/FindFFTW.cmake)
endif()

add_subdirectory(maths)
add_subdirectory(maths/fourier)
add_subdirectory (raytracing) #
add_subdirectory (examples/ot) #
add_subdirectory (examples/paths) #
add_subdirectory (examples/layers) #
add_subdirectory (examples/axicon) #
# add_subdirectory (examples/mctest) # 
add_subdirectory (examples/usptest) #
add_subdirectory (examples/sphere) #
add_subdirectory (examples/kirchhoff3d) #
## add_subdirectory (examples/math)
add_subdirectory (examplesXML) #
add_subdirectory (examples/lenses) # 
add_subdirectory (XML) #
add_subdirectory (test)

