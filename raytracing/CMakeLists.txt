#CMAKE configuration file for ray library
message ("Preparing library goat_raytracing")
set (
	RAYTRACING_SRCS box.cpp detector.cpp ellipsoid.cpp iray.cpp  
                     lightsrc.cpp lightsrc_mc.cpp misc.cpp objectshape.cpp octree.cpp 
		     raytrace.cpp ray_pow.cpp resutil.cpp surface.cpp triangle.cpp tubedray.cpp 
		     triangle_box_intersection.cpp  superarray.cpp 
		     error.cpp mem.cpp raytrace_inel.cpp raytrace_usp.cpp inel_calc.cpp 
		     plane.cpp pulsecalculation.cpp
     "fft.cpp")  

set (
	RAYTRACING_HDRS box.h detector.h ellipsoid.h iray.h 
                     lightsrc.h lightsrc_mc.h misc.h objectshape.h octree.h 
		     raybase.h raytrace.h ray_pow.h resutil.h surface.h triangle.h tubedray.h 
		     triangle_box_intersection.h  superarray.h supergrid.h
		     error.h mem.h raytrace_inel.h raytrace_usp.h inel_calc.h 
		     refractive_index_functions.h plane.h  fft.h gridentry.h pulsecalculation.h
     "constants.h")

add_compile_definitions (WITH_SUPERGITTER)
add_library (goat_raytracing ${RAYTRACING_SRCS} ${RAYTRACING_HDRS} )

if (WITHOUT_OPENMP) 
	message (STATUS "OpenMP support skipped due to parameter \"WITHOUT_OPENMP\"")
else()
	find_package(OpenMP)
	if (OPENMP_FOUND)
 		message (STATUS "OpenMP found")
	        message ("OpenMP compiler flags:" ${OpenMP_CXX_FLAGS})
	        target_compile_options (goat_raytracing PUBLIC ${OpenMP_CXX_FLAGS})
	else()
		message (WARNING "OpenMP not found - Calculation without parallelization")
	endif()
endif()
target_compile_definitions (goat_raytracing PUBLIC WITH_OCTREE)
target_include_directories (goat_raytracing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths ../ticpp)
