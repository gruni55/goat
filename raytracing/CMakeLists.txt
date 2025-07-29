#CMAKE configuration file for ray library
message ("Preparing library goat_raytracing")
set (
	RAYTRACING_SRCS box.cpp detector.cpp ellipsoid.cpp iray.cpp  
                     lightsrc.cpp lightsrc_mc.cpp misc.cpp objectshape.cpp octree.cpp 
		     raytrace.cpp ray_pow.cpp resutil.cpp surface.cpp triangle.cpp tubedray.cpp 
		     triangle_box_intersection.cpp  superarray.cpp 
		     error.cpp mem.cpp raytrace_inel.cpp raytrace_usp.cpp inel_calc.cpp raytrace_field.cpp 
		     raytrace_usp.cpp raytrace_field_usp.cpp raytrace_usp_rt.cpp plane.cpp pulsecalculation_field.cpp pulsecalculation.cpp pulsecalculation_rt.cpp gridentry.cpp asphericLens.cpp
		     sphericLens.cpp lens.cpp fft.cpp  cone.cpp cylinder.cpp vortex.cpp refractive_index_functions.cpp)  

set (
	RAYTRACING_HDRS box.h detector.h ellipsoid.h iray.h 
                     lightsrc.h lightsrc_mc.h misc.h objectshape.h octree.h 
		     raybase.h raytrace.h ray_pow.h resutil.h surface.h triangle.h tubedray.h 
		     triangle_box_intersection.h  superarray.h supergrid.h
		     error.h mem.h raytrace_inel.h raytrace_usp.h raytrace_field.h inel_calc.h 
		     raytrace_usp_rt.h refractive_index_functions.h plane.h  fft.h gridentry.h pulsecalculation.h 
		     pulsecalculation_rt.h asphericLens.h sphericLens.h lens.h constants.h cone.h cylinder.h vortex.h)


if (WITH_OCTREE)
	message(STATUS "Octree calculation enabled")
else()
	message (STATUS "Calculation WITHOUT octree")
endif()

add_compile_definitions (WITH_SUPERGITTER)
# option(WITH_OPENMP "Enable OpenMP parallelization" OFF)
if (WITH_OPENMP) 
	find_package(OpenMP)
        if (OpenMP_FOUND)   
#	if (OPENMP_CXX_FOUND)
 		message (STATUS "OpenMP found ==> creating goat_raytracing_mp library")
                add_library (goat_raytracing_mp  ${RAYTRACING_SRCS} ${RAYTRACING_HDRS} )
	        target_compile_options (goat_raytracing_mp PUBLIC ${OpenMP_CXX_FLAGS})
                target_compile_definitions (goat_raytracing_mp PUBLIC WITH_OCTREE WITH_OPENMP)
				target_link_libraries(goat_raytracing_mp PUBLIC OpenMP::OpenMP_CXX) # <- WICHTIG!
                target_include_directories (goat_raytracing_mp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths ../ticpp)  
		install(TARGETS goat_raytracing_mp 
                RUNTIME DESTINATION lib/$<CONFIG>
                LIBRARY DESTINATION lib/$<CONFIG>
                ARCHIVE DESTINATION lib/$<CONFIG>
                )
	else()
		message (WARNING "OpenMP not found - Calculation without parallelization")
	endif()
else()
message (STATUS "OpenMP support skipped due to parameter \"WITHOUT_OPENMP\"")
endif()
message (STATUS "Write libraries to:  ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}" )
add_library (goat_raytracing ${RAYTRACING_SRCS} ${RAYTRACING_HDRS}  )
target_compile_definitions (goat_raytracing PUBLIC WITH_OCTREE)
target_include_directories (goat_raytracing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths ../ticpp)
install(TARGETS goat_maths goat_raytracing  
                RUNTIME DESTINATION lib/$<CONFIG>
                LIBRARY DESTINATION lib/$<CONFIG>
                ARCHIVE DESTINATION lib/$<CONFIG>
                )
