#CMAKE configuration file for ray library

set (RAYTRACING_SRCS box.cpp detector.cpp ellipsoid.cpp iray.cpp lightsrc.cpp lightsrc_mc.cpp objectshape.cpp octree.cpp raytrace.cpp ray_pow.cpp surface.cpp triangle.cpp tubedray.cpp triangle_box_intersection.cpp resutil.cpp misc.cpp superarray.cpp error.cpp mem.cpp raytrace_inel.cpp inel_calc.cpp plane.cpp grid.cpp)  

set (RAYTRACING_HDRS box.h detector.h ellipsoid.h iray.h lightsrc.h lightsrc_mc.h misc.h objectshape.h octree.h  raybase.h raytrace.h ray_pow.h resutil.h surface.h triangle.h tubedray.h triangle_box_intersection.h superarray.h error.h mem.h raytrace_inel.cpp inel_calc.h plane.h grid.h)

add_compile_definitions (WITH_SUPERGITTER)
add_library (goat_raytracing ${RAYTRACING_SRCS} ${RAYTRACING_HDRS} )
target_compile_definitions (goat_raytracing PUBLIC WITH_OCTREE)
target_include_directories (goat_raytracing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths )
