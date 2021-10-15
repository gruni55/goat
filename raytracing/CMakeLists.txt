#CMAKE configuration file for ray library

set (RAYTRACING_SRCS box.cpp detector.cpp ellipsoid.cpp iray.cpp lightsrc.cpp objectshape.cpp octree.cpp plane.cpp Raytrace.cpp ray_pow.cpp surface.cpp triangle.cpp tubedray.cpp triangle_box_intersection.cpp resutil.cpp misc.cpp)

set (RAYTRACING_HDRS box.h Detector.h ellipsoid.h iray.h LightSrc.h misc.h objectshape.h octree.h  plane.h raybase.h raytrace.h ray_pow.h resutil.h surface.h triangle.h tubedray.h triangle_box_intersection.h)

add_library (goat_raytracing ${RAYTRACING_SRCS} ${RAYTRACING_HDRS} )
target_compile_definitions (goat_raytracing PUBLIC WITH_OCTREE)
target_include_directories (goat_raytracing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths )
