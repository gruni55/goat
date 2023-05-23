set (XML_HDRS xml.h)
set (XML_SRCS xml.cpp)

add_library (goat_xml ${XML_SRCS} ${XML_HDRS} )
target_include_directories (goat_xml PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths ../tinyxml2 ../raytracing)
