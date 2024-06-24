set (XML_HDRS xml.h "xmltoken.h")
set (XML_SRCS xml.cpp "xmltoken.cpp")

add_library (goat_xml ${XML_SRCS} ${XML_HDRS} )
target_compile_definitions (goat_xml PUBLIC WITH_OCTREE WITH_SUPERGITTER)
target_include_directories (goat_xml PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths ../tinyxml2 ../raytracing)
install(TARGETS goat_xml
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
                )

