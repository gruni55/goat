set (XML_HDRS xml.h xmltoken.h tinyxml2.h storageBackend.h)
set (XML_SRCS xml.cpp xmltoken.cpp tinyxml2.cpp storageBackend.cpp)

add_library (goat_xml ${XML_SRCS} ${XML_HDRS} )

#if (MSVC)  # Microsoft Visual C++ Compiler (Windows)
#    target_compile_options(goat_xml PRIVATE
#        /fp:precise   # Floating Point Verhalten exakt, keine unsicheren Optimierungen
#        /O2           # Optimierungslevel 2
#    )
#else()  # GCC, Clang (Linux, macOS)
#    target_compile_options(goat_xml PRIVATE
# #       -O2
#        -fno-fast-math
#        -fno-unsafe-math-optimizations
#        -ffp-contract=off
#    )
#endif()

target_compile_definitions (goat_xml PUBLIC WITH_OCTREE WITH_SUPERGITTER)
target_include_directories (goat_xml PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../maths ../tinyxml2 ../raytracing)
install(TARGETS goat_xml
                RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/${CMAKE_BUILD_TYPE}
                LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}
                )
