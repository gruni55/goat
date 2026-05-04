if(TARGET goat::fftw)
    return()
endif()

set(FFTW3_ROOT $ENV{FFTW3_ROOT})
if(NOT FFTW3_ROOT)
    set(FFTW3_ROOT "" CACHE PATH "Root directory of FFTW3")
endif()

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS ${FFTW3_ROOT}
    PATH_SUFFIXES "" include
)

find_library(FFTW3_IMPLIB
    NAMES fftw3 fftw3-3 libfftw3-3
    HINTS ${FFTW3_ROOT}
    PATH_SUFFIXES "" lib
)

find_file(FFTW3_DLL
    NAMES libfftw3-3.dll fftw3.dll fftw3-3.dll
    HINTS ${FFTW3_ROOT}
    PATH_SUFFIXES "" bin lib
)

message(STATUS "FFTW3_ROOT        = ${FFTW3_ROOT}")
message(STATUS "FFTW3_INCLUDE_DIR = ${FFTW3_INCLUDE_DIR}")
message(STATUS "FFTW3_IMPLIB      = ${FFTW3_IMPLIB}")
message(STATUS "FFTW3_DLL         = ${FFTW3_DLL}")

if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_IMPLIB)
    message(FATAL_ERROR "FFTW3 for Windows not found")
endif()

add_library(goat::fftw SHARED IMPORTED)

if(FFTW3_DLL)
    set_target_properties(goat::fftw PROPERTIES
        IMPORTED_IMPLIB "${FFTW3_IMPLIB}"
        IMPORTED_LOCATION "${FFTW3_DLL}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    )
else()
    set_target_properties(goat::fftw PROPERTIES
        IMPORTED_IMPLIB "${FFTW3_IMPLIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    )
endif()