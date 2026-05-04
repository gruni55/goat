if(TARGET goat::fftw)
    return()
endif()

set(FFTW3_ROOT "$ENV{FFTW3_ROOT}" CACHE PATH "Root directory of FFTW3")

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS "${FFTW3_ROOT}"
    PATH_SUFFIXES include ""
)

find_library(FFTW3_LIBRARY
    NAMES fftw3 fftw3-3 libfftw3-3
    HINTS "${FFTW3_ROOT}"
    PATH_SUFFIXES lib lib64 bin ""
)

find_file(FFTW3_DLL
    NAMES libfftw3-3.dll fftw3.dll fftw3-3.dll
    HINTS "${FFTW3_ROOT}"
    PATH_SUFFIXES bin lib ""
)

message(STATUS "FFTW3_ROOT        = ${FFTW3_ROOT}")
message(STATUS "FFTW3_INCLUDE_DIR = ${FFTW3_INCLUDE_DIR}")
message(STATUS "FFTW3_LIBRARY     = ${FFTW3_LIBRARY}")
message(STATUS "FFTW3_DLL         = ${FFTW3_DLL}")

if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
    message(FATAL_ERROR
        "FFTW3 not found.\n"
        "On Ubuntu install it with:\n"
        "  sudo apt install libfftw3-dev\n"
        "On Windows set FFTW3_ROOT to your FFTW directory."
    )
endif()

add_library(goat::fftw UNKNOWN IMPORTED)

set_target_properties(goat::fftw PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
)

if(WIN32 AND FFTW3_DLL)
    set_target_properties(goat::fftw PROPERTIES
        IMPORTED_IMPLIB "${FFTW3_LIBRARY}"
        IMPORTED_LOCATION "${FFTW3_DLL}"
    )
else()
    set_target_properties(goat::fftw PROPERTIES
        IMPORTED_LOCATION "${FFTW3_LIBRARY}"
    )
endif()
