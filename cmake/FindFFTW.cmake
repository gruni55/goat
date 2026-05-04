# - Find FFTW
# This module defines:
#   FFTW_FOUND
#   FFTW_INCLUDE_DIR
#   FFTW_LIB
#   FFTW_OMP_LIB (optional)
#
# It also defines imported targets:
#   FFTW::fftw3
#   FFTW::fftw3_omp (optional)

# ------------------------------------------------------------
# Include dir
# ------------------------------------------------------------

message(STATUS "Using GOAT FindFFTW.cmake")
find_path(FFTW_INCLUDE_DIR
    NAMES fftw3.h
    PATHS
        $ENV{FFTW3_ROOT}/include
        /usr/include
        /usr/local/include
        /opt/homebrew/include       # macOS (Apple Silicon)
        /usr/local/opt/fftw/include # macOS (Homebrew Intel)
)

# ------------------------------------------------------------
# Platform-specific library search
# ------------------------------------------------------------
if(WIN32)

    find_library(FFTW_LIB
        NAMES libfftw3-3 fftw3
        PATHS
            $ENV{FFTW3_ROOT}/lib
            $ENV{FFTW3_ROOT}
    )

    find_library(FFTW_OMP_LIB
        NAMES libfftw3_omp-3 fftw3_omp
        PATHS
            $ENV{FFTW3_ROOT}/lib
            $ENV{FFTW3_ROOT}
    )

elseif(APPLE)

    find_library(FFTW_LIB
        NAMES fftw3
        PATHS
            /opt/homebrew/lib
            /usr/local/lib
            /usr/local/opt/fftw/lib
            $ENV{FFTW3_ROOT}/lib
    )

    find_library(FFTW_OMP_LIB
        NAMES fftw3_omp
        PATHS
            /opt/homebrew/lib
            /usr/local/lib
            /usr/local/opt/fftw/lib
            $ENV{FFTW3_ROOT}/lib
    )

else() # Linux

    find_library(FFTW_LIB
        NAMES fftw3
        PATHS
            /usr/lib/x86_64-linux-gnu
            /usr/lib
            /usr/local/lib
            $ENV{FFTW3_ROOT}/lib
    )

    find_library(FFTW_OMP_LIB
        NAMES fftw3_omp
        PATHS
            /usr/lib/x86_64-linux-gnu
            /usr/lib
            /usr/local/lib
            $ENV{FFTW3_ROOT}/lib
    )

endif()

# ------------------------------------------------------------
# Handle result
# ------------------------------------------------------------
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_INCLUDE_DIR FFTW_LIB
)

# ------------------------------------------------------------
# Create imported targets
# ------------------------------------------------------------
if(FFTW_FOUND)

    if(NOT TARGET FFTW::fftw3)
        add_library(FFTW::fftw3 UNKNOWN IMPORTED)
        set_target_properties(FFTW::fftw3 PROPERTIES
            IMPORTED_LOCATION "${FFTW_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
        )
    endif()

    if(FFTW_OMP_LIB AND NOT TARGET FFTW::fftw3_omp)
        add_library(FFTW::fftw3_omp UNKNOWN IMPORTED)
        set_target_properties(FFTW::fftw3_omp PROPERTIES
            IMPORTED_LOCATION "${FFTW_OMP_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
        )
    endif()

endif()

# ------------------------------------------------------------
# Debug output (optional)
# ------------------------------------------------------------
message(STATUS "FFTW_INCLUDE_DIR: ${FFTW_INCLUDE_DIR}")
message(STATUS "FFTW_LIB: ${FFTW_LIB}")
message(STATUS "FFTW_OMP_LIB: ${FFTW_OMP_LIB}")
