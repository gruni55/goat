# Fourier submodule


set(FOURIER_HDRS
    fft2d.h fftutils.h
)

set(FOURIER_SRCS
    fft2d.cpp fftutils.cpp
)

add_library(goat_fourier ${FOURIER_SRCS} ${FOURIER_HDRS})

target_include_directories(goat_fourier PUBLIC $ENV{FFTW3_ROOT} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_directories(goat_fourier PUBLIC $ENV{FFTW3_ROOT})
target_link_libraries(goat_fourier
    PUBLIC goat_maths
    FFTW::fftw3
#    PRIVATE goat::fftw
)

target_compile_definitions(goat_fourier PRIVATE WITH_FFTW)

goat_install(TARGETS goat_fourier)


# =========================
# TESTPROGRAMM 
# =========================

list(APPEND CMAKE_MODULE_PATH "$ENV{GOATDir}/cmake")

find_package (FFTW REQUIRED)


add_executable(testFFTW fftwtest.cpp)

# target_link_libraries(testFFTW PRIVATE goat::fftw)

if(WIN32)
    add_custom_command(TARGET testFFTW POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$ENV{FFTW3_ROOT}/libfftw3-3.dll"
            $<TARGET_FILE_DIR:testFFTW>
    )
endif()
target_link_libraries(testFFTW PRIVATE goat_fourier)
target_include_directories(testFFTW PUBLIC "$ENV{FFTW3_ROOT}")
target_link_directories(testFFTW PUBLIC $ENV{FFTW3_ROOT})


#

add_executable(test_fft2d testfft2d.cpp)

target_link_libraries(test_fft2d PRIVATE
    goat_fourier
)

# DLL kopieren (Windows)
if(WIN32)
    add_custom_command(TARGET test_fft2d POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$ENV{FFTW3_ROOT}/libfftw3-3.dll"
            $<TARGET_FILE_DIR:test_fft2d>
    )
endif()

