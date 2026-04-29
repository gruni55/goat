#pragma once

#include <complex>
#include <vector>
#include <stdexcept>
#include <cstddef>

#include "fftw3.h"
#include "vector.h"

namespace GOAT
{
    namespace maths
    {
        namespace fourier
        {

            enum class fieldComponent
            {
                Ex = 0,
                Ey = 1,
                Ez = 2
            };

            using complexd = std::complex<double>;
            using complexVector = GOAT::maths::Vector<complexd>;
            using vectorField2D = std::vector<std::vector<complexVector>>;

            void copy(const vectorField2D& src,
                fftw_complex* dst,
                std::size_t nx,
                std::size_t ny,
                fieldComponent component);

            void copy(const fftw_complex* src,
                vectorField2D& dst,
                std::size_t nx,
                std::size_t ny,
                fieldComponent component);

            void copy(
                const fftw_complex* src,
                std::complex<double>* dst,
                std::size_t n);

            void copy(const std::complex<double>* src,
                fftw_complex* dst,
                std::size_t n);
   
        }
    }
}

