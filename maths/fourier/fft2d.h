#pragma once

#include <complex>
#include <memory>
#include "fftutils.h"
#include "vector.h"

/**
 * @defgroup angularSpectrum Angular spectrum propagation
 * @brief Functions for forward and inverse transformations in the angular spectrum framework.
 *
 * These routines provide discrete transformations of complex-valued fields
 * between spatial domain and spectral domain, as well as scaling operations.
 *
 * The functions operate component-wise on complex fields. Depending on the
 * variant, results are either written to the output field or added to an
 * existing field.
 */

namespace GOAT
{
    namespace maths
    {
        /**
        * @brief The fft2D class provides an interface to compute the 2D Fast Fourier Transform (FFT) and its inverse.
         * It uses the FFTW library for efficient computation. The class manages the necessary buffers and plans for the FFT operations.
         * 
         * @note The class is non-copyable but movable. It provides methods for forward and inverse FFT, as well as scaling the inverse transform.
         * 
		 * @see https://www.fftw.org/ for more information on the FFTW library.
        */
        namespace fourier
        {
      
            class fft2D
            {
            public:
		
                fft2D(std::size_t nx, std::size_t ny);
                ~fft2D();

                fft2D(const fft2D&) = delete;
                fft2D& operator=(const fft2D&) = delete;

                fft2D(fft2D&& other) noexcept;
                fft2D& operator=(fft2D&& other) noexcept;

				std::size_t nx() const noexcept; ///< returns the size of the first dimension (x) of the FFT
				std::size_t ny() const noexcept; ///< returns the size of the second dimension (y) of the FFT
				std::size_t size() const noexcept; ///< returns the total size of the FFT (nx * ny)


                /**
          * @brief Performs forward transformation of a complex field.
          * @ingroup angularSpectrum
          *
          * Transforms a complex field from spatial domain to spectral domain.
          * The output field is overwritten.
          *
          * This operation typically corresponds to a discrete Fourier transform.
          *
          * @param in Input field in spatial domain.
          * @param out Output field in spectral domain.
          */
                void forward(const std::complex<double>* in, std::complex<double>* out);
                void forward(const vectorField2D& in, vectorField2D& out, fieldComponent component);
				void forward(const vectorField2D& in, fftw_complex* out, fieldComponent component);
                void forward(const vectorField2D& in, vectorField2D& out);
                void forward(const complexd* in, fftw_complex* out);
                void forward(fftw_complex* in, fftw_complex* out);

                /**
         * @brief Performs inverse transformation of a complex field.
         * @ingroup angularSpectrum
         *
         * Transforms a complex field from spectral domain back to spatial domain.
         * The output field is overwritten.
         *
         * This operation typically corresponds to an inverse discrete Fourier transform.
         *
         * @param in Input field in spectral domain.
         * @param out Output field in spatial domain.
         */
                void inverse(const std::complex<double>* in, std::complex<double>* out);
                void inverse(const vectorField2D& in, vectorField2D& out, fieldComponent component);
				void inverse(fftw_complex* in, vectorField2D& out, fieldComponent component);
                void inverse(const vectorField2D& in, vectorField2D& out);
                void inverse(const fftw_complex* in, complexd* out);
                void inverse(fftw_complex* in, fftw_complex* out);

              
            private:
                void scaleInverse(std::complex<double>* data) const;
                static void scaleInverse(std::complex<double>* data, std::size_t n);
                void scaleInverse(vectorField2D& data, fieldComponent component) const;
                void scaleInverse(vectorField2D& data) const;
                class impl;
                std::unique_ptr<impl> pImpl;
                fftw_plan m_forwardPlan = nullptr;

            };
        }
    }
}

