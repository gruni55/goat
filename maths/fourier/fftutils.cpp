#include "fftutils.h"
#include <cmath>

namespace GOAT
{
    namespace maths
    {
        namespace fourier
        {


            void copy(const std::complex<double>* src,
                fftw_complex* dst,
                std::size_t n)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    dst[i][0] = src[i].real();
                    dst[i][1] = src[i].imag();
                }
            }

            void copy(
                const fftw_complex* src,
                std::complex<double>* dst,
                std::size_t n)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    dst[i] = std::complex<double>(src[i][0], src[i][1]);
                }
            }



            void copy(const vectorField2D& src,
                fftw_complex* dst,
                std::size_t nx,
                std::size_t ny,
                fieldComponent component)
            {
                if (dst == nullptr)
                    throw std::invalid_argument("copyToFftwBuffer: dst is null");

                if (src.size() != nx)
                    throw std::runtime_error("copyToFftwBuffer: wrong x dimension");

                for (std::size_t x = 0; x < nx; ++x)
                {
                    if (src[x].size() != ny)
                        throw std::runtime_error("copyToFftwBuffer: wrong y dimension");
					int comp = static_cast<int>(component);
                    for (std::size_t y = 0; y < ny; ++y)
                    {
                        const std::size_t idx = y * nx + x;
                        const complexd& v = src[x][y][comp];
                        dst[idx][0] = v.real();
                        dst[idx][1] = v.imag();
                    }
                }
            }

            void copy(const fftw_complex* src,
                vectorField2D& dst,
                std::size_t nx,
                std::size_t ny,
                fieldComponent component)
            {
                if (src == nullptr)
                    throw std::invalid_argument("copyFromFftwBuffer: src is null");


                if (dst.size() != nx)
                    dst.resize(nx);

				int comp = static_cast<int>(component);
                for (std::size_t x = 0; x < nx; ++x)
                {
                    if (dst[x].size() != ny)
                        dst[x].resize(ny);

                    for (std::size_t y = 0; y < ny; ++y)
                    {
                        const std::size_t idx = y * nx + x;
                        dst[x][y][comp] = complexd(src[idx][0], src[idx][1]);
                    }
                }
            }



        } // namespace fourier
    } //    namespace maths
}// namespace GOAT
