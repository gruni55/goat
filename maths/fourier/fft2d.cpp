#include "fft2d.h"

#include <stdexcept>
#include <vector>
#include <utility>
#include "fftutils.h"

namespace GOAT
{
    namespace maths
    {
        namespace fourier
        {
            class fft2D::impl
            {
            public:
                impl(std::size_t nx_, std::size_t ny_)
                    : m_nx(nx_), m_ny(ny_)
                {
                    if (m_nx == 0 || m_ny == 0)
                        throw std::invalid_argument("fft2D: nx and ny must be > 0");

                    const std::size_t n = m_nx * m_ny;

                    m_inBuffer = reinterpret_cast<fftw_complex*>(
                        fftw_malloc(sizeof(fftw_complex) * n));
                    m_outBuffer = reinterpret_cast<fftw_complex*>(
                        fftw_malloc(sizeof(fftw_complex) * n));

                    if (m_inBuffer == nullptr || m_outBuffer == nullptr)
                    {
                        cleanup();
                        throw std::runtime_error("fft2D: failed to allocate FFTW buffers");
                    }

                    m_forwardPlan = fftw_plan_dft_2d(
                        static_cast<int>(m_ny),
                        static_cast<int>(m_nx),
                        m_inBuffer,
                        m_outBuffer,
                        FFTW_FORWARD,
                        FFTW_MEASURE);

                    m_inversePlan = fftw_plan_dft_2d(
                        static_cast<int>(m_ny),
                        static_cast<int>(m_nx),
                        m_inBuffer,
                        m_outBuffer,
                        FFTW_BACKWARD,
                        FFTW_MEASURE);

                    if (m_forwardPlan == nullptr || m_inversePlan == nullptr)
                    {
                        cleanup();
                        throw std::runtime_error("fft2D: failed to create FFTW plans");
                    }
                }

                ~impl()
                {
                    cleanup();
                }


				//  ---------------- forward methods ------------------

                void forward(const std::complex<double>* in, std::complex<double>* out)
                {
                    if (in == nullptr || out == nullptr)
                        throw std::invalid_argument("fft2D::forward: null pointer");

                    const std::size_t n = m_nx * m_ny;
                    copy(in, m_inBuffer, n);
                    fftw_execute(m_forwardPlan);
                    copy(m_outBuffer, out, n);                
                }

                void forward(const vectorField2D& in, vectorField2D& out, fieldComponent component)
                {
                    copy(in, m_inBuffer, m_nx, m_ny, component);
                    fftw_execute(m_forwardPlan);
                    copy(m_outBuffer, out, m_nx, m_ny, component);
                }

                void forward(const complexd* in, fftw_complex* out)
                {
                    if (in == nullptr || out == nullptr)
                        throw std::invalid_argument("fft2D::forward: null pointer");
                    const std::size_t n = m_nx * m_ny;
                    copy(in, m_inBuffer, n);
                    fftw_execute_dft(m_forwardPlan,m_inBuffer, out);
				}

                void forward(const vectorField2D& in, fftw_complex* out, fieldComponent component)
                {
					if (out == nullptr)
                        throw std::invalid_argument("fft2D::forward: null pointer");
                    
                    if (m_inBuffer == nullptr)
                        throw std::invalid_argument("fft2D::forward: m_inBuffer: null pointer");
                    
                    copy(in, m_inBuffer, m_nx, m_ny, component);
                    fftw_execute_dft(m_forwardPlan, m_inBuffer, out); 
                }

                void forward(fftw_complex* in, fftw_complex* out)
                {
                    if (in == nullptr || out == nullptr)
                        throw std::invalid_argument("fft2D::forward: null pointer");

                    fftw_execute_dft(m_forwardPlan, in, out);
                }



                void inverse(fftw_complex* in, fftw_complex* out)
                {
                    if (in == nullptr || out == nullptr)
                        throw std::invalid_argument("fft2D::inverse: null pointer");

                    fftw_execute_dft(m_inversePlan, in, out);
                }

                void inverse(const fftw_complex* in, complexd* out)
                {
                    if (in == nullptr || out == nullptr)
                        throw std::invalid_argument("fft2D::inverse: null pointer");
                    const std::size_t n = m_nx * m_ny;
                    fftw_execute_dft(m_inversePlan, const_cast<fftw_complex*>(in), m_outBuffer);
                    copy(m_outBuffer, out, n);
				}


                void inverse(const std::complex<double>* in, std::complex<double>* out)
                {
                    if (in == nullptr || out == nullptr)
                        throw std::invalid_argument("fft2D::inverse: null pointer");

                    const std::size_t n = m_nx * m_ny;
                    copy(in, m_inBuffer, n);
                    fftw_execute(m_inversePlan);
                    copy(m_outBuffer, out, n);
                }

                void inverse(const vectorField2D& in, vectorField2D& out, fieldComponent component)
                {
                    copy(in, m_inBuffer, m_nx, m_ny, component);
                    fftw_execute(m_inversePlan);
                    copy(m_outBuffer, out, m_nx, m_ny, component);
                }

                void inverse(fftw_complex* in,vectorField2D& out, fieldComponent component)
                {
                     if (in == nullptr) throw std::invalid_argument("fft2D::inverse: null pointer");
                    fftw_execute_dft(m_inversePlan,in,m_outBuffer);
					copy(m_outBuffer, out, m_nx, m_ny, component); }
                void inverse(const vectorField2D& in, vectorField2D& out)
                {
                    inverse(in, out, fieldComponent::Ex);
                    inverse(in, out, fieldComponent::Ey);
                    inverse(in, out, fieldComponent::Ez);
                }


                std::size_t nx() const noexcept
                {
                    return m_nx;
                }

                std::size_t ny() const noexcept
                {
                    return m_ny;
                }

                std::size_t size() const noexcept
                {
                    return m_nx * m_ny;
                }

            private:
                void cleanup() noexcept
                {
                    if (m_forwardPlan != nullptr)
                    {
                        fftw_destroy_plan(m_forwardPlan);
                        m_forwardPlan = nullptr;
                    }

                    if (m_inversePlan != nullptr)
                    {
                        fftw_destroy_plan(m_inversePlan);
                        m_inversePlan = nullptr;
                    }

                    if (m_inBuffer != nullptr)
                    {
                        fftw_free(m_inBuffer);
                        m_inBuffer = nullptr;
                    }

                    if (m_outBuffer != nullptr)
                    {
                        fftw_free(m_outBuffer);
                        m_outBuffer = nullptr;
                    }
                }

            private:
                std::size_t m_nx = 0;
                std::size_t m_ny = 0;
                fftw_complex* m_inBuffer = nullptr;
                fftw_complex* m_outBuffer = nullptr;

                fftw_plan m_forwardPlan = nullptr;
                fftw_plan m_inversePlan = nullptr;
            };

            fft2D::fft2D(std::size_t nx, std::size_t ny)
                : pImpl(std::make_unique<impl>(nx, ny))
            {}

          
            fft2D::~fft2D() = default;

            fft2D::fft2D(fft2D&& other) noexcept = default;
            fft2D& fft2D::operator=(fft2D&& other) noexcept = default;

            std::size_t fft2D::nx() const noexcept
            {
                return pImpl->nx();
            }

            std::size_t fft2D::ny() const noexcept
            {
                return pImpl->ny();
            }

            std::size_t fft2D::size() const noexcept
            {
                return pImpl->size();
            }


			// ------------------ forward methods ------------------
            void fft2D::forward(const vectorField2D& in, vectorField2D& out)
            {
                forward(in, out, fieldComponent::Ex);
                forward(in, out, fieldComponent::Ey);
                forward(in, out, fieldComponent::Ez);
            }

            void fft2D::forward(const std::complex<double>* in, std::complex<double>* out)
            {
                pImpl->forward(in, out);
            }

            void fft2D::forward(const vectorField2D& in, vectorField2D& out, fieldComponent component)
            {
                pImpl->forward(in, out, component);
            }

            void fft2D::forward(const vectorField2D& in, fftw_complex* out, fieldComponent component)
            {
                pImpl->forward(in, out, component);
            }

            void fft2D::forward(const complexd* in, fftw_complex* out)
            {
				pImpl->forward(in, out);
            }

            void fft2D::forward(fftw_complex* in, fftw_complex* out)
            {
				pImpl->forward(in, out);
            }


			// ------------------ inverse methods ------------------

            void fft2D::inverse(fftw_complex* in, fftw_complex* out)
            {
                pImpl->inverse(in, out);

                const double scale =
                    1.0 / static_cast<double>(pImpl->size());

                for (std::size_t k = 0; k < pImpl->size(); ++k)
                {
                    out[k][0] *= scale;
                    out[k][1] *= scale;
                }
            }

            void fft2D::inverse(const std::complex<double>* in, std::complex<double>* out)
            {
                pImpl->inverse(in, out);
				scaleInverse(out);
            }

            void fft2D::inverse(const vectorField2D& in, vectorField2D& out, fieldComponent component)
            {
                pImpl->inverse(in, out, component);
				scaleInverse(out, component);
            }

            void fft2D::inverse(fftw_complex* in, vectorField2D& out, fieldComponent component)
            {
				pImpl->inverse(in, out, component);
				scaleInverse(out, component);
            }

            void fft2D::inverse(const vectorField2D& in, vectorField2D& out)
            {
                pImpl->inverse(in, out);
				scaleInverse(out);
            }

            void fft2D::inverse(const fftw_complex* in, complexd* out)
            {
				pImpl->inverse(in, out);
            }


            void fft2D::scaleInverse(std::complex<double>* data) const
            {
                scaleInverse(data, size());
            }

            void fft2D::scaleInverse(std::complex<double>* data, std::size_t n)
            {
                if (data == nullptr)
                    throw std::invalid_argument("fft2D::scaleInverse: null pointer");

                if (n == 0)
                    return;

                const double s = 1.0 / static_cast<double>(n);

                for (std::size_t i = 0; i < n; ++i)
                {
                    data[i] *= s;
                }
            }
            void fft2D::scaleInverse(vectorField2D& data, fieldComponent component) const
            {
                if (data.empty())
                    return;
                
                const double s = 1.0 / static_cast<double>(size());
                for (size_t ix = 0; ix < data.size(); ++ix)
                    for (size_t iy = 0; iy < data[ix].size(); iy++)
                    {
                        data[ix][iy][static_cast<size_t>(component)] *= s;
                    }
            }
            void fft2D::scaleInverse(vectorField2D& data) const
            {
                if (data.empty() )
                    return;
                
                const double s = 1.0 / static_cast<double>(size());
                for (size_t ix = 0; ix < data.size(); ++ix)
                    for (size_t iy = 0; iy < data[ix].size(); iy++)
                            data[ix][iy]*= s;
						
            }
		} // namespace fourier
	} // namespace maths
} // namespace GOAT
    
