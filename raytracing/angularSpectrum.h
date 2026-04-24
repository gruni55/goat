#pragma once
#include "propagator.h"
#include "fourier/fft2d.h"
#include "fourier/fftutils.h"
namespace GOAT
{
	namespace raytracing
	{
		class AngularSpectrum : public Propagator
		{
			public:
			AngularSpectrum(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
			AngularSpectrum(double wvl, maths::Vector<double>P, maths::Vector<double> n, double d, int N);
			~AngularSpectrum() 
			{
				cleanup();
			}
			
			void calcOne(DetectorPlane* det, bool clear);

		private:
			GOAT::maths::fourier::fft2D fft;
			double dx = 0.0;
			double dy = 0.0;

			void propagateComponent(DetectorPlane* det, maths::fourier::fieldComponent component, double dz);
			void AngularSpectrum::applyTransferFunction(fftw_complex* spec, double dz);
			void addField(const maths::fourier::vectorField2D& field, maths::fourier::fieldComponent component);
			fftw_complex* spec;
			void initSpec()
			{
				const std::size_t nx = this->N1();
				const std::size_t ny = this->N2();
				spec = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nx * ny));
				if (spec == nullptr)
					throw std::runtime_error("AngularSpectrum: failed to allocate spec buffer");
			}
			void cleanup()
			{
				if (spec != nullptr)
				{
					fftw_free(spec);
					spec = nullptr;
				}
			}
		};
	}
}