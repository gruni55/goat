#include "angularspectrum.h"
namespace GOAT
{
	namespace raytracing
	{
		AngularSpectrum::AngularSpectrum(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2)
			: Propagator(wvl, P, e1, e2, n1, n2), fft(static_cast<std::size_t>(n1), static_cast<std::size_t>(n2))
		{
			type = DETECTOR_ANGULAR_SPECTRUM;
			dx = abs(e1) / static_cast<double>(n1);
			dy = abs(e2) / static_cast<double>(n2);
			initSpec();
		}

		AngularSpectrum::AngularSpectrum(double wvl, maths::Vector<double> P, maths::Vector<double> n, double d, int N) 
			: Propagator(wvl, P, n, d, N), fft(static_cast<std::size_t>(n1), static_cast<std::size_t>(n2))
		{
			type = DETECTOR_ANGULAR_SPECTRUM;
			dx = abs(e1) / static_cast<double>(n1);
			dy = abs(e2) / static_cast<double>(n2);
			initSpec();
		}

		void AngularSpectrum::calcOne(DetectorPlane* det, bool clear)
		{
			if (det == nullptr)
				throw std::invalid_argument("AngularSpectrumPropagator::calcOne: det is null");

			// 1. Dimensionen prüfen
			if (det->N1() != this->N1() || det->N2() != this->N2())
				throw std::runtime_error("AngularSpectrumPropagator::calcOne: source and target grid size differ");

			// 2. Abstand bestimmen
			maths::Vector<double> dP = P - det->position();
			maths::Vector<double> n = this->norm();   // oder det->getNormal(), wenn parallel
			double dz = dP * n;

			// 3. Ebenenparallelität prüfen
			maths::Vector<double> nSrc = det->norm();
			maths::Vector<double> nDst = n;

			if ((abs(nSrc * nDst) - 1.0) > 1e-8)
				throw std::runtime_error("AngularSpectrumPropagator::calcOne: source and target plane are not parallel");

			// 4. Für jede Komponente propagieren
			propagateComponent(det, maths::fourier::fieldComponent::Ex, dz);
			propagateComponent(det, maths::fourier::fieldComponent::Ey, dz);
			propagateComponent(det, maths::fourier::fieldComponent::Ez, dz);
		}



		void AngularSpectrum::propagateComponent(DetectorPlane* det, maths::fourier::fieldComponent component, double dz)
		{
			const std::size_t nx = this->N1();
			const std::size_t ny = this->N2();

			maths::fourier::vectorField2D tmp;

			// 1. Vorwärts-FFT: Quelldetektor -> Frequenzraum
			fft.forward(det->D, spec, component);

			// 2. Transferfunktion anwenden
			applyTransferFunction(spec, dz);

			// 3. Rücktransformation
			fft.inverse(spec, tmp, component);

			// 4. Ergebnis auf Ziel addieren
			addField(tmp, component);
		}

		void GOAT::raytracing::AngularSpectrum::addField(const maths::fourier::vectorField2D& field, maths::fourier::fieldComponent component)
		{
			for (std::size_t x = 0; x < field.size(); ++x)
				for (std::size_t y = 0; y < field[x].size(); ++y)
					D[x][y][static_cast<size_t>(component)] += field[x][y][static_cast<size_t>(component)];
		}

		int shiftedIndex(std::size_t i, std::size_t N)
		{
			return (i < N / 2)
				? static_cast<int>(i)
				: static_cast<int>(i) - static_cast<int>(N);
		}

		void AngularSpectrum::applyTransferFunction(fftw_complex* spec, double dz)
		{
			if (spec == nullptr)
				throw std::invalid_argument("AngularSpectrum::applyTransferFunction: spec is null");

			const std::size_t nx = this->n1;   // oder getN1()
			const std::size_t ny = this->n2;   // oder getN2()

			const double k0 = 2.0 * M_PI / getWavelength();

			for (std::size_t y = 0; y < ny; ++y)
			{
				const int my = shiftedIndex(y, ny);
				const double ky = 2.0 * M_PI * static_cast<double>(my)
					/ (static_cast<double>(ny) * dy);

				for (std::size_t x = 0; x < nx; ++x)
				{
					const int mx = shiftedIndex(x, nx);
					const double kx = 2.0 * M_PI * static_cast<double>(mx)
						/ (static_cast<double>(nx) * dx);

					const double kz2 = k0 * k0 - kx * kx - ky * ky;
					const std::size_t idx = y * nx + x;

					if (kz2 < 0.0)
					{
						// evaneszente Anteile für den ersten Ansatz abschneiden
						spec[idx][0] = 0.0;
						spec[idx][1] = 0.0;
						continue;
					}

					const double kz = std::sqrt(kz2);
					const double phase = kz * dz;

					const double cr = std::cos(phase);
					const double ci = std::sin(phase);

					const double ar = spec[idx][0];
					const double ai = spec[idx][1];

					// spec *= exp(i*kz*dz)
					spec[idx][0] = ar * cr - ai * ci;
					spec[idx][1] = ar * ci + ai * cr;
				}
			}
		}

	} // raytracing
} // GOAT