#include "angularspectrum.h"
#include <omp.h>
#include "fourier/fft2d.h"
#include "fourier/fftutils.h"
namespace GOAT
{
	namespace raytracing
	{

		int shiftedIndex(std::size_t i, std::size_t N)
		{
			return (i < N / 2)
				? static_cast<int>(i)
				: static_cast<int>(i) - static_cast<int>(N);
		}

		struct AngularSpectrum::Impl
		{
			fftw_complex* spec = nullptr;
			GOAT::maths::fourier::fft2D fft;
			AngularSpectrum& self;

			Impl(AngularSpectrum& self)
				: self(self), fft(self.n1, self.n2)
			{
				spec = fftw_alloc_complex(self.n1 * self.n2);
				initSpec();
			}

			~Impl()
			{
				cleanup();
			}

			void addField(const maths::fourier::vectorField2D& field, maths::fourier::fieldComponent component)
			{
				for (std::size_t x = 0; x < field.size(); ++x)
					for (std::size_t y = 0; y < field[x].size(); ++y)
						self.D[x][y][static_cast<size_t>(component)] += field[x][y][static_cast<size_t>(component)];
			}

			void applyTransferFunction(fftw_complex* spec, double dz)
			{
				if (spec == nullptr)
					throw std::invalid_argument("AngularSpectrum::applyTransferFunction: spec is null");

				const std::size_t nx = self.n1;   // oder getN1()
				const std::size_t ny = self.n2;   // oder getN2()

				const double k0 = 2.0 * M_PI / self.wvl;

				for (std::size_t y = 0; y < ny; ++y)
				{
					const int my = shiftedIndex(y, ny);
					const double ky = 2.0 * M_PI * static_cast<double>(my)
						/ (static_cast<double>(ny) * self.dy);

					for (std::size_t x = 0; x < nx; ++x)
					{
						const int mx = shiftedIndex(x, nx);
						const double kx = 2.0 * M_PI * static_cast<double>(mx)
							/ (static_cast<double>(nx) * self.dx);

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
					} // for x
				} // for y
			}
				void addRoi(const maths::fourier::vectorField2D & tmp, DetectorPlane * src, maths::fourier::fieldComponent component)
				{
					int ix0 = 0;
					int iy0 = 0;

					if (!self.computeRoiOffset(src, ix0, iy0))
						throw std::runtime_error("AngularSpectrum::addRoi: target ROI does not fit source grid");

					const std::size_t c = static_cast<std::size_t>(component);

					int n1 = static_cast<int>(self.N1());
					int n2 = static_cast<int>(self.N2());

					for (std::size_t xT = 0; xT < n1; ++xT)
					{
						for (std::size_t yT = 0; yT < n2; ++yT)
						{
							std::size_t xS = static_cast<std::size_t>(ix0) + xT;
							std::size_t yS = static_cast<std::size_t>(iy0) + yT;

							self.D[xT][yT][c] += tmp[xS][yS][c];
						}
					}
				}

			

			void propagateComponent(DetectorPlane* det, maths::fourier::fieldComponent component, double dz)
			{
				const std::size_t nx = static_cast<int>(self.n1);
				const std::size_t ny = static_cast<int>(self.n2);

				maths::fourier::vectorField2D tmp;

				// 1. Vorwärts-FFT: Quelldetektor -> Frequenzraum
				fft.forward(det->D, spec, component);

				// 2. Transferfunktion anwenden
				applyTransferFunction(spec, dz);

				// 3. Rücktransformation
				fft.inverse(spec, tmp, component);

				// 4. Ergebnis auf Ziel addieren
				// addField(tmp, component);
				addRoi(tmp, det, component);
			}

			void initSpec()
			{
				const std::size_t nx = static_cast<int>(self.n1);
				const std::size_t ny = static_cast<int>(self.n2);
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

		AngularSpectrum::AngularSpectrum(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2)
			: Propagator(wvl, P, e1, e2, n1, n2), impl(std::make_unique<Impl>(*this))
		{
			type = DETECTOR_ANGULAR_SPECTRUM;
			dx = d1 / static_cast<double>(n1);
			dy = d2 / static_cast<double>(n2);
		}

		AngularSpectrum::AngularSpectrum(double wvl, maths::Vector<double> P, maths::Vector<double> n, double d, int N) 
			: Propagator(wvl, P, n, d, N), impl(std::make_unique<Impl>(*this))
		{
			type = DETECTOR_ANGULAR_SPECTRUM;
			dx = d1 / static_cast<double>(n1);
			dy = d2 / static_cast<double>(n2);
		}

		void AngularSpectrum::calc(bool clear)
		{
			if (clear)
        clean();

    if (numberOfSources() == 0)
        throw std::runtime_error("AngularSpectrum::calc: no source detector available");

			calcOne(getSources().front(), false);
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
			maths::fourier::fieldComponent comp;
#pragma omp parallel for private(comp)
			for (int c = 0; c < 3; ++c)
			{
				comp = static_cast<maths::fourier::fieldComponent>(c);
				impl->propagateComponent(det, comp, dz);
			}
		}

		void AngularSpectrum::adaptTargetGridToSource(DetectorPlane* src)
		{
			const double srcDx = abs(src->gete1()) / static_cast<double>(src->N1());
			const double srcDy = abs(src->gete2()) / static_cast<double>(src->N2());

			dx = srcDx;
			dy = srcDy;

			maths::Vector<double> e1Hat = this->e1 / abs(this->e1);
			maths::Vector<double> e2Hat = this->e2 / abs(this->e2);

			double LxDesired = abs(this->e1);
			double LyDesired = abs(this->e2);

			int newN1 = static_cast<int>(std::round(LxDesired / dx));
			int newN2 = static_cast<int>(std::round(LyDesired / dy));

			if (newN1 <= 0 || newN2 <= 0)
				throw std::runtime_error("AngularSpectrum::adaptTargetGridToSource: invalid target grid size");

			this->n1 = newN1;
			this->n2 = newN2;

			this->e1 = static_cast<double>(newN1) * dx * e1Hat;
			this->e2 = static_cast<double>(newN2) * dy * e2Hat;

			clean();
		}

		bool AngularSpectrum::computeRoiOffset(DetectorPlane* src,
			int& ix0,
			int& iy0) const
		{
			maths::Vector<double> e1Hat = src->gete1() / abs(src->gete1());
			maths::Vector<double> e2Hat = src->gete2() / abs(src->gete2());

			maths::Vector<double> srcCorner = src->position() - 0.5 * src->gete1() - 0.5 * src->gete2();

			maths::Vector<double> targetCorner = P - 0.5 * this->e1 - 0.5 * this->e2;

			maths::Vector<double> shift = targetCorner - srcCorner;

			double sx = shift * e1Hat;
			double sy = shift * e2Hat;

			double fx = sx / dx;
			double fy = sy / dy;

			ix0 = static_cast<int>(std::round(fx));
			iy0 = static_cast<int>(std::round(fy));

			const double eps = 1e-6;

			if (std::abs(fx - static_cast<double>(ix0)) > eps)
				return false;

			if (std::abs(fy - static_cast<double>(iy0)) > eps)
				return false;

			if (ix0 < 0 || iy0 < 0)
				return false;

			if (ix0 + static_cast<int>(n1) > static_cast<int>(src->N1()))
				return false;

			if (iy0 + static_cast<int>(n2) > static_cast<int>(src->N2()))
				return false;

			return true;
		}

		AngularSpectrum::~AngularSpectrum() = default;

	} // raytracing
} // GOAT