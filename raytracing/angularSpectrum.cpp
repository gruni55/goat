#include "angularSpectrum.h"
#include <omp.h>
#include <array>
#include <algorithm>
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
			
			struct Workspace
			{
				fftw_complex* spec = nullptr;
				fftw_complex* spatial = nullptr;
				std::unique_ptr<GOAT::maths::fourier::fft2D> fft;

				int nx = 0;
				int ny = 0;

				~Workspace()
				{
					cleanup();
				}

				void cleanup()
				{
					if (spec)
					{
						fftw_free(spec);
						spec = nullptr;
					}

					fft.reset();
					nx = 0;
					ny = 0;
				}

				void prepare(int newNx, int newNy)
				{
					if (newNx == nx && newNy == ny && spec && fft)
						return;

					cleanup();

					nx = newNx;
					ny = newNy;

					std::size_t N =
						static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);

					spec = fftw_alloc_complex(N);
					spatial = fftw_alloc_complex(N);
					fft = std::make_unique<GOAT::maths::fourier::fft2D>(nx, ny);
				}
			};

			std::array<Workspace, 3> ws;

			
			AngularSpectrum& self;

			Impl(AngularSpectrum& self)
				: self(self)
			{
			}

			void cleanup()
			{
				for (auto& w : ws)
					w.cleanup(); 
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

			void applyTransferFunction(fftw_complex* spec, double dz, maths::Vector<double> shift)
			{
				if (spec == nullptr)
					throw std::invalid_argument("AngularSpectrum::applyTransferFunction: spec is null");

				const std::size_t nx = self.n1;   // oder getN1()
				const std::size_t ny = self.n2;   // oder getN2()

				const double k0 = 2.0 * M_PI / self.wvl;
				const double xShift = shift * self.gete1();
				const double yShift = shift * self.gete2();

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
						const double phase = kz * dz - xShift * kx - yShift * ky;

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

				void prepare(int n1, int n2)
				{
					for (auto& w : ws)
						w.prepare(n1, n2);
				}

			/*void propagateComponent(DetectorPlane* det, maths::fourier::fieldComponent component, double dz)
			{
				const std::size_t nx = static_cast<int>(self.n1);
				const std::size_t ny = static_cast<int>(self.n2);

				maths::fourier::vectorField2D tmp;
				int c = static_cast<std::size_t>(component);
				std::vector<std::vector<GOAT::maths::Vector<std::complex<double>>>> D;
				det->holographicField(D);
				
				// 1. Vorwärts-FFT: Quelldetektor -> Frequenzraum
				ws[c].fft->forward(D, ws[c].spec, component);

				// 2. Transferfunktion anwenden
				applyTransferFunction(ws[c].spec, dz);

				// 3. Rücktransformation
				ws[c].fft->inverse(ws[c].spec, tmp, component);

				// 4. Ergebnis auf Ziel addieren
			    addField(tmp, component);
				// addRoi(tmp, det, component);
			}*/

				void copySlmComponentToFftwInput(
					DetectorPlane* src,
					AngularSpectrum* dst,
					const maths::fourier::vectorField2D& slmField,
					maths::fourier::fieldComponent component,
					fftw_complex* in)
				{
					const int n1SLM = src->N1();
					const int n2SLM = src->N2();

					const int n1AS = dst->N1();
					const int n2AS = dst->N2();

					const double d1 = src->D1() / static_cast<double>(n1SLM);
					const double d2 = src->D2() / static_cast<double>(n2SLM);

					auto e1 = dst->gete1();
					auto e2 = dst->gete2();

					auto deltaP = src->position() - dst->position();

					const double slmMin1 = deltaP * e1 - 0.5 * src->D1();
					const double slmMin2 = deltaP * e2 - 0.5 * src->D2();

					const double asMin1 = -0.5 * dst->D1();
					const double asMin2 = -0.5 * dst->D2();

					const int offset1 =
						static_cast<int>(std::round((slmMin1 - asMin1) / d1));

					const int offset2 =
						static_cast<int>(std::round((slmMin2 - asMin2) / d2));

					/*if (offset1 < 0 || offset2 < 0 ||
						offset1 + n1SLM > n1AS ||
						offset2 + n2SLM > n2AS)
					{
						throw std::runtime_error("SLM field does not fit into AS field.");
					}*/

					const int c = static_cast<int>(component);

					for (int j = 0; j < n2SLM; ++j)
					{
						for (int i = 0; i < n1SLM; ++i)
						{
							const int iAS = offset1 + i;
							const int jAS = offset2 + j;

							const int idxAS = jAS * n1AS + iAS;

							const std::complex<double> E = slmField[i][j][c];

							in[idxAS][0] = E.real();
							in[idxAS][1] = E.imag();
						}
					}
				}

				void zeroFftwBuffer(fftw_complex* buffer, int n)
				{
					for (int k = 0; k < n; ++k)
					{
						buffer[k][0] = 0.0;
						buffer[k][1] = 0.0;
					}
				}

				void copyFftwToAsField(fftw_complex* in,
					maths::fourier::fieldComponent component)
				{
					const int n1 = static_cast<int>(self.N1());
					const int n2 = static_cast<int>(self.N2());

					const std::size_t c = static_cast<std::size_t>(component);

					const double scale =
						1.0 / static_cast<double>(n1 * n2);

					for (int j = 0; j < n2; ++j)
					{
						for (int i = 0; i < n1; ++i)
						{
							const int idx = j * n1 + i;

							// self.field(i, j)[c] =
							self.D[i][j][c]=
								std::complex<double>(
									in[idx][0] * scale,
									in[idx][1] * scale
								);
						}
					}
				}

			void propagateComponent(DetectorPlane* det, maths::fourier::fieldComponent component, double dz, maths::Vector<double> shift)
			{
				int c = static_cast<std::size_t>(component);
				// 1. SLM-Feld holen
				maths::fourier::vectorField2D slmField;
				det->holographicField(slmField);

				// 2. fftwInput auf AS-Größe nullen
				zeroFftwBuffer(ws[c].spec, self.n1 * self.n2);

				// 3. SLM-Komponente mit Offset in ws[c].in kopieren
				copySlmComponentToFftwInput(det, &self, slmField, component, ws[c].spatial);

				// 4. FFT
				ws[c].fft->forward(ws[c].spatial, ws[c].spec);

				// 5. Transferfunktion
				applyTransferFunction(ws[c].spec, dz, shift);

				// 6. IFFT
				ws[c].fft->inverse(ws[c].spec, ws[c].spatial);

				// 7. Ergebnis in AS-Feld schreiben
				copyFftwToAsField(ws[c].spatial, component);
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
			if (clear)
				clean();

			if (det == nullptr)
				throw std::invalid_argument("AngularSpectrumPropagator::calcOne: det is null");

		/*	// 1. Dimensionen prüfen
			if (det->N1() != this->N1() || det->N2() != this->N2())
				throw std::runtime_error("AngularSpectrumPropagator::calcOne: source and target grid size differ");
			*/	
			// 2. Abstand bestimmen
			maths::Vector<double> dP = P - det->position();
			maths::Vector<double> n = this->norm();   // oder det->getNormal(), wenn parallel

			if ((dP * n) < 0.0) n = -n;

			double dz = dP * n;

			// 3. Ebenenparallelität prüfen
			maths::Vector<double> nSrc = det->norm();
			maths::Vector<double> nDst = n;
			maths::Vector<double> shift = dP - dz * n;
			if ((abs(nSrc * nDst) - 1.0) > 1e-8)
				throw std::runtime_error("AngularSpectrumPropagator::calcOne: source and target plane are not parallel");

			// 4. Für jede Komponente propagieren
			maths::fourier::fieldComponent comp;
			impl->prepare(n1, n2);

		//	impl->prepare(det->N1(), det->N2());
			std::cout << "OpenMP aus" << std::endl;
// #pragma omp parallel for private(comp)
			for (int c = 0; c < 3; ++c)
			{
				comp = static_cast<maths::fourier::fieldComponent>(c);
				impl->propagateComponent(det, comp, dz, shift);
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
