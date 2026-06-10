#include "kirchhoff.h"
#include <omp.h>
namespace GOAT
{
	namespace raytracing
	{
		Kirchhoff::Kirchhoff(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2) 
			: Propagator(wvl, P, e1, e2, n1, n2)
		{		
			type = DETECTOR_KIRCHHOFF;
		}

		Kirchhoff::Kirchhoff(double wvl, maths::Vector<double> P, maths::Vector<double> n, double d, int N) : Propagator(wvl, P, n, d, N)
		{
			type = DETECTOR_KIRCHHOFF;
		}


		void smoothField(std::vector<std::vector<maths::Vector<std::complex<double>>>>& E)
		{
			int n1 = E.size();
			if (n1 == 0) return;

			int n2 = E[0].size();

			auto tmp = E;

			for (int i = 1; i < n1 - 1; ++i)
			{
				for (int j = 1; j < n2 - 1; ++j)
				{
					tmp[i][j] =
						(E[i - 1][j - 1] + 2.0 * E[i][j - 1] + E[i + 1][j - 1]
							+ 2.0 * E[i - 1][j] + 4.0 * E[i][j] + 2.0 * E[i + 1][j]
							+ E[i - 1][j + 1] + 2.0 * E[i][j + 1] + E[i + 1][j + 1]) / 16.0;
				}
			}

			// Ränder unverändert übernehmen
			for (int i = 1; i < n1 - 1; ++i)
			{
				tmp[i][0] = E[i][0];
				tmp[i][n2 - 1] = E[i][n2 - 1];
			}

			for (int j = 0; j < n2; ++j)
			{
				tmp[0][j] = E[0][j];
				tmp[n1 - 1][j] = E[n1 - 1][j];
			}

			E.swap(tmp);
		}

		void Kirchhoff::calcOne(DetectorPlane* det, bool clear)
		{
			if (clear) clean();
			std::cout << "Calculating Kirchhoff with wavelength " << wvl << " for detector " << det->getID() << std::endl;
			std::cout << "det->D1() = " << det->D1()	<< "   det->N1() = " << det->N1() << std::endl;
			// Kirchhoff - Ebene
			maths::Vector<double> P,Pc;
			maths::Vector<double> d1, d2;
			auto& Dref = D;
			int n1, n2;
			n1 = N1();
			n2 = N2();
			auto e1 = this->e1 / abs(this->e1);
			auto e2 = this->e2 / abs(this->e2);
			
			double l1, l2;
			l1 = D1();
			l2 = D2();
			Pc = position();
		/*	const double invN1 = (n1 > 0) ? 1.0 / n1 : 0.0;
			const double invN2 = (n2 > 0) ? 1.0 / n2 : 0.0;
		*/

			// smoothField(det->D);
 #pragma omp parallel for collapse(2) schedule(static) default(none) shared(n1, n2, e1, e2, l1, l2, Pc, det, Dref) num_threads(noThreads)
			for (int i1=0; i1<n1; i1++)
				for (int i2 = 0; i2 < n2; i2++)
				{
					//P = Pc + (i1 / (double)n1 - 0.5) * d1 + (i2 / (double)n2 - 0.5) * d2; // Point at the Kirchhoff array
					maths::Vector<double> P =
						Pc + (i1 / (double)n1 - 0.5) * l1 * e1
						   + (i2 / (double)n2  - 0.5) * l2 * e2;
					Dref[i1][i2] += point(det, P, wvl);
				}
		}

		

		maths::Vector<std::complex<double>> point(DetectorPlane* det, maths::Vector<double> P, double wvl)
		{
			using GOAT::maths::Vector;

			const int n1 = det->N1();
			const int n2 = det->N2();

			const double l1 = det->D1();
			const double l2 = det->D2();

			const double dU = l1 / static_cast<double>(n1);
			const double dV = l2 / static_cast<double>(n2);
			const double dA = dU * dV;

			Vector<double> pc = det->position();

			Vector<double> e1s = det->gete1();
			e1s /= abs(e1s);

			Vector<double> e2s = det->gete2();
			e2s /= abs(e2s);

			Vector<double> n = det->norm();
			n /= abs(n);

			const double k = 2.0 * M_PI / wvl;

			std::complex<double> I(0.0, 1.0);
			Vector<std::complex<double>> E(0.0, 0.0, 0.0);

			// 2x2 Subsampling innerhalb jedes Pixels
			// relative Offsets innerhalb der Pixelzelle
			constexpr double subOffsets[2] = { -0.25, +0.25 };
			constexpr double subWeight = 0.25; // 4 Unterpunkte -> je 1/4

			for (int i1 = 0; i1 < n1; ++i1)
			{
				for (int i2 = 0; i2 < n2; ++i2)
				{
					const Vector<std::complex<double>>& field = det->D[i1][i2];

					// Pixelzentrum
					const double uCenter = (((static_cast<double>(i1) + 0.5) / static_cast<double>(n1)) - 0.5) * l1;
					const double vCenter = (((static_cast<double>(i2) + 0.5) / static_cast<double>(n2)) - 0.5) * l2;

					for (double a : subOffsets)
					{
						for (double b : subOffsets)
						{
							const double u = uCenter + a * dU;
							const double v = vCenter + b * dV;

							Vector<double> R = pc + u * e1s + v * e2s;

							Vector<double> rv = P - R;
							const double r = abs(rv);
							if (r < 1e-12)
								continue;

							rv /= r;

							double cosTheta = rv * n;
							if (cosTheta <= 0.0)
								continue;

							// Falls du den transversalen Anteil testen willst:
							// Vector<std::complex<double>> s = rv * (field * rv);
							// Vector<std::complex<double>> fieldTrans = field - s;

							std::complex<double> phase = std::exp(-I * k * r);

							E += field * (subWeight * dA * cosTheta) * phase / (I * wvl * r);
						}
					}
				}
			}

			return E;
		}
		
		Kirchhoff3D::Kirchhoff3D(Box* box, INDEX_TYPE nn)
		{		
			field3D = raytracing::SuperArray<maths::Vector<std::complex<double>>>(box->getr0(), nn, nn, nn);
			box->setActive(true);
			field3D.addInc(box);
			fieldInitialized = true;
			this->box = box;
		}

		void Kirchhoff3D::addDetector(DetectorPlane* det)
		{
			sources.push_back(det);
		}

		void Kirchhoff3D::addDetectorList(std::vector<DetectorPlane*> detList)
		{
			for (auto det : detList)
				sources.push_back(det);
		}

		void Kirchhoff3D::setR0(double r0)
		{
			if (r0 != box->getr0())
			{
				box->setr0(r0);
				std::cout << "[setR9]" << "\ttype=" << field3D.type << std::endl;

				// field3D = raytracing::SuperArray<maths::Vector<std::complex<double>>>(box->getr0(), field3D.n[0][0], field3D.n[0][1], field3D.n[0][2]);
				box->setActive(true);
				field3D.reinit(r0, field3D.nges[0], field3D.nges[1], field3D.nges[2]);
			}
		}

		void Kirchhoff3D::setNN(INDEX_TYPE nn)
		{
			std::cout << "[setNN]" << "\t nn=" << nn << "\ttype=" << field3D.type << std::endl;
			field3D.setNumberOfCellsPerDirection(nn);
			fieldInitialized = true;
		}

		void Kirchhoff3D::setSpatialResolution(double res)
		{
			INDEX_TYPE nn = (INDEX_TYPE)ceil(2.0 * box->getr0() / res);
			setNN(nn);
		}

		void Kirchhoff3D::calc(double wvl, int noThreads)
		{
			double k = 2.0 * M_PI / wvl;
			for (auto det : sources)
				calc(det, wvl,noThreads, false);
		}

		void Kirchhoff3D::calc(DetectorPlane* det, double wvl, int noThreads, bool clear)
		{
			if (!fieldInitialized)
			{
				field3D = raytracing::SuperArray<maths::Vector<std::complex<double>>>(box->getr0(), field3D.n[0][0], field3D.n[0][1], field3D.n[0][2]);
				fieldInitialized = true;
			}
			if (clear) field3D.fill(maths::czero);
			double d = box->getr0() * 2.0 / (double)(field3D.nges[0] - 1);
			maths::Vector<double> hd = box->d;
			maths::Vector<double> Pc = box->getPos();
			auto* field = &field3D;
			for (auto det : sources)
            #pragma omp parallel for collapse(3) schedule(static) default(none) shared(field,Pc,wvl,d,hd,det) num_threads(noThreads)
			for (INDEX_TYPE ix=0; ix < field->n[0][0]; ix++)
				for (INDEX_TYPE iy = 0; iy < field->n[0][1]; iy++)
				{
					for (INDEX_TYPE iz = 0; iz < field->n[0][2]; iz++)
					{
						maths::Vector<double> P = Pc + maths::Vector<double>((ix / (double)field->n[0][0] - 0.5) * hd[0], (iy / (double)field->n[0][1] - 0.5) * hd[1], (iz / (double)field->n[0][2] - 0.5) * hd[2]);
						(*field)(0, ix, iy, iz) += point(det, P, wvl);
					}
				}
		}

	} // namespace raytracing
} // namespace GOAT