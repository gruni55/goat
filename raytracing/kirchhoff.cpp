#include "kirchhoff.h"
#include <omp.h>
namespace GOAT
{
	namespace raytracing
	{
		Kirchhoff::Kirchhoff(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2) : DetectorPlane(P, e1, e2, n1, n2)
		{
			k = 2.0 * M_PI / wvl;
			this->e1 = e1 / abs(e1);
			this->e2 = e2 / abs(e2);
			this->wvl = wvl;
			std::cout << "e1=" << this->e1 << "\t e2=" << this->e2 << std::endl;
		}

		void Kirchhoff::calc(DetectorPlane* det, bool clear)
		{
			if (clear) clean();
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
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(n1, n2, e1, e2, l1, l2, Pc, det, Dref)
			for (int i1=0; i1<n1; i1++)
				for (int i2 = 0; i2 < n2; i2++)
				{
					//P = Pc + (i1 / (double)n1 - 0.5) * d1 + (i2 / (double)n2 - 0.5) * d2; // Point at the Kirchhoff array
					maths::Vector<double> P =
						Pc + (i1 / (double)n1 - 0.5) * l1 * e1
						   + (i2 / (double)n2  - 0.5) * l2 * e2;
					Dref[i1][i2] += point(det, P);
				}
		}

		void Kirchhoff::calc(std::vector<DetectorPlane*> detList)
		{
			for (auto det : detList)
				calc(det);
		}

		maths::Vector<std::complex<double> > Kirchhoff::point(DetectorPlane* det, maths::Vector<double> P)
		{
			maths::Vector<double> R, Pc;
			maths::Vector<double> d1, d2;
			int n1, n2;
			n1 = det->N1();
			n2 = det->N2();
			d1 = det->gete1() * det->D1();
			d2 = det->gete2() * det->D2();
			double l1, l2;
			l1 = det->D1();
			l2 = det->D2();
			auto e1s = det->gete1();
			e1s /= abs(e1s);
			auto e2s = det->gete2();
			e2s /= abs(e2s);
			
		//	std::cout << "n1=" << n1 << "\tn2=" << n2 << "\te1 = " << e1 << "\t d1 = " << d1 << "\t e2 = " << e2 << "\t d2 = " << d2 << "\tD1 = " << D1() << "\t D2 = " << D2() << std::endl;

			Pc = det->position();
			maths::Vector<std::complex<double> > E;
			maths::Vector<double> rv;
			maths::Vector<std::complex<double> > s;
			double r;			
			GOAT::maths::Vector<std::complex<double> >cex(1, 0, 0);
			for (int i1 = 0; i1 < n1; i1++)
				for (int i2 = 0; i2 < n2; i2++)
				{
					R = Pc + (i1 / (double)n1 - 0.5) * l1 * e1s + (i2 / (double)n2 - 0.5) * l2 * e2s; // Point at raytracing area
					rv = P - R; 
					r = abs(rv);
					rv /= r;
					s = rv*(det->D[i1][i2] * rv);
				    // E += (det->D[i1][i2] - s) * exp(-I * k * r) / (I * wvl * r);
					 E += det->D[i1][i2]  * exp(-I * k * r) / (I * wvl * r);				
				}
			return E;
		}
	}
}