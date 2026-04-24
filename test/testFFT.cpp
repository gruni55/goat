#include "angularSpectrum.h"

int main()
{
	GOAT::maths::Vector<double> P(0, 0, 0);
	GOAT::maths::Vector<double> e1(800, 0, 0);
	GOAT::maths::Vector<double> e2(0, 800, 0);
	int n1 = 1024;
	int n2 = 1024;
	double wvl = 1;
	GOAT::raytracing::DetectorPlane det(P, e1, e2, n1, n2);
	int m = 5;  // kleine Mode

	double radius = 20.0; // in deinen Längeneinheiten
	GOAT::maths::Vector<double> center = P + 0.5 * e1 + 0.5 * e2;

	for (size_t x = 0; x < det.N1(); ++x)
	{
		for (size_t y = 0; y < det.N2(); ++y)
		{
			double px = (static_cast<double>(x) + 0.5) / det.N1();
			double py = (static_cast<double>(y) + 0.5) / det.N2();

			GOAT::maths::Vector<double> r = P + px * e1 + py * e2;
			double rho = abs(r - center);

			if (rho <= radius)
				det.D[x][y] = GOAT::maths::Vector<std::complex<double>>({ 1.0, 0.0 }, 0.0, 0.0);
			else
				det.D[x][y] = GOAT::maths::Vector<std::complex<double>>(0.0, 0.0, 0.0);
		}
	}

	GOAT::raytracing::AngularSpectrum propagator(wvl, P + 1000 * GOAT::maths::ez, e1, e2, n1, n2);
	propagator.addDetector(&det);
	propagator.calc();
	propagator.save("C:\\tmp\\testFFT.dat");
	return 0;
}