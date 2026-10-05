#include "angularSpectrum.h"
#include "xml.h"
#include "raytrace.h"

int main()
{
	GOAT::maths::Vector<double> P(0, 0, 0);
	GOAT::maths::Vector<double> e1(1600, 0, 0);
	GOAT::maths::Vector<double> e2(0, 1600, 0);
	int n1 = 512;
	int n2 = 512;
	double wvl = 0.5;
	GOAT::raytracing::DetectorPlane det(P, e1, e2, n1, n2);
	// det.setNorm(-GOAT::maths::ez);
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
	det.setID("source");

	GOAT::raytracing::AngularSpectrum propagator(wvl, P + 5000 * GOAT::maths::ez, e1, e2, n1, n2);

	propagator.addDetector(&det);
	propagator.calc();
	propagator.save("C:\\tmp\\testFFT.dat");
	GOAT::raytracing::Scene S;
	S.addDetector(&det);
	S.addDetector(&propagator);
	GOAT::XML::xmlWriter writer(S);
	writer.write("C:\\tmp\\testFFT.xml");

	return 0;
}