#include "kirchhoff.h"
#include "raytrace.h"

void generatePinhole(double r,GOAT::raytracing::DetectorPlane *det)
{
	
	GOAT::maths::Vector<std::complex<double> > E(0, 1, 0);
	int N = det->N1();
	double d = det->D1() / (double)N;
	double D = det->D1();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double x = -D / 2. + d / 2. + i * d;
			double y = -D / 2. + d / 2. + j * d;
			if (sqrt(x * x + y * y) < r)
			{
				(*det)(i, j)=E;
			}
		}
	}
	
}

int main(int argv, char** argc)
{
	GOAT::maths::Vector<double> P(0, 0, 0);
	GOAT::maths::Vector<double> n(0, 0, 1);
	double D = 100;
	int N = 250;
	GOAT::raytracing::DetectorPlane* det = new GOAT::raytracing::DetectorPlane(P, n, D, 100);
	generatePinhole(10.0, det);
	GOAT::raytracing::Kirchhoff K(1.0, P + 1000.0 * GOAT::maths::ez, n, 10.0 * D, N);
	K.addDetector(det);
	K.setNumberOfThreads(20);
	K.calc();
	K.save("c:\\tmp\\kirchhofftest_100.dat");
	return 0;
}
