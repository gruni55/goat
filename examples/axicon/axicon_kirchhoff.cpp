#include "kirchhoff.h"
#include "lightsrc_mc.h"
#include "raytrace.h"
#include <chrono>

void makeThinAxicon(GOAT::raytracing::DetectorPlane& P)
{
	int n1 = P.N1();
	int n2 = P.N2();
	double L = P.D1();                       // Seitenlänge (µm), hier 250
	double dx = (n1 > 1) ? L / (n1 - 1) : 0;  // Schrittweite

	// Parameter
	double lambda = 1.0;                      // µm
	double k0 = 2.0 * M_PI / lambda;
	double n = 1.5;
	double alpha = 10.0 * M_PI / 180.0;       // in rad
	double kr = k0 * (n - 1.0) * std::tan(alpha);  // ≈ 0.554

	double cx = 0.5 * (n1 - 1);
	double cy = 0.5 * (n2 - 1);

	for (int i1 = 0; i1 < n1; ++i1)
	{
		for (int i2 = 0; i2 < n2; ++i2)
		{
			double x = (i1 - cx) * dx;
			double y = (i2 - cy) * dx;
			double r = std::sqrt(x * x + y * y);

			// konstante Amplitude, reine Axicon-Phase
			P.D[i1][i2] = GOAT::maths::ex * std::exp(I * kr * r);
		}
	}
}


void makeAxicon(GOAT::raytracing::DetectorPlane& P)
{
	double d = P.D1();
	double kr = 0.558;
	int n1 = P.N1();
	int n2 = P.N2();
	double dx = P.D1() / (double)n1;
	double x, y,r;
	for (int i1 = 0; i1<n1; i1++)
		for (int i2 = 0; i2 < n2; i2++)
		{
			x = i1 * dx - d / 2.0;
			y = i2 * dx - d / 2.0;
			r = sqrt(x * x + y * y);
			P.D[i1][i2]=GOAT::maths::ex * exp(I * kr * r);
		}
}


void makeSlit(GOAT::raytracing::DetectorPlane &P, double width)
{
	double d = P.D1();
	
	int n1 = P.N1();
	int dn = floor(width / d * n1);
	for (int i1 = n1/2 - dn / 2; i1 <= n1/2 + dn / 2; i1++)
		for (int i2 = 0; i2 < P.N2(); i2++)
			P.D[i1][i2] = GOAT::maths::Vector<std::complex<double> >(0.0, 1.0, 0.0);
}

void makeHole(GOAT::raytracing::DetectorPlane& P, double radius)
{
	double d = P.D1();
	int n1 = P.N1();
	int n2 = P.N2();
	int rn = floor(radius / d * (double)n1);
	int ri2;
	int rn2 = rn * rn;
	int ix, iy;
#pragma omp for collapse(6)
	for (int i1 = n1 / 2 - rn; i1 <= n1 / 2 + rn; i1++)
		for (int i2 = n2 / 2 - rn; i2 <= n2 / 2 + rn; i2++)
		{
			ix = i1 - n1 / 2;
			iy = i2 - n2 / 2;
			ri2 = ix * ix + iy * iy;
			if (ri2<=rn2) P.D[i1][i2] = GOAT::maths::Vector<std::complex<double> >(0.0, 1.0, 0.0);
		}
}


using namespace GOAT;
int main(int argc, char** argv)
{

	// ---- Light source ----
	maths::Vector<double> LSPos = -10000.0 * maths::ez;
	int numRays = 10000000;
	// numRays = 10;
	double wvl = 1;
	raytracing::LightSrcRing_mc LS(LSPos, numRays, wvl, 0, 125);
	LS.setk(maths::ez);
	LS.setPol(maths::Vector<std::complex<double>>(1.0, 0.0, 0.0));
	LS.setNumRays(numRays);

	// ---- Object ----
	maths::Vector<double> conePos = maths::dzero;
	double radius_um = 3000.0;
	double alpha = 10.0 / 180.0 * M_PI;
	double height_um = tan(alpha) * radius_um;
	std::complex<double> nCone = 1.5;
	raytracing::Cone c(conePos, radius_um, height_um, nCone);
	c.setConeAngle(80.0 / 180.0 * M_PI);
	// ---- Detector (faces source) ----
	const double eps_um = 1.0;
	maths::Vector<double> detPos(0, 0, height_um + eps_um);
	maths::Vector<double> detNorm(0, 0, -1);
	double detSize = 250.0;  int detGridsize = 250;
	raytracing::DetectorPlane det(detPos, detNorm, detSize, detGridsize);

	// ---- Scene ----
	raytracing::Scene S; S.setnS(1.0); S.setr0(1E+8);
	S.addLightSource(&LS); S.addObject(&c); S.addDetector(&det);
	std::cout << "raytracing...";
	raytracing::Raytrace_pure rp(S); 
	rp.setNumReflex(0);
	auto start = std::chrono::high_resolution_clock::now();  // Startzeitpunkt
	rp.trace();
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "done. (" << duration.count() << "ms)" << std::endl;
	// makeThinAxicon(det);
	// makeAxicon(det);
	// makeHole(det, 12.5);
	// makeSlit(det, 20);
	/*std::cout << "initialization...";
	makeHole(det, 50);
	std::cout << "done" << std::endl; */
	 // std::cout << "Tracing...\n"; rp.trace(); std::cout << "Trace done.\n";
	// S.Det[0]->D[50][50] = maths::Vector<std::complex<double>>(0, 1.0, 0);
	//----- now the Kirchhoff stuff ----

	maths::Vector<double> P;
	int n = 250; // number of cells/direction
	double l = 250; // edge length
	// maths::Vector<double> Pc(0, 0, 14000);
	 maths::Vector<double> Pc(0, 0, height_um+1000);

	raytracing::Kirchhoff kh(wvl, Pc, maths::ex * l, maths::ey * l, n, n);
	
	det.save("C:\\tmp\\detector.dat");
	kh.calc((raytracing::DetectorPlane *)&det);
	auto end2 = std::chrono::high_resolution_clock::now();
	auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end);
	std::cout << "Kirchhoff calculation took: " << duration2.count() << "ms" << std::endl;
	kh.save("C:\\tmp\\axicon.dat");
	
}