#include "raytrace.h"

int main(int argc, char** argv)
{
	LightSrcPlane LS(-200 * ex, 10, 1.0,500);
	double alpha = 5.0 / 180.0 * M_PI;
	Vector<double> k(cos(alpha), sin(alpha), 0.0);
	std::cout << "k=" << k << std::endl;
	LS.setk(k);
	double l = 10;
	Vector<double> d(l, 100, 100);
	Box Layer=Box(dzero, d, 1.5);

	Scene S;
	S.setr0(1000.0);
	S.addLightSource(&LS);
    S.addObject(&Layer);
	
	int N = 100;
	DetectorPlane D(20.0 * ex, -ex, 100, N);
	S.addDetector(&D);

	Raytrace_pure RP(S);
	RP.setNumReflex(20);
	
	Vector<std::complex<double> > E1;
	std::ofstream os("C:\\Users\\Thomas\\tmp\\test.dat");
	for (double wvl = 1.0; wvl <= 5.0; wvl += 0.01)
	{
		RP.S.LS[0]->wvl = wvl;
		RP.S.cleanAllDetectors();
		RP.trace();
		E1 = czero;
		for (int ix = 0; ix < N; ix++)
			for (int iy = 0; iy < N; iy++)
			{
				E1 = E1 + RP.S.Det[0]->D[ix][iy];
			}
		os << wvl << "\t" << abs2(E1) << std::endl;
		std::cout << wvl << "\t" << abs2(E1)  << std::endl;
	}
	os.close();


	return 0;
}
