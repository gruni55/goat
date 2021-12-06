#include "raytrace.h"

int main(int argc, char** argv)
{
	LightSrcPlane LS(-200 * ex, 300, 1.0,50.0,I*ey);
	LS.setk(ex);
	double l = 5;
	Vector<double> d(l, 300, 300);
	Box Layer(dzero, d, 1.5);
	double alpha = 30.0 / 180.0 * M_PI;
	Layer.setBeta(alpha);

	Scene S;
	
	
	S.addLightSource(&LS);
	S.setRaytype(LIGHTSRC_RAYTYPE_IRAY);
    S.addObject(&Layer);
	S.setr0(1000.0);
	S.setnS(1.0);
	
	int N = 1;
	DetectorPlane D(200.0 * ex, -ex,500, N);
	DetectorPlane D0(-180.0 * ex, -ex, 400, N);
	S.addDetector(&D);
	S.addDetector(&D0);
	Raytrace_pure RP(S);
	RP.setNumReflex(2);
	
	Vector<std::complex<double> > E1, E2;
	std::ofstream os("C:\\Users\\Thomas\\tmp\\test.dat");
	 for (double d = 1.0; d <= 3.0; d += 0.01)
	// for (double wvl = 4.0; wvl <= 20.0; wvl += 0.01)
	{
		((Box *)RP.S.Obj[0])->setD(Vector<double> (d, 300, 300));
		// RP.S.LS[0]->wvl = wvl;
		RP.S.cleanAllDetectors();
		RP.trace();
		E1 = czero;
		E2 = czero;
		for (int ix = 0; ix < N; ix++)
			for (int iy = 0; iy < N; iy++)
			{
				E1 = E1 + RP.S.Det[0]->D[ix][iy];
				E2 = E2 + RP.S.Det[1]->D[ix][iy];
			}
		os << d << "\t" << abs2(E1) << "\t" << abs2(E2) << std::endl;
		std::cout << d << "\t" << abs2(E1)  << "\t" << abs2(E2) << std::endl;
	}
	os.close();


	return 0;
}
