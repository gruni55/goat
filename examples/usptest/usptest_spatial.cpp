#include <chrono>
#include "pulsecalculation.h"

std::complex<double> none(double wvl)
{
	return 1.0;
}

/**
* Wavelength dependent refractive index for glass (BK7) taken from refractiveindex.info
*/

std::complex<double> nGlass_BK7(double wvl)
{
	double wvl2 = wvl * wvl;
	std::complex<double> n = sqrt(1.03961212*wvl2/(wvl2- 0.00600069867)+ 0.231792344*wvl2/(wvl2- 0.0200179144)+ 1.01046945*wvl2/(wvl2- 103.560653)+1.0);
	return n;
}

std::complex<double> nGlass(double wvl)
{
	return 1.5075;
}



int main(int argc, char** argv)
{
	GOAT::maths::Vector<double> LSPos(-5000, 0, 0);
	int numRays = 1000000;
	double LSSize = 100;
	GOAT::maths::Vector<double> LSDir(1.0, 0.0, 0.0);
	// GOAT::raytracing::LightSrcGauss LS(LSPos, numRays, 1.0, 0.1, GOAT::maths::dzero);
	 GOAT::raytracing::LightSrcPlane_mc LS(LSPos, numRays, 1.0, LSSize);
	LS.setk(LSDir);
	//LS.setNA(0.7);

	GOAT::maths::Vector<double> ObjPos = GOAT::maths::dzero;
	GOAT::maths::Vector<double> ObjDim(100, 100, 5);
	GOAT::raytracing::Box Obj(ObjPos, ObjDim, 1.0);

	GOAT::raytracing::Scene S;
	S.setnS(1.0);
	S.setr0(10000.0);
	S.addLightSource(&LS);
	S.addObject(&Obj);

	std::vector<std::function<std::complex<double>(double) > > nList;
	nList.push_back(none);
	nList.push_back(none);

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(100E-15);
	pc.setRefractiveIndexFunctions(nList);
	double tref = 1E6;
	pc.setReferenceTime(tref);

	double t = (1E+6 + 50)/GOAT::raytracing::C_LIGHT_MU;
	pc.field(t);
	GOAT::raytracing::saveEyPol (pc.trafo.SAres, "H:\\data\\field.dat");
	GOAT::raytracing::saveabsE(pc.trafo.SAres, "H:\\data\\fieldabs.dat");
	return 0;
}
