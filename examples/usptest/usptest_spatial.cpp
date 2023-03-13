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
	int numRays = 100;
	double LSSize = 100;
	GOAT::maths::Vector<double> LSDir(1.0, 0.0, 0.0);
	GOAT::raytracing::LightSrcPlane LS(LSPos, numRays, 1.0, LSSize);
	LS.setk(LSDir);
	

	GOAT::maths::Vector<double> ObjPos = GOAT::maths::dzero;
	GOAT::maths::Vector<double> ObjDim(2, 100, 100);
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
	pc.setPulseWidth(1E-15);
	pc.setRefractiveIndexFunctions(nList);
	pc.field(0.0);
	saveabsEbin(pc.trafo.SAres, "H:\\data\\field.dat");
	return 0;
}
/*int main(int argc, char** argv)
{
	// Define light source 
	GOAT::maths::Vector<double> LSPos = -300 * GOAT::maths::ex;
	GOAT::maths::Vector<double> LSDir = GOAT::maths::ex;
	GOAT::maths::Vector<std::complex<double> > LSPol(0.0, 1.0, 0.0);
	double LSd = 100;
	int LSNumRays = 1000000;
	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, LSNumRays, 1.0, LSd,LSPol);
	LS.setk(LSDir); 

	// Define Object
	GOAT::maths::Vector<double> ObjPos = GOAT::maths::dzero;
	GOAT::maths::Vector<double> ObjR = GOAT::maths::Vector<double>(30, 30, 30);
	GOAT::raytracing::Ellipsoid Obj(ObjPos, ObjR, 1.5);
	
	GOAT::raytracing::Scene S;
	S.addLightSource(&LS);
	S.addObject(&Obj);
	S.setr0(1000);

	GOAT::raytracing::Raytrace_Inel rt(S, 6000);
	rt.setExcitationFieldOnly();

	GOAT::raytracing::RRTParms rrtparms;
	rrtparms.coherency = GOAT::raytracing::INEL_RADIATION_INCOHERENT;
	rrtparms.P= 300 * GOAT::maths::ex;
	rrtparms.n = -GOAT::maths::ex;

	auto start = std::chrono::high_resolution_clock::now();
	rt.trace(rrtparms);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

	rt.exportExcitation("h:\\data\\usptest_bin");
  return 0;
}
*/