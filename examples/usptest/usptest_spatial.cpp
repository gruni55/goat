#include <chrono>
#include "pulsecalculation.h"
#include "refractive_index_functions.h"



int main(int argc, char** argv)
{
	GOAT::maths::Vector<double> LSPos(-5000, 0, 0);
	int numRays = 100000;
	double LSSize = 100;
	GOAT::maths::Vector<double> LSDir(1.0, 0.0, 0.0);
	// GOAT::raytracing::LightSrcGauss LS(LSPos, numRays, 1.0, 0.1, GOAT::maths::dzero);
	 GOAT::raytracing::LightSrcPlane_mc LS(LSPos, numRays, 1.0, LSSize);
	LS.setk(LSDir);
	LS.setD(100.0, 5.0);
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
	nList.push_back(GOAT::raytracing::n_Air);
	nList.push_back(GOAT::raytracing::n_BK7);

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(100E-15);
	pc.setRefractiveIndexFunctions(nList);
	double tref = 1E6;
	pc.setReferenceTime(tref);

	double t = (1E+6 + 50)/GOAT::raytracing::C_LIGHT_MU;
	pc.field(t);
	GOAT::raytracing::saveEyPol (pc.trafo.SAres, "field.dat");
	GOAT::raytracing::saveabsE(pc.trafo.SAres, "fieldabs.dat");
	return 0;
}
