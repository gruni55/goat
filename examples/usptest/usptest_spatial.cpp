#include <chrono>
#include "pulsecalculation.h"
#include "refractive_index_functions.h"


std::complex<double> n_lin(double wvl)
{
	double n0 = 1.0;
	double m = 5E-4;
	double c = n0 - m;
	return m * wvl + c;
}


int main(int argc, char** argv)
{
	GOAT::maths::Vector<double> LSPos(-5E+5, 0, 0);
	double r0 = 1E6;
	int numRays = 100;
	double LSSize = 100;
	GOAT::maths::Vector<double> LSDir(1.0, 0.0, 0.0);
	// GOAT::raytracing::LightSrcGauss LS(LSPos, numRays, 1.0, 0.1, GOAT::maths::dzero);
	//LS.setNA(0.7);
   
	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, numRays, 1.0, LSSize);
	LS.setk(LSDir);
	LS.setD(5.0, 5.0);
	

	// GOAT::maths::Vector<double> ObjPos = GOAT::maths::dzero;
	GOAT::maths::Vector<double> ObjPos(0.0, 0.0, 0.0);
	GOAT::maths::Vector<double> ObjDim(100, 5, 5);
	GOAT::raytracing::Box Obj(ObjPos, ObjDim, 1.0);

	GOAT::raytracing::Scene S;
	S.setnS(1.5);
	S.setr0(r0);
	S.addLightSource(&LS);
	S.addObject(&Obj);

	std::vector<std::function<std::complex<double>(double) > > nList;
	/*nList.push_back(GOAT::raytracing::n_Air);*/
	/*nList.push_back(GOAT::raytracing::n_Air);
	nList.push_back(GOAT::raytracing::n_Air);*/
	/*nList.push_back(GOAT::raytracing::n_Vacuum);
	nList.push_back(GOAT::raytracing::n_Vacuum);*/

	nList.push_back(n_lin);
	nList.push_back(n_lin);

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(100);
	pc.setRefractiveIndexFunctions(nList);
	double tref = 1E+6;
	pc.setReferenceTime(tref);
	pc.setSpatialResolution(1.0);

	double t = (5E+5-250)/GOAT::raytracing::C_LIGHT_MU_FS*1.0;
	pc.field(t);
	GOAT::raytracing::saveEyPol (pc.trafo.SAres, "field.dat");
	GOAT::raytracing::saveabsE(pc.trafo.SAres, "fieldabs.dat");
	return 0;
}
