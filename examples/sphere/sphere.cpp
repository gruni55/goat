#include "refractive_index_functions.h"
#include "pulsecalculation.h"

int main(int argc, char** argv)
{
	double r0 = 1E+5;

	double r = 5000.0;
	GOAT::maths::Vector<double> Pos(0, 0, 0);
	GOAT::maths::Vector<double> dim(r, r, r);
	GOAT::raytracing::Ellipsoid E(Pos, dim, 1.5, r0);
	E.setActive(true);

	GOAT::maths::Vector<double> LSPos(-1.2*r, 0, 0);
	int numRays = 10;

	GOAT::raytracing::LightSrcPlane LS(LSPos, numRays, 1.0, r);
	LS.setk(GOAT::maths::ex);

	GOAT::raytracing::Scene S;
	S.addLightSource(&LS);
	S.addObject(&E);
	S.setr0(r0);

	GOAT::raytracing::Raytrace_Path rp(S);
	rp.trace("test.dat"); 

	std::vector<std::function<std::complex<double>(double) > > nList;
	nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_Air);

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(100);
	pc.setRefractiveIndexFunctions(nList);
	pc.setSpatialResolution(1000.0);
	pc.setReferenceTime(0.0);
	
	double t = 3400.0;
	pc.field(t);
	GOAT::raytracing::saveabsE(pc.trafo.SAres, "fieldabs.dat", 0);
	return 0;
}
