#include "refractive_index_functions.h"
#include "pulsecalculation.h"

int main(int argc, char **argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(100000);

	// Add object
	GOAT::raytracing::surface boxObj;
	// boxObj.createsurface("H:\\Dokumente\\srf-files\\box_30_10_10.srf");	
	boxObj.importBinSTL("H:\\data\\axicon_r6.35mm_10grad.stl");
	boxObj.setActive(false);	
	 S.addObject(&boxObj);

	// Add calculation space
	GOAT::maths::Vector<double> calcBoxPos(25000, 0, 0);
	GOAT::maths::Vector<double> calcBoxDim(50000, 10000, 10000);
	GOAT::raytracing::Box calcBox(calcBoxPos, calcBoxDim, 1.5);
	calcBox.setActive(true);
	S.addObject(&calcBox);

	// Add light source
	GOAT::maths::Vector<double> LSPos(-50, 0, 0);
	int numRays = 10000;
	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, numRays, 1.0, 6350);
	LS.setk(GOAT::maths::Vector<double>(1.0, 0.0, 0.0));

	S.addLightSource(&LS);

	// --- Refractive index functions ---
	std::vector < std::function<std::complex<double>(double) > > nList;
	nList.push_back(GOAT::raytracing::n_BK7);    // for the box 
	nList.push_back(GOAT::raytracing::n_Vacuum); // for the calculation space
    nList.push_back(GOAT::raytracing::n_Vacuum); // for the surrounding medium

	// Now set the pulse calculation parameters...

	GOAT::raytracing::pulseCalculation pc(S);
	double pulseWidth = 50.0;
	double refTime = 0.0;
	double spatialRes = 100;
	
	pc.setSpatialResolution(spatialRes);
	pc.setPulseWidth(pulseWidth);
	pc.setReferenceTime(refTime);
	pc.setRefractiveIndexFunctions(nList);

	double time = 20;
	pc.field(time);
	GOAT::raytracing::saveFullE(pc.trafo.SAres,"/home/weigel/data/ellipsoid.dat",0);
  	GOAT::raytracing::saveFullE(pc.trafo.SAres,"/home/weigel/data/ellipsoid2.dat",1);

}
