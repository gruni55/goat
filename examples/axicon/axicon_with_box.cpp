#include "refractive_index_functions.h"
#include "pulsecalculation.h"

int main(int argc, char** argv)
{
	double r0 = 1E+5;  // radius of the calculation space
	GOAT::raytracing::Scene S;

	S.setr0(r0);

	// --- Objects ---
	GOAT::maths::Vector<double> axiconPos(0, 0, 0);
	GOAT::raytracing::surface axicon(axiconPos, 1.5);
	axicon.importBinSTL("axicon_100_10.stl");
	axicon.setActive(false);

	GOAT::maths::Vector<double> boxPos(115, 0, 0);
	GOAT::maths::Vector<double> boxDim(300, 100, 100);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.5);
	box.setActive(true);

	S.addObject(&axicon);
	S.addObject(&box);

	// --- Light source ---
	GOAT::maths::Vector<double> lsPos(-10, 0, 0);
	int numRays = 3000000;
	GOAT::raytracing::LightSrcPlane_mc ls(lsPos, numRays, 1.0, 200.0);
	ls.setk(GOAT::maths::ex);
	ls.setPol(GOAT::maths::Vector<std::complex<double> >(0.0,1.0,0.0));

	S.addLightSource(&ls);

	// --- Refractive index functions ---
	std::vector < std::function<std::complex<double>(double) > > nList;
	nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_Vacuum);
	nList.push_back(GOAT::raytracing::n_Vacuum);

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(50);
	pc.setRefractiveIndexFunctions(nList);
	pc.setSpatialResolution(1.0);
	pc.setReferenceTime(0.0);

	pc.field(1000.0);
	GOAT::raytracing::saveFullE(pc.trafo.SAres, "fieldtotal.dat", 1);
}
