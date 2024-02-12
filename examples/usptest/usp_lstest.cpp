#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>

int main(int argc, char** argv)
{
	std::string ergFName = "C:\\users\\thomas\\data\\lstest.dat";
	GOAT::maths::Vector<double> LSPos(0, -10.0, 0);
	int nRays = 1;
	double wvl = 0.8;
	double r0 = 1E+10;
	GOAT::raytracing::LightSrcPlane LS(LSPos, nRays, wvl, 1E-10);
	LS.setPol(GOAT::maths::Vector<std::complex<double> >(1.0, 0.0, 0.0));
	LS.setk(GOAT::maths::ey);

	// ----------- Detector object --------------
	GOAT::maths::Vector<double> detPos(0, 5000, 0);
	GOAT::maths::Vector<double> detDim(2, 10000, 2);
	GOAT::raytracing::Box det(detPos, detDim, 1.0);
	det.setActive(true);

	// ------- refractive index list ------
	std::vector< std::function< std::complex< double >(double) > > nList;
	nList.push_back(GOAT::raytracing::n_Vacuum);
	nList.push_back(GOAT::raytracing::n_Vacuum);

	// ------------- Scene ---------------
	GOAT::raytracing::Scene S;
	S.setr0(r0);
	S.addLightSource(&LS);	
	S.addObject(&det);

	// ----------- parameters for pulse calculation ------------
	double pulseWidth = 1000;
	double refTime = 15000;
	double spatialRes = 0.25;

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setBandwidth(0.02);
	pc.setSpatialResolution(spatialRes);
	pc.setRefractiveIndexFunctions(nList);
	pc.setSpectralRanges(2000);
	pc.setNumWavelengthsPerRange(10);
	pc.setCenterWavelength(wvl);
	pc.setNumReflex(0);

	double time = 17000;
	pc.field(time);

	GOAT::raytracing::saveFullE(pc.trafo.SAres, ergFName, 0);

	// std::ofstream os("/home/weigel/data/lengths.dat");


	std::ofstream os;
	os.open("C:\\Users\\Thomas\\data\\lengths.dat");
	os << pc.rt.SA[0] << std::endl;
	os.close();
}