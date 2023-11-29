#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>

int main(int argc, char** argv)
{
	std::string prismFName;
	// prismFName = "/home/weigel/data/prism.srf";
	prismFName = "H:\\data\\prism.srf";
	std::string ergFName;
	// ergFName = "/home/weigel/data/test2.dat";
	ergFName = "C:\\users\\thomas\\data\\test.dat";


	double sf = 10000;
	double r0 = 1E+10;
	double wvl = 0.5;
	// ------------ Light source ------------
	GOAT::maths::Vector<double> LSPos(0, -50 * sf, 0);
	int nRays = 1;
	GOAT::raytracing::LightSrcPlane LS(LSPos, nRays, wvl, 1E-5);
	LS.setPol(GOAT::maths::Vector<std::complex<double> >(1.0, 0.0, 0.0));
	LS.setk(GOAT::maths::ey);

	double alpha = 44.3;

	// ------------ Prism 1 -------------
	GOAT::raytracing::surface prism1;
	prism1.setn(1.5);
	prism1.createsurface(prismFName);
	prism1.setGamma((90.0 - alpha) / 180.0 * M_PI);
	prism1.setPos(3 * sf, -1, -1 * sf);
	prism1.setActive(false);
	//  std::cout << "%->"<<  prism1.Tree.BBox << std::endl;

	// ------------- Prism 2 -------------
	GOAT::raytracing::surface prism2;
	prism2.setn(1.5);
	prism2.createsurface(prismFName);
	prism2.setGamma((270.0 - alpha) / 180.0 * M_PI);
	prism2.setPos(37.0 * sf, 65 * sf, -1 * sf);
	prism2.setActive(false);

	// ------------- Prism 3 -------------
	GOAT::raytracing::surface prism3;
	prism3.setn(1.5);
	prism3.createsurface(prismFName);
	prism3.setGamma((270.0 + alpha) / 180.0 * M_PI);
	prism3.setPos(37.0 * sf, 100 * sf, -1 * sf);
	prism3.setActive(false);

	// ------------- Prism 4 -------------

	GOAT::raytracing::surface prism4;
	prism4.setn(1.5);
	prism4.createsurface(prismFName);
	prism4.setGamma((90.0 + alpha) / 180.0 * M_PI);
	prism4.setPos(3 * sf, 165 * sf, -1 * sf);
	prism4.setActive(false);

	// ----------- Detector object --------------
	GOAT::maths::Vector<double> detPos(0, 2E+6, 0);
	GOAT::maths::Vector<double> detDim(2, 10000, 2);
	GOAT::raytracing::Box det(detPos, detDim, 1.0);
	det.setActive(true);


	// ------------- Scene ---------------
	GOAT::raytracing::Scene S;
	S.setr0(r0);
	S.addLightSource(&LS);
	S.addObject(&prism1);
	S.addObject(&prism2);
	S.addObject(&prism3);
	S.addObject(&prism4);
	S.addObject(&det);


	// ------- refractive index list ------
	std::vector< std::function< std::complex< double >(double) > > nList;
	nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_Vacuum);
	nList.push_back(GOAT::raytracing::n_Vacuum);

	// ----------- parameters for pulse calculation ------------
	double pulseWidth = 1000;
	double refTime = 9e+6;
	double spatialRes = 0.25;

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(pulseWidth);
	pc.setSpatialResolution(spatialRes);
	pc.setRefractiveIndexFunctions(nList);
	pc.setSpectralRanges(1);
	pc.setNumWavelengthsPerRange(1);
	/*  pc.setSpectralRanges (400);
	  pc.setNumWavelengthsPerRange(100);*/
	pc.setCenterWavelength(wvl);
	pc.setNumReflex(0);

	// ------------ pulse calculation --------------
	// double time=1.2E+6-1.5E+4; 
	  // double time = 136930;
	 // double time = 1718233.107+120000;
	double time = 9098230 + 30000;
	pc.field(time);

	GOAT::raytracing::saveFullE(pc.trafo.SAres, ergFName, 4);

	// std::ofstream os("/home/weigel/data/lengths.dat");


	std::ofstream os;
	os.open("C:\\Users\\weigetz9\\data\\lengths.dat");
	// os.open("/home/weigel/data/lengths2.dat");
	os << pc.rt.SA[0] << std::endl;
	os.close();

}
