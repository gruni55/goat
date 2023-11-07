#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>

int main (int argc, char **argv)
{
	std::string prismFName;
	// prismFName = "/home/weigel/data/prism.srf";
	prismFName = "H:\\data\\prism.srf";

 double sf=10000;
 double r0=10000*sf;
 // ------------ Light source ------------
 GOAT::maths::Vector<double> LSPos (0,-5*sf,0);
 int nRays=3;
 GOAT::raytracing::LightSrcPlane LS(LSPos, nRays,1,2);
 LS.setPol(GOAT::maths::Vector<std::complex<double> > (1.0,0.0,0.0));
 LS.setk(GOAT::maths::ey);

 // ------------ Prism 1 -------------
 GOAT::raytracing::surface prism1;
 prism1.setn(1.5);
 prism1.createsurface (prismFName);
 prism1.setGamma(90.0/180.0*M_PI);
 prism1.setPos(2.5*sf,0,-1*sf);
 prism1.setActive(false);

// ------------- Prism 2 -------------
 GOAT::raytracing::surface prism2;
 prism2.setn(1.5);
 prism2.createsurface (prismFName);
 prism2.setGamma(270.0/180.0*M_PI);
 prism2.setPos(1*sf,5*sf,-1*sf);
 prism2.setActive(false);

// ------------- Prism 3 -------------
 GOAT::raytracing::surface prism3;
 prism3.setn(1.5);
 prism3.createsurface (prismFName);
 prism3.setGamma(270.0/180.0*M_PI);
 prism3.setPos(1*sf,10*sf,-1*sf);
 prism3.setActive(false);

// ------------- Prism 4 -------------

 GOAT::raytracing::surface prism4;
 prism4.setn(1.5);
 prism4.createsurface (prismFName);
 prism4.setGamma(90.0/180.0*M_PI);
 prism4.setPos(2.5*sf,15*sf,-1*sf);
 prism4.setActive(false);

// ----------- Detector object --------------
GOAT::maths::Vector<double> detPos(0,25*sf,0);
GOAT::maths::Vector<double> detDim(3,12000,3);
GOAT::raytracing::Box det(detPos,detDim,1.0);
det.setActive(true); 

 
// ------------- Scene ---------------
GOAT::raytracing::Scene S;
S.addLightSource(&LS);
S.addObject(&prism1);
S.addObject(&prism2);
S.addObject(&prism3);
S.addObject(&prism4);
S.addObject(&det);
S.setr0(r0);

// ------- refractive index list ------
std::vector< std::function< std::complex< double >(double) > > nList;
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_Vacuum);
nList.push_back (GOAT::raytracing::n_Vacuum);

// ----------- parameters for pulse calculation ------------
  double pulseWidth = 100;
  double refTime = 120000;
  double spatialRes = 0.8;

  GOAT::raytracing::pulseCalculation pc(S);
  pc.setPulseWidth (pulseWidth);
  pc.setSpatialResolution (spatialRes);
  pc.setRefractiveIndexFunctions(nList);
  pc.setSpectralRanges (100);
  pc.setNumWavelengthsPerRange(100);

// ------------ pulse calculation --------------
double time=120000.0; 
pc.field (time);

GOAT::raytracing::saveFullE(pc.trafo.SAres,"C:\\Users\\Thomas\\data\\test1.dat",4); 
// std::ofstream os("/home/weigel/data/lengths.dat");

/*
std::ofstream os("H:\\data\\lengths.dat");
os << pc.rt.SA[0] << std::endl;
os.close();
*/
}