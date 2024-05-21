#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>
#include <strstream>
#include <string>
/**
* This program provides a calculation of 
*/
int main (int argc, char **argv)
{
	std::string prismFName;
	 prismFName = "/home/weigel/data/prism.srf";
//    prismFName = "H:\\data\\prism_60.srf";
	prismFName = "C:\\users\\thomas\\data\\prism.srf";
     std::string ergFName;
//     ergFName = "C:\\users\\weigetz9\\data\\test.dat";
       	 // ergFName = "/home/weigel/data/test.dat";
	ergFName = "C:\\users\\thomas\\data\\test2.dat";


 double sf=10000; // scaling factor (to convert �m into cm)
 double r0=1E+10;
 double wvl=0.5;
 double alpha = 30.0; // Brewster-Winkel
 // ------------ Light source ------------
 GOAT::maths::Vector<double> LSPos (0,0,0);
 int nRays=5000;
 double r=0.2*sf; // radius of the light source in cm
 // GOAT::raytracing::LightSrcPlane_mc LS(LSPos, nRays,wvl,1E-10);
GOAT::raytracing::LightSrcGauss_mc LS(LSPos, nRays, wvl, r, GOAT::maths::ex * 100 * sf,2*r);
 LS.setPol(GOAT::maths::Vector<std::complex<double> > (0.0,1.0,0.0));
 LS.setk(GOAT::maths::ex);
 
 GOAT::maths::Vector<double> P1(30  * sf,  -3 * sf, -1 * sf);         // Position prism 1
 GOAT::maths::Vector<double> P2(60 * sf,  -12 * sf, -1 * sf);  // Position prism 2
 GOAT::maths::Vector<double> P3(90 * sf,  -12 * sf, -1 * sf); // Position prism 3
 GOAT::maths::Vector<double> P4(120  * sf, -3 * sf, -1 * sf);  // Position prism 4
 GOAT::maths::Vector<double> k1 = norm(P2 - P1);
 GOAT::maths::Vector<double> k2 = norm(P3 - P4);

 // ------------ Prism 1 -------------
 GOAT::raytracing::surface prism1;
 prism1.setn(1.5);
 prism1.createsurface (prismFName);
 prism1.setGamma((0.0-alpha)/180.0*M_PI);
 prism1.setPos(P1);
 prism1.setActive(false);

// ------------- Prism 2 -------------
 GOAT::raytracing::surface prism2;
 prism2.setn(1.5);
 prism2.createsurface (prismFName);
 prism2.setGamma((180-alpha+0.5)/180.0*M_PI);
 prism2.setPos(P2);
 prism2.setActive(false);

// ------------- Prism 3 -------------
 GOAT::raytracing::surface prism3;
 prism3.setn(1.5);
 prism3.createsurface (prismFName);
 prism3.setGamma((180+alpha)/180.0*M_PI);
 prism3.setPos(P3);
 prism3.setActive(false);

// ------------- Prism 4 -------------

 GOAT::raytracing::surface prism4;
 prism4.setn(1.5);
 prism4.createsurface (prismFName);
 prism4.setGamma((alpha)/180.0*M_PI);
 prism4.setPos(P4);
 prism4.setActive(false);

// ----------- Detector object --------------
GOAT::maths::Vector<double> detPos(150 * sf,-0.225*sf,0);
GOAT::maths::Vector<double> detDim(10,0.25*sf,0.25*sf);
GOAT::raytracing::Box det(detPos,detDim,1.0);
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

nList.push_back (GOAT::raytracing::n_BK7);
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_BK7);
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_Vacuum);
nList.push_back (GOAT::raytracing::n_Vacuum);

// ----------- parameters for pulse calculation ------------
  double pulseWidth = 50;
  double refTime = 1.7e+6;
  double spatialRes = 5;

  GOAT::raytracing::pulseCalculation pc(S);  
  pc.setPulseWidth (pulseWidth);  
  pc.setSpatialResolution (spatialRes);
  pc.setRefractiveIndexFunctions(nList);
 
  pc.setSpectralRanges(100);
  pc.setNumWavelengthsPerRange(1);
  pc.setCenterWavelength(wvl);
  pc.setNumReflex(0);  

// ------------ pulse calculation --------------

  double time;
  std::ofstream os;
  double fwhm;
  std::size_t fwhms;
  std::vector<double> d;
  std::vector<double> avgd;




// ------------- move prisms -----------------
int counter=0;
std::string fname="/home/weigel/data/misalign_hr.dat";

   LS.setNumRays(1);
   LS.setD(1);
   time = pc.findHitTime(4);

  std::cout << "% estimated time:" << time << std::endl;
//   pc.setReferenceTime(time);
  LS.setNumRays(nRays);
  LS.setD(2*r);

  // time=5.33925e+06;
  double D;
  std::ofstream oserr;
  oserr.open("c:\\users\\thomas\\data\\testerr.dat");
  do
  {
	  D = pc.field(time, GOAT::raytracing::PULSECALCULATION_NOT_CLEAR_RESULT);
	  oserr <<  D << std::endl;
	  GOAT::raytracing::saveFullE(pc.trafo.SAres, fname, 4);
  } while (D > 0.1);
 std::cout << "% current filename: " << fname << std::endl;
}
