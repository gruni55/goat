#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>
#include <strstream>
/**
* This program provides a calculation of 
*/
int main (int argc, char **argv)
{
	std::string prismFName;
//	 prismFName = "/home/weigel/data/prism.srf";
    prismFName = "H:\\data\\prism_60.srf";
//	prismFName = "C:\\users\\thomas\\data\\prism.srf";
     std::string ergFName;
     ergFName = "C:\\users\\weigetz9\\data\\test.dat";
       	 // ergFName = "/home/weigel/data/test.dat";
//	ergFName = "C:\\users\\thomas\\data\\test2.dat";


 double sf=10000; // scaling factor (to convert µm into cm)
 double r0=1E+10;
 double wvl=0.5;
 double alpha = 30.0; // Brewster-Winkel
 // ------------ Light source ------------
 GOAT::maths::Vector<double> LSPos (0,-2*sf,0);
 int nRays=1;
 GOAT::raytracing::LightSrcPlane LS(LSPos, nRays,wvl,1E-10);
 LS.setPol(GOAT::maths::Vector<std::complex<double> > (0.0,1.0,0.0));
 LS.setk(GOAT::maths::ex);
 GOAT::maths::Vector<double> P1(30  * sf,  -4 * sf, -1 * sf);         // Position prism 1
 GOAT::maths::Vector<double> P2(60 * sf,  -23 * sf, -1 * sf);  // Position prism 2
 GOAT::maths::Vector<double> P3(90 * sf,  -23 * sf, -1 * sf); // Position prism 3
 GOAT::maths::Vector<double> P4(120  * sf, -4 * sf, -1 * sf);  // Position prism 4
 GOAT::maths::Vector<double> k1 = norm(P2 - P1);
 GOAT::maths::Vector<double> k2 = norm(P3 - P4);
 /*P2 = P2 - 5 * sf * k1;
 P3 = P3 - 5 * sf * k2;*/

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
 prism2.setGamma((180-alpha)/180.0*M_PI);
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
GOAT::maths::Vector<double> detPos(150 * sf,-2*sf,0);
GOAT::maths::Vector<double> detDim(1*sf,2,2);
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
/*
nList.push_back(GOAT::raytracing::n_Vacuum);
nList.push_back(GOAT::raytracing::n_Vacuum);
nList.push_back(GOAT::raytracing::n_Vacuum);
nList.push_back(GOAT::raytracing::n_Vacuum);
*/

nList.push_back (GOAT::raytracing::n_BK7);
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_BK7);
nList.push_back (GOAT::raytracing::n_BK7); 
nList.push_back (GOAT::raytracing::n_Vacuum);
nList.push_back (GOAT::raytracing::n_Vacuum);

// ----------- parameters for pulse calculation ------------
double pulseWidth = 100;
  double refTime = 1.7e+6;
  double spatialRes = 0.25;

  GOAT::raytracing::pulseCalculation pc(S);  
  pc.setPulseWidth (pulseWidth);  
  pc.setSpatialResolution (spatialRes);
  pc.setRefractiveIndexFunctions(nList);
 
  pc.setSpectralRanges(500);
  pc.setNumWavelengthsPerRange(10);
  pc.setCenterWavelength(wvl);
  pc.setNumReflex(0);  

// ------------ pulse calculation --------------

  double time = 136930;
  std::ofstream os;
  double fwhm;
  std::size_t fwhms;
  std::string fname;
  std::stringstream ss;
  os.open("C:\\users\\weigetz9\\data\\data_wvl3.dat");
  std::vector<double> d;
  std::vector<double> avgd;
  std::ofstream osall("C:\\users\\weigetz9\\data\\wvl3.dat");

  
  for (double wvl = 0.2; wvl <= 1.0; wvl += 0.005)
  {
    

      pc.setCenterWavelength(wvl);
      time = pc.findHitTime(4);
      pc.setReferenceTime(time);
      std::cout << "estimated time:" << time << std::endl;
      pc.field(time + 40000);
      
      for (int i = 0; i < 40000; i++)
          d.push_back(abs2(pc.trafo.SAres.G[4][i][4][4]));

      avgd = GOAT::maths::movingAvg(d, 50);
      for (int i = 0; i < 40000; i++)
      {
          osall << d[i] << "\t";
      }
      osall << std::endl << std::flush;
      std::size_t maxIndex;
      GOAT::maths::findmax(avgd, maxIndex);
      fwhms = GOAT::maths::FWHM(avgd, maxIndex);
      fwhm = fwhms * 0.25 / 0.3;
      os << wvl << "\t" << fwhm << "\t" << fwhms << std::endl;
      std::cout << "wvl: " << wvl << "µm\tFWHM: " << fwhm << std::endl;
      d.clear();
      avgd.clear();
  }
  os.close();
  osall.close();
  /*os.open("C:\\users\\weigetz9\\data\\data3.dat");
  for (alpha = 26; alpha < 27; alpha += 0.05)
  {
      prism1.setGamma((0.0 - alpha) / 180.0 * M_PI);
      prism2.setGamma((180 - alpha) / 180.0 * M_PI);
      prism3.setGamma((180 + alpha) / 180.0 * M_PI);
      prism4.setGamma((alpha) / 180.0 * M_PI);
      time = pc.findHitTime(4);
      std::cout << "estimated time:" << time << std::endl;
      pc.field(time + 40000);

      //  GOAT::raytracing::saveFullE(pc.trafo.SAres,ergFName,4); 

      std::vector<double> d;
      std::vector<double> avgd;
      for (int i = 0; i < 40000; i++)
          d.push_back(abs2(pc.trafo.SAres.G[4][i][4][4]));

      avgd = GOAT::maths::movingAvg(d, 50);
      std::size_t maxIndex;
      GOAT::maths::findmax(avgd, maxIndex);
      fwhms = GOAT::maths::FWHM(avgd, maxIndex);
      fwhm = fwhms * 0.25 / 0.3;
      os << alpha << "\t" << fwhm << "\t" << fwhms << std::endl;
      std::cout << "alpha: " << alpha << "°\tFWHM: " << fwhm << std::endl; 
  }*/
  os.close();
  /*
  alpha = 26.45;
  prism1.setGamma((0.0 - alpha) / 180.0 * M_PI);
  prism2.setGamma((180 - alpha) / 180.0 * M_PI);
  prism3.setGamma((180 + alpha) / 180.0 * M_PI);
  prism4.setGamma((alpha) / 180.0 * M_PI);*/

  /*
  time = pc.findHitTime(4);
  
  std::cout << "estimated time:" << time << std::endl;
  pc.field(time+20000);
  std::vector<double> d;
  std::vector<double> avgd;
  //GOAT::raytracing::saveFullE(pc.trafo.SAres, ergFName, 0);
  for (int i = 0; i < 40000; i++)
      d.push_back(abs2(pc.trafo.SAres.G[4][i][4][4]));
  avgd = GOAT::maths::movingAvg(d, 50);
  os.open("C:\\users\\weigetz9\\data\\test.dat");            

  for (int i = 0; i < 40000; i++)
  {
      os << d[i] << "\t" << avgd[i] << std::endl;      
  }
  os.close();
  */
}
