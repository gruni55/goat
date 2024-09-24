#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>
#include <strstream>
#include <string>


int main (int argc, char **argv)
{

 double r0=1E+10;
 GOAT::raytracing::Scene S;
 S.setr0(r0);

 GOAT::maths::Vector<double> LSPos;
 int nrays=1;

 GOAT::raytracing::LightSrcPlane LS(LSPos,nrays,0.5,1E-10);
 LS.setk(GOAT::maths::ex);

 double cm=1E+4;
 double mm=1E+3;
 GOAT::maths::Vector<double> boxPos(3000.0*cm,0,0);
 GOAT::maths::Vector<double> boxDim(1.0*cm,0.25,0.25);
 GOAT::raytracing::Box box(boxPos,boxDim,1.0);
 box.setActive(true);

 GOAT::maths::Vector<double> objPos(1000.0*cm,0,0);
 GOAT::maths::Vector<double> objDim(100.0*cm,100.0,100.0);
 GOAT::raytracing::Box obj(objPos,objDim,1.0);
 obj.setActive(false);

 S.addLightSource(&LS);
 S.addObject(&obj);
 S.addObject(&box);


 std::vector< std::function< std::complex< double >(double) > > nList;
 nList.push_back (GOAT::raytracing::n_BK7);
 nList.push_back (GOAT::raytracing::n_Vacuum);
 nList.push_back (GOAT::raytracing::n_Vacuum);

  double pulseWidth = 50;
  double refTime = 1.7e+6;
  double spatialRes = 1;
  double wvl = 0.5;

  GOAT::raytracing::pulseCalculation pc(S);
  pc.setPulseWidth (pulseWidth);
  pc.setSpatialResolution (spatialRes);
  pc.setRefractiveIndexFunctions(nList);

  pc.setSpectralRanges(1000);
  pc.setNumWavelengthsPerRange(1);
  pc.setCenterWavelength(wvl);
  pc.setNumReflex(0);

  double time = pc.findHitTime(1);
//  std::cout << "estimated time: " << time << "fs" << std::endl;
  pc.field(1.0002E8);
  GOAT::raytracing::saveFullE(pc.trafo.SAres, "/home/weigel/data/testusp.dat", 1);
 return 0;
}

