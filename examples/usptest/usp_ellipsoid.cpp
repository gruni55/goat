#include "pulsecalculation.h"
#include "refractive_index_functions.h"
int main (int argc, char **argv)
{

 /*
  GOAT::maths::Vector<double> EllPos(0,0,0);
  GOAT::maths::Vector<double> EllR (10,50,50);
  GOAT::raytracing::Ellipsoid Ell(EllPos,EllR, 1.5);
  Ell.setActive (false);
*/
  GOAT::maths::Vector<double> conePos(0,0,0);
  double coneRadius=50;
  GOAT::raytracing::Cone coneObj(conePos, coneRadius, 1, 1.5);
  coneObj.setConeAngle (10.0/180.0*M_PI);
  coneObj.setActive(false); 
  

  GOAT::maths::Vector<double> DetBoxPos(0,0,150);
  GOAT::maths::Vector<double> DetBoxDim(30,30,100);
  GOAT::raytracing::Box DetBox(DetBoxPos,DetBoxDim, 1.0);
  DetBox.setActive (true);
  
   
  
  GOAT::maths::Vector<double> LSPos(0,0,-15);
  double LSD=100.0;
  int LSnumRays=50100;
  double LSwvl=1.0;
//  GOAT::raytracing::LightSrcPlane LS(LSPos, LSnumRays, LSwvl, LSD);
GOAT::raytracing::LightSrcRing_mc LS(LSPos,LSnumRays,LSwvl,0,2);
  LS.setk(GOAT::maths::Vector<double> (0.0,0.0,1.0));
  GOAT::raytracing::Scene S;
  
  S.addObject(&coneObj);
  S.addObject(&DetBox);  
  S.addLightSource (&LS);
  S.setr0 (500.0);
  S.setnS(1.0);
  
  std::vector< std::function< std::complex< double >(double) > > nList;
  nList.push_back (GOAT::raytracing::n_BK7);
  nList.push_back (GOAT::raytracing::n_Vacuum);
  nList.push_back (GOAT::raytracing::n_Vacuum);

  double pulseWidth = 100;
  double refTime = 100.0;
  double spatialRes = 0.4;
  
  GOAT::raytracing::pulseCalculation pc(S);
  pc.setPulseWidth (pulseWidth);
  pc.setSpatialResolution (spatialRes);
  pc.setRefractiveIndexFunctions(nList);
  pc.setSpectralRanges (1);
  pc.setNumWavelengthsPerRange(1);
  
  double time = 400;
  pc.field (time); 
//   GOAT::raytracing::saveFullE(pc.trafo.SAres,"/home/weigel/data/axicon1.dat",0);
   GOAT::raytracing::saveFullE(pc.trafo.SAres,"/home/weigel/data/test.dat",1);

/*
 GOAT::raytracing::Raytrace_Path rp(S);
 rp.setNumReflex(0);
 rp.trace("/home/weigel/data/rays.dat");
*/
}
