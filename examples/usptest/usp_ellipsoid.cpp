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
  

  GOAT::maths::Vector<double> DetBoxPos(0,0,100);
  GOAT::maths::Vector<double> DetBoxDim(30,30,100);
  GOAT::raytracing::Box DetBox(DetBoxPos,DetBoxDim, 1.0);
  DetBox.setActive (true);
  
   
  
  GOAT::maths::Vector<double> LSPos(0,0,-15);
  double LSD=100.0;
  int LSnumRays=600000.0;
  double LSwvl=1.0;
  GOAT::raytracing::LightSrcPlane_mc LS(LSPos, LSnumRays, LSwvl, LSD);
  LS.setk(GOAT::maths::Vector<double> (0.0,0.0,1.0));
  GOAT::raytracing::Scene S;
  
  S.addObject(&coneObj);
  S.addObject(&DetBox);  
  S.addLightSource (&LS);
  S.setr0 (200.0);
  S.setnS(1.0);
  
  std::vector< std::function< std::complex< double >(double) > > nList;
  nList.push_back (GOAT::raytracing::n_BK7);
  nList.push_back (GOAT::raytracing::n_Vacuum);
  nList.push_back (GOAT::raytracing::n_Vacuum);

  double pulseWidth = 100.0;
  double refTime = 0.0;
  double spatialRes = 0.2;
  
  GOAT::raytracing::pulseCalculation pc(S);
  pc.setPulseWidth (pulseWidth);
  pc.setSpatialResolution (spatialRes);
  pc.setRefractiveIndexFunctions(nList);
  
  double time = 300+135;
  pc.field (time); 
//   GOAT::raytracing::saveFullE(pc.trafo.SAres,"/home/weigel/data/ellipsoid.dat",0);
   GOAT::raytracing::saveFullE(pc.trafo.SAres,"/home/weigel/data/test1.dat",1);

/*
 GOAT::raytracing::Raytrace_Path rp(S);
 rp.setNumReflex(0);
 rp.trace("/home/weigel/data/rays.dat");
*/
}
