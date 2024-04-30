#include "raytrace.h"

int main()
{
  GOAT::raytracing::Scene S;
  
  GOAT::maths::Vector<double> conePos (0,0,0);
  double coneRadius = 20;
  double coneHeight = 10;
  GOAT::raytracing::Cone coneObj(conePos, coneRadius, coneHeight,1.5);
  S.addObject(&coneObj);
  
  int numRays=10000000;
  GOAT::maths::Vector<double> LSPos (0,0,-20);
  GOAT::raytracing::LightSrcPlane_mc LS(LSPos, numRays, 1, 40);
  LS.setk(GOAT::maths::ez);


  S.addLightSource (&LS);
  S.setRaytype(LIGHTSRC_RAYTYPE_RAY);

  S.setr0(100);
  S.setnS(1.0);

  double width=50;
  int numCells=301;
  GOAT::maths::Vector<double> detPos (0,0,30);
  GOAT::raytracing::DetectorPlane det(detPos,GOAT::maths::ez,width,numCells);
  S.addDetector (&det); 
  
  GOAT::raytracing::Raytrace_pure rp (S);
  rp.trace();
  // S.Det[0]->saveabs("/home/weigel/data/testf1.dat");
  S.Det[0]->save("C:\\Users\\Thomas\\data\\testf.dat");
  

/*
  GOAT::raytracing::Raytrace_Path rp;
        rp.setShowOutgoingRays(true);
        rp.setNumReflex(0);
        rp.setScene(S);
        rp.trace("/home/weigel/data/test.dat");
  */
  return 0;
 }
