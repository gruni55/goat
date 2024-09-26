#include "raytrace.h"

int main (int argc, char **argv)
{
 int numRays=100;
 
 GOAT::raytracing::Scene S;
 GOAT::maths::Vector<double> LSPos (-100,0,0);
 GOAT::raytracing::LightSrcPlane LS(LSPos,numRays,1.0,20);
 LS.setk (GOAT::maths::ex);
 S.addLightSource(&LS);

 GOAT::maths::Vector<double> mirrorPos (0,0,0);
 GOAT::maths::Vector<double> mirrorDim (1,50,50);
 std::complex<double> nMirror(1.4,0);
 GOAT::raytracing::Box mirror(mirrorPos, mirrorDim, nMirror);
 mirror.setBeta(45.0/180.0*M_PI); 
 
 GOAT::maths::Vector<double> D1Pos(-75, 0, 0);
 GOAT::raytracing::DetectorPlane D1(D1Pos, -GOAT::maths::ex, 20, 100);

 GOAT::maths::Vector<double> D2Pos(0, 0, 100);
 GOAT::raytracing::DetectorPlane D2(D2Pos, -GOAT::maths::ez, 20, 100);

 GOAT::maths::Vector<double> D3Pos(75, 0, 0);
 GOAT::raytracing::DetectorPlane D3(D3Pos, -GOAT::maths::ex, 20, 100);

 S.addDetector(&D1);
 S.addDetector(&D2);
 S.addDetector(&D3);

 S.addObject(&mirror);
 S.setr0(200);

 GOAT::raytracing::Raytrace_Path rp(S);
 rp.setNumReflex(1);
 rp.setShowOutgoingRays(true);
 // rp.trace ("/home/weigel/data/mirror.dat");
 rp.trace("H:\\data\\mirror.dat");
 D1.saveabs("H:\\data\\d1.dat");
 D2.saveabs("H:\\data\\d2.dat");
 D3.saveabs("H:\\data\\d3.dat");
 return 0;
}
