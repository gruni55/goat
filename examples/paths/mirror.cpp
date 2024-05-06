#include "raytrace.h"

int main (int argc, char **argv)
{
 int numRays=10;
 
 GOAT::raytracing::Scene S;
 GOAT::maths::Vector<double> LSPos (-100,0,0);
 GOAT::raytracing::LightSrcPlane LS(LSPos,numRays,1.0,20);
 LS.setk (GOAT::maths::ex);
 S.addLightSource(&LS);

 GOAT::maths::Vector<double> mirrorPos (0,0,0);
 GOAT::maths::Vector<double> mirrorDim (1,50,50);
 std::complex<double> nMirror(0.04,7.1155);
 GOAT::raytracing::Box mirror(mirrorPos, mirrorDim, nMirror);
 mirror.setBeta(45.0/180.0*M_PI); 
 
 S.addObject(&mirror);
 S.setr0(200);
 GOAT::raytracing::Raytrace_Path rp(S);
 rp.setNumReflex(1);
 rp.setShowOutgoingRays(true);
 rp.trace ("/home/weigel/data/mirror.dat");
 
 return 0;
}
