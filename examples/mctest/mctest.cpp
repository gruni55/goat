#include "raytrace.h"

int main(int argc, char* argv[])
{
    int nRays = 10;
    GOAT::raytracing::Scene S;
    GOAT::maths::Vector<double> P = GOAT::maths::dzero;
    GOAT::maths::Vector<double> r(5, 5, 5);
    GOAT::raytracing::LightSrcGauss_mc LS(-GOAT::maths::ey * 100, nRays, 1.0,0.1,GOAT::maths::dzero,200);
    LS.setk(GOAT::maths::ey);
   
    int nDet = 400;
    GOAT::raytracing::DetectorPlane D(GOAT::maths::Vector<double>(0, 1, 0), GOAT::maths::ey,5, nDet);
    S.addLightSource(&LS);
    S.setnS(1.0);
    S.setr0(500);
    GOAT::raytracing::Raytrace_Path rp;
    rp.setScene(S);
    rp.setNumReflex(1);
    rp.trace("rays.dat");

    S.LS[0]->setNumRays(1000000);
    S.addDetector(&D);
    GOAT::raytracing::Raytrace_pure rp2;
    rp2.setScene(S);
    rp2.setNumReflex(1);
    rp2.trace();
    D.saveabs("test.dat");
    S.Det[0]->saveabs("test2.dat");
    S.Det[0]->savePhase("test2p.dat", 2);
    return 0;
}