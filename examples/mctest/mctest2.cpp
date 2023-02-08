#include "raytrace_inel.h"
#include <chrono>



int main(int argc, char* argv[])
{

    

    auto start = std::chrono::steady_clock::now();    

    int nRays = 10;
    GOAT::raytracing::Scene S;
    GOAT::maths::Vector<double> P = GOAT::maths::dzero;
    GOAT::maths::Vector<std::complex<double> > Pol(1.0,0.0,0.0);
    
   // GOAT::raytracing::LightSrcGauss_mc LS(-GOAT::maths::ey * 100, nRays, 1.0,0.1,GOAT::maths::dzero,200);
   GOAT::raytracing::LightSrcPlane_mc LS(-GOAT::maths::ey * 100,nRays,1.0,40.0,Pol);
    LS.setk(GOAT::maths::ey);
   
    GOAT::raytracing::Ellipsoid O(GOAT::maths::Vector<double>(0,0,0),GOAT::maths::Vector<double>(20,20,20),1.5);
    S.addObject(&O);
    S.addLightSource(&LS);
    S.setnS(1.0);
    S.setr0(500);
    S.LS[0]->setNumRays(100000);
    GOAT::raytracing::Raytrace_Inel ri(S,2500);    
    ri.setNumReflex(2);
    ri.setExcitationFieldOnly();
    GOAT::raytracing::RRTParms rrtparms;
    ri.trace(rrtparms);
    auto end = std::chrono::steady_clock::now();
    double us=std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    double seconds = floor(us / 1E+6);
    double ms= floor((us - seconds*1E+6) / 1000.0);
    double restus = us - seconds*1E+6 - ms * 1000.0;
    std::cout << "Calculation time: "<< std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << "µs" ;
    std::cout << "= " << seconds << "s :" << ms << "ms : " << restus << "µs" << std::endl;
  /*  GOAT::raytracing::saveabsE(*ri.SGE, "test.dat");
    GOAT::raytracing::saveExPhase(*ri.SGE,"testp.dat");*/
    return 0;
}
