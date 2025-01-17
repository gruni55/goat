#include "raytrace.h"

int main (int argc, char **argv)
{
    GOAT::raytracing::Scene S;
    S.setr0(1000);

    GOAT::maths::Vector<double> PosLS(0,0,-10);
    GOAT::raytracing::LightSrcPlane LS(PosLS,2,1.0,2);
    double alpha=std::atof(argv[1])/180.0*M_PI;
    GOAT::maths::Vector<double> k(sin(alpha),0,cos(alpha));

    LS.setk(k);
    GOAT::maths::Vector<std::complex<double> > Pol(1,0,0);
    LS.setPol(Pol);

    S.addLightSource(&LS);
    

    GOAT::maths::Vector<double> PosDet(0,0,0);
    GOAT::raytracing::DetectorPlane det(PosDet,-GOAT::maths::ez,100,10);
    
    S.addDetector(&det);

    GOAT::raytracing::Raytrace_pure rt(S);
    rt.trace();

    S.Det[0]->save("field.dat");

    return 0;
}