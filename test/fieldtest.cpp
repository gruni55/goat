#include "raytrace_inel.h"

int main (int argc, char **argv)
{
    GOAT::raytracing::Scene S;

    int nRays=1000;
    GOAT::maths::Vector<double> LSPos(-100,0,0);
    GOAT::raytracing::LightSrcPlane LS(LSPos,nRays,0.5,100.0);
    LS.setk(GOAT::maths::ex);

    S.addLightSource(&LS);

    GOAT::maths::Vector<double> ellPos (0,0,0);
    GOAT::maths::Vector<double> ellDim (50,50,50);
    GOAT::raytracing::Ellipsoid ell(ellPos,ellDim,1.5);
    ell.setActive(true);
    S.addObject(&ell);

    S.setr0(250);

    GOAT::raytracing::RRTParms parms;
    parms.e1=GOAT::maths::ex;
    parms.e2=GOAT::maths::ey;

    GOAT::raytracing::Raytrace_Inel rt(S,1000);
    rt.setExcitationFieldOnly();
    rt.trace(parms);

    rt.exportExcitation ("/home/weigel/data/test");
  return 0;
}
