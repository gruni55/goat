#include "raytrace_inel.h"
#include "pulsecalculation.h"
#include "refractive_index_functions.h"

int main (int argc, char **argv)
{
    GOAT::raytracing::Scene S;

    int nRays=20000;
    GOAT::maths::Vector<double> LSPos(-100,0,0);
    // GOAT::raytracing::LightSrcPlane_mc LS(LSPos,nRays,0.5,100.0);
    GOAT::raytracing::LightSrcGauss_mc LS(LSPos, nRays, 10, 40, GOAT::maths::dzero, 100);
    
    S.addLightSource(&LS);

    /*GOAT::maths::Vector<double> ellPos(0, 0, 0);
    GOAT::maths::Vector<double> ellDim (50,50,50);
    GOAT::raytracing::Ellipsoid ell(ellPos,ellDim,1.5);
    ell.setActive(true);
    S.addObject(&ell);*/

    GOAT::maths::Vector<double> boxPos(0, 0, 0);
    GOAT::maths::Vector<double> boxDim(100, 100, 100);
    GOAT::raytracing::Box box(boxPos, boxDim, 1.0);
    S.addObject(&box);
    S.setr0(250);
    std::vector< std::function< std::complex< double >(double) > > nList;
    
    nList.push_back(GOAT::raytracing::n_Vacuum);
    nList.push_back(GOAT::raytracing::n_Vacuum);

    // ----------- parameters for pulse calculation ------------
   /* double pulseWidth = 100;
    double refTime = 1000;
    double spatialRes = 0.25;
    double wvl = 0.5;

    GOAT::raytracing::pulseCalculation pc(S);
    pc.setPulseWidth(pulseWidth);
    pc.setSpatialResolution(spatialRes);
    pc.setRefractiveIndexFunctions(nList);

    pc.setSpectralRanges(50);
    pc.setNumWavelengthsPerRange(100);
    pc.setCenterWavelength(wvl);
    pc.setNumReflex(0);

    double time = pc.findHitTime(0) + 100;
    pc.field(time);
    GOAT::raytracing::saveFullE(pc.trafo.SAres, "c:\\users\\weigetz9\\data\\testbox.dat", 1);
    */
    
    GOAT::raytracing::RRTParms parms;
    parms.e1=GOAT::maths::ex;
    parms.e2=GOAT::maths::ey;

    GOAT::raytracing::Raytrace_Inel rt(S,500);
    rt.setExcitationFieldOnly();
    rt.trace(parms);

    rt.exportExcitation ("C:\\users\\weigetz9\\data\\blubb",GOAT::raytracing::INEL_EXPORT_EXCITATION_FIELD_VECTOR);
  return 0;
}
