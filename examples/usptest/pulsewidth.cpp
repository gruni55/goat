#include "pulsecalculation.h"
#include "refractive_index_functions.h"

int main (int argc, char **argv)
{
    double m=1E+6;
    double cm=1E+4;
    double mm=1E+3;

    GOAT::maths::Vector<double> LSPos(0,0,0);
    GOAT::raytracing::LightSrcPlane LS(LSPos, 1,1E-5,1);
    LS.setk(GOAT::maths::ex);


    GOAT::maths::Vector<double> objPos(1.0*m,0,0);
    GOAT::maths::Vector<double> objDim(10.0*cm,10,10);
    GOAT::raytracing::Box obj(objPos,objDim,1.5);
    obj.setActive(false);

    int n=3000;
    GOAT::maths::Vector<double> detPos(2.0*m,0,0);
    GOAT::maths::Vector<double> detDim(n,1,1);
    GOAT::raytracing::Box det(detPos,detDim,1.5);
    det.setActive(true);

    GOAT::raytracing::Scene S;
    S.addLightSource(&LS);
    S.addObject(&obj);
    S.addObject(&det);
    S.setr0(10.0E+6);

    std::vector< std::function< std::complex< double >(double) > > nList;
    nList.push_back (GOAT::raytracing::n_BK7);
    nList.push_back (GOAT::raytracing::n_Vacuum);
    nList.push_back (GOAT::raytracing::n_Vacuum);

    GOAT::raytracing::pulseCalculation pc(S);
    double wvl=0.5;
    double spatialResolution = 0.5;

    pc.setNumWavelengthsPerRange(1);
    pc.setNumReflex(0);
    pc.setSpectralRanges(500);
    pc.setSpatialResolution(spatialResolution);
    pc.setPulseWidth(50);
    pc.setRefractiveIndexFunctions(nList);
    pc.setCenterWavelength(wvl);
    std::vector<double> d(n);
    double fwhms, fwhm;
    std::size_t maxIndex;
    double time;
    std::ofstream os("/home/weigel/data/test3.dat");

    std::string prefix="/home/weigel/data/pulsewidth/test_";
    std::string fname;
    std::ofstream osh;
    int i=0;
     wvl=0.5;
    for (double wvl=1.57; wvl<1.7; wvl+=0.005)
    {
    pc.setCenterWavelength(wvl);
    time=pc.findHitTime(1);
    pc.setReferenceTime(time-500);
    pc.field(time);
    //fname=prefix + std::to_string(i) + ".dat";
    //osh.open(fname);
    for (int i=0; i<n; i++)
    {
        d[i]=abs2(pc.trafo.SAres.G[1][i][2][2]);
    //    osh << d[i] << std::endl;
    }
   // osh.close();
    i++;
    GOAT::maths::findmax(d, maxIndex);
      fwhms = GOAT::maths::FWHM(d, maxIndex);
      fwhm = fwhms * spatialResolution / 0.3;
    os << time << "\t" << wvl << "\t" << pc.dWvl << "\t" << fwhm << std::endl;
    // saveFullE(pc.trafo.SAres,"/home/weigel/data/feld2.dat",1);
    }
    os.close();
    return 0;
}
