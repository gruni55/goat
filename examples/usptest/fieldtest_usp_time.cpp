#include "raytrace_field_usp.h"
#include "pulsecalculation_field.h"
#include "refractive_index_functions.h"
#include <string>


int main(int argc, char** argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(1E+3);
	S.setnS(1.0);

	// ------------ Light source --------------
	int numRays = 500;
	GOAT::raytracing::LightSrcPlane LS(-100 * GOAT::maths::ex+0* GOAT::maths::ey, numRays, 0.5, 60.0);
	LS.setk(GOAT::maths::ex);
	S.addLightSource(&LS);

	// --------------- Object : sphere with 50um radius -------------
	GOAT::maths::Vector<double> ellPos(-20,10,0);
	GOAT::maths::Vector<double> ellDim(30, 30, 30);
	GOAT::raytracing::Ellipsoid ell(ellPos, ellDim, 1.5);
    S.addObject(&ell);

       	GOAT::maths::Vector<double> boxObjPos(0, 0, 0);
	GOAT::maths::Vector<double> boxObjDim(30, 4, 4);
	GOAT::raytracing::Box boxObj(boxObjPos, boxObjDim, 1.0);
// S.addObject(&boxObj);

	// -------------- Box detector ---------------------
	GOAT::maths::Vector<double> boxPos(0, 0, 0);
	GOAT::maths::Vector<double> boxDim(100, 60, 1);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);

        // ------- refractive index list ------
	std::vector< std::function< std::complex< double >(double) > > nList;

	nList.push_back (GOAT::raytracing::n_BK7);
	nList.push_back (GOAT::raytracing::n_Vacuum);


        GOAT::raytracing::pulseCalculation_Field pc(S);
	pc.setPulseWidth(50);
	pc.setSpatialResolution(1);
	pc.setSpectralRanges(100);
	pc.setNumWavelengthsPerRange(1);
	pc.setCenterWavelength(0.5);
	pc.setNumReflex(2);
 	pc.setReferenceTime(0);
	pc.setRefractiveIndexFunctions(nList);
	pc.addBoxDetector(&box);

        std::string fname;
        std::string prefix="/home/weigel/data/field_uspr_";
        double dt=5;
        double time;
/*        pc.field (450); 
        pc.saveIntensity("/home/weigel/data/field_usp.dat", 0);*/

        for (int i=0; i<=100; i++)
        {
        time=0+i*dt;
        std::cout << "current time: " << time << "fs" << std::endl;
 	pc.field (time);
        fname=prefix + std::to_string(0+i*5) + ".dat"; 
        std::cout << "save intensity to file:" <<  fname << std::endl;
	pc.saveIntensity(fname, 0);
        }
	return 0;
}


