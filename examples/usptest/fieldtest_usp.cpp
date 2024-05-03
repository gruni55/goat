#include "raytrace_field_usp.h"
#include "pulsecalculation_field.h"
#include "refractive_index_functions.h"



int main(int argc, char** argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(1E+3);
	S.setnS(1.0);

	// ------------ Light source --------------
	int numRays = 20000;
	GOAT::raytracing::LightSrcPlane_mc LS(-100 * GOAT::maths::ex+0* GOAT::maths::ey, numRays, 0.5, 60.0);
	LS.setk(GOAT::maths::ex);
	S.addLightSource(&LS);

	// --------------- Object : sphere with 50um radius -------------
	GOAT::maths::Vector<double> ellPos;
	GOAT::maths::Vector<double> ellDim(30, 30, 30);
	GOAT::raytracing::Ellipsoid ell(ellPos, ellDim, 1.5);
   S.addObject(&ell);

       	GOAT::maths::Vector<double> boxObjPos(0, 0, 0);
	GOAT::maths::Vector<double> boxObjDim(30, 4, 4);
	GOAT::raytracing::Box boxObj(boxObjPos, boxObjDim, 1.0);
// S.addObject(&boxObj);

	// -------------- Box detector ---------------------
	GOAT::maths::Vector<double> boxPos(0, 0, 0);
	GOAT::maths::Vector<double> boxDim(100, 10, 10);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);

        // ------- refractive index list ------
	std::vector< std::function< std::complex< double >(double) > > nList;

	nList.push_back (GOAT::raytracing::n_BK7);
	nList.push_back (GOAT::raytracing::n_Vacuum);


        GOAT::raytracing::pulseCalculation_Field pc(S);
	pc.setPulseWidth(50);
	pc.setSpatialResolution(0.25);
	pc.setSpectralRanges(500);
	pc.setNumWavelengthsPerRange(1);
	pc.setCenterWavelength(0.5);
	pc.setNumReflex(0);
 	pc.setReferenceTime(0);
	pc.setRefractiveIndexFunctions(nList);
	pc.addBoxDetector(&box);

 	pc.field (200);
	pc.saveIntensity("h:\\data\\field_usp.dat", 0);
    //    GOAT::raytracing::saveFullE (pc.trafo.SAres,"H:\\data\\blubb.dat",0);
	


//	GOAT::raytracing::saveFullE (pc.trafo.SAres,"/home/weigel/data/blubb.dat",0);
	save (pc.rt.SA[0],"C:\\users\\weigetz9\\data\\test.dat");

	
	/*GOAT::raytracing::Raytrace_field_usp rf(S);
	rf.addBoxDetector(&box);
	rf.setResolution(0.25);
	rf.setNumReflex(2);
	rf.trace();
	
	GOAT::raytracing::saveFullE(rf.SE, "C:\\users\\weigetz9\\data\\field_usp.dat");*/
	return 0;
}


