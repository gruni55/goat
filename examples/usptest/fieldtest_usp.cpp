#include "raytrace_field_usp.h"
#include "pulsecalculation_field.h"
#include "refractive_index_functions.h"

int main(int argc, char** argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(1E+3);
	S.setnS(1.0);

	// ------------ Light source --------------
	int numRays = 25000;
	GOAT::raytracing::LightSrcPlane_mc LS(-100 * GOAT::maths::ex, numRays, 1.0, 120.0);
	LS.setk(GOAT::maths::ex);
	S.addLightSource(&LS);

	// --------------- Object : sphere with 50um radius -------------
	GOAT::maths::Vector<double> ellPos;
	GOAT::maths::Vector<double> ellDim(50, 50, 50);
	GOAT::raytracing::Ellipsoid ell(ellPos, ellDim, 1.5);
   S.addObject(&ell);

	// -------------- Box detector ---------------------
	GOAT::maths::Vector<double> boxPos(60, 0, 0);
	GOAT::maths::Vector<double> boxDim(120, 120, 120);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);

        // ------- refractive index list ------
	std::vector< std::function< std::complex< double >(double) > > nList;

	nList.push_back (GOAT::raytracing::n_BK7);
	nList.push_back (GOAT::raytracing::n_Vacuum);


        GOAT::raytracing::pulseCalculation_Field pc(S);
	pc.setPulseWidth(50);
	pc.setSpatialResolution(1);
	pc.setSpectralRanges(200);
	pc.setNumWavelengthsPerRange(1);
	pc.setCenterWavelength(0.5);
	pc.setNumReflex(0);
 	pc.setReferenceTime(400);
	pc.setRefractiveIndexFunctions(nList);
	pc.addBoxDetector(&box);
 	pc.field (450);
        GOAT::raytracing::saveFullE (pc.trafo.SAres,"C:\\users\\weigetz9\\data\\blubb.dat",0);


	
	/*GOAT::raytracing::Raytrace_field_usp rf(S);
	rf.addBoxDetector(&box);
	rf.setResolution(0.25);
	rf.setNumReflex(2);
	rf.trace();
	
	GOAT::raytracing::saveFullE(rf.SE, "C:\\users\\weigetz9\\data\\field_usp.dat");*/
	return 0;
}


