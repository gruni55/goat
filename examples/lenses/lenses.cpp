#include "raytrace.h"
#include "raytrace_inel.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation.h"
#include "refractive_index_functions.h"

int main(int argc, char** argv)
{
	int nRays = 100;
	GOAT::raytracing::Scene S;

	GOAT::maths::Vector<double> LSPos(0, 0, -50);
	// GOAT::raytracing::LightSrcRing_mc LS(LSPos, nRays, 1.0, 0, 5);
	 GOAT::raytracing::LightSrcPlane LS(LSPos, nRays, 1.0,20.0);
	LS.setk(GOAT::maths::ez);

	
	S.addLightSource(&LS);

	GOAT::raytracing::lensParms lensParms;
	lensParms.left.R = 30;
	lensParms.left.curvature = GOAT::raytracing::flat;

	lensParms.right.R = 15;
	lensParms.right.curvature = GOAT::raytracing::convex;
	
	lensParms.offset = 0.0;
	lensParms.radius = 20.0;

	GOAT::maths::Vector<double> lensPos;
	GOAT::raytracing::sphericLens Lens(lensPos,1.5,lensParms);
	Lens.setActive(false);

	GOAT::maths::Vector<double> boxPos(0, 0, 110);
	GOAT::maths::Vector<double> boxDim(40, 40, 150);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);
	box.setActive(true);
	S.addObject(&box);
	S.setr0(500);

	std::vector< std::function< std::complex< double >(double) > > nList;
	nList.push_back(GOAT::raytracing::n_Vacuum);
	nList.push_back(GOAT::raytracing::n_Vacuum);

	// ----------- parameters for pulse calculation ------------
	double pulseWidth = 100;
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

	double time=pc.findHitTime(0)+100;
	pc.field(time);
	GOAT::raytracing::saveFullE(pc.trafo.SAres, "h:\\data\\test.dat", 0);

	/*
	GOAT::raytracing::Raytrace_Inel rt(S, 1000);
	rt.setExcitationFieldOnly();
	GOAT::raytracing::RRTParms parms;
	
	rt.trace(parms);
	rt.exportExcitation("h:\\data\\test", GOAT::raytracing::INEL_EXPORT_EXCITATION_FIELD_ABS);
	*/
	/*GOAT::raytracing::Raytrace_Path rt(S);
	rt.setNumReflex(0);
	rt.trace("H:\\data\\path.dat");
	*/

	return 0;
}