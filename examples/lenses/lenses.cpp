#include "raytrace.h"
#include "raytrace_inel.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation.h"
#include "refractive_index_functions.h"

int main(int argc, char** argv)
{
	constexpr double mm = 1000.0;
	constexpr double cm = 10000.0;

	int nRays = 50;
	GOAT::raytracing::Scene S;

	GOAT::maths::Vector<double> LSPos(0, 0, -2*cm);
	// GOAT::raytracing::LightSrcRing_mc LS(LSPos, nRays, 1.0, 0, 5);
	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, nRays, 1.0,4.0 * mm);
	LS.setk(GOAT::maths::ez);	
	S.addLightSource(&LS);
	S.setr0(10*cm);
	S.setnS(1.0);

	

	GOAT::raytracing::lensParms lensParms;
	lensParms.left.R = 5;
	lensParms.left.curvature = GOAT::raytracing::flat;

	lensParms.right.R = 5.2 * mm;
	lensParms.right.curvature = GOAT::raytracing::convex;
	
	lensParms.offset = 1.5 * mm;
	lensParms.radius = 3.0 * mm;

	GOAT::maths::Vector<double> lensPos;
	GOAT::raytracing::sphericLens Lens(lensPos,1.5,lensParms);
	Lens.setActive(false);
	S.addObject(&Lens);

	GOAT::raytracing::Raytrace_Path rt(S);
	rt.setNumReflex(0);
	rt.trace("h:\\data\\test.dat");




	GOAT::maths::Vector<double> boxPos(0, 0, 2*cm);
	GOAT::maths::Vector<double> boxDim(2*mm, 2*mm, 150);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);
	box.setActive(true);

	
	S.addObject(&box);
	S.setr0(500);

	std::vector< std::function< std::complex< double >(double) > > nList;
    nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_Vacuum);
	nList.push_back(GOAT::raytracing::n_Vacuum);

	// ----------- parameters for pulse calculation ------------
	double pulseWidth = 100;
	double refTime = 1000;
	double spatialRes = 1;
	double wvl = 0.5;

	GOAT::raytracing::pulseCalculation pc(S);
	pc.setPulseWidth(pulseWidth);
	pc.setSpatialResolution(spatialRes);
	pc.setRefractiveIndexFunctions(nList);

	pc.setSpectralRanges(500);
	pc.setNumWavelengthsPerRange(10);
	pc.setCenterWavelength(wvl);
	pc.setNumReflex(0);

	double time=pc.findHitTime(1)+100;
	pc.field(time);
	GOAT::raytracing::saveFullE(pc.trafo.SAres, "/home/weigel/data/testlens.dat", 1);

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
