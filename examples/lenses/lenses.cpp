#include "raytrace.h"
#include "raytrace_inel.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation_field.h"
#include "refractive_index_functions.h"

int main(int argc, char** argv)
{
	constexpr double mm = 1000.0;
	constexpr double cm = 10000.0;

	int nRays = 1000000;
	GOAT::raytracing::Scene S;

	GOAT::maths::Vector<double> LSPos(0, 0, -2*cm);
	GOAT::raytracing::LightSrcRing_mc LS(LSPos, nRays, 1.0, 0, 2.0*mm);
	// GOAT::raytracing::LightSrcPlane_mc LS(LSPos, nRays, 1.0,15.0 * mm);
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

	GOAT::maths::Vector<double> detBoxPos(0,0,1*mm);
	GOAT::maths::Vector<double> detBoxDim(1,500,20*mm);
	GOAT::raytracing::Box detBox(detBoxPos,detBoxDim,1);
	detBox.setActive(true);

	GOAT::raytracing::pulseCalculation_Field pf(S);
	pf.setCenterWavelength(0.5);
	std::vector< std::function< std::complex< double >(double) > > nList;
    nList.push_back(GOAT::raytracing::n_BK7);
	nList.push_back(GOAT::raytracing::n_Vacuum);
	pf.setRefractiveIndexFunctions(nList);
    pf.addBoxDetector(&detBox);
	pf.setNumWavelengthsPerRange(1);
	pf.setSpatialResolution(1);
	pf.setSpectralRanges(1);
	pf.setPeriod(2.0*1666.7);
	pf.field(0);
	pf.saveIntensity("/home/weigel/data/linse.dat",0);




	/*GOAT::maths::Vector<double> detPos(0,0,1*mm);
	GOAT::raytracing::DetectorPlane det(detPos,-GOAT::maths::ez,2*mm,500);*/

	/*GOAT::maths::Vector<double> det2Pos(0,0,-1*cm);
	GOAT::raytracing::DetectorPlane det2(det2Pos,-GOAT::maths::ez,make -j
	5*mm,500);*/

/*
	GOAT::raytracing::Raytrace_Path rt(S);
	rt.setNumReflex(0);
	rt.trace("/home/weigel/data/test.dat");

	// S.addDetector (&det);
        // S.addDetector (&det2);

	LS.setNumRays(10000000);
	GOAT::raytracing::Raytrace_pure rp(S);
	rp.setNumReflex(0);
	rp.trace();

	det.save("/home/weigel/data/lensfield_ring.dat");
	// det2.save("/home/weigel/data/lensfield2_ring.dat");


/*

GOAT::maths::Vector<double> boxPos(0, 0, 0.5*cm);
	GOAT::maths::Vector<double> boxDim(200, 200, 150);
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

	pc.setSpectralRanges(200);
	pc.setNumWavelengthsPerRange(10);
	pc.setCenterWavelength(wvl);
	pc.setNumReflex(0);

        LS.setD(1);
        LS.setNumRays(1);  
	double time=pc.findHitTime(1);
        std::cout << "& estimated time=" << time << std::endl;
     
        LS.setD(4.0*mm);
        LS.setNumRays(10000);  
     
        double D;
        std::ofstream oserr("/home/weigel/data/testerr.dat");
        do
        { 
	 D=pc.field(time,GOAT::raytracing::PULSECALCULATION_NOT_CLEAR_RESULT);
	GOAT::raytracing::saveFullE(pc.trafo.SAres, "/home/weigel/data/testlens.dat", 1);
        oserr << D << std::endl;
        std::cout << D << std::endl;

        }  
        while (D>1E-10);

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
