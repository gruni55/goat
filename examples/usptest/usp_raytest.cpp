#include "raytrace.h"
#include "cone.h"


int main(int argc, char** argv)
{

	GOAT::maths::Vector<double> conePos(0, 0, 0);
	double coneHeight = 3;
	double coneRadius = 50;
	GOAT::raytracing::Cone ConeObj(conePos, coneRadius, coneHeight, 1.5);
	ConeObj.setConeAngle(10.0 / 180.0 * M_PI);
	// ConeObj.setBeta(M_PI_2);


	
	GOAT::maths::Vector<double> LSPos(0, 0, -10);
	double LSD =10.0;
	int LSnumRays = 10000000;
	double LSwvl = 1.0;
	// GOAT::raytracing::LightSrcPlane_mc LS(LSPos, LSnumRays, LSwvl, LSD);

	GOAT::raytracing::LightSrcRing_mc LS(LSPos, LSnumRays, LSwvl, 0, 50.0);
	
	LS.setk(GOAT::maths::Vector<double>(0.0, 0.0, 1.0));

	GOAT::maths::Vector<double> detPos(0, 0, 1000);
	GOAT::raytracing::DetectorPlane dp(detPos, GOAT::maths::ez, 100, 200);


	GOAT::raytracing::Scene S; 
	S.addObject(&ConeObj);
	
	//S.addObject(&Ell);	
	//S.addObject(&DetBox);
	//S.setRaytype(LIGHTSRC_RAYTYPE_RAY);
	S.addDetector(&dp);
	S.addLightSource(&LS);
	S.setr0(2000.0);
	S.setnS(1.0);
	
	// GOAT::raytracing::Raytrace_Path rp(S);
	 GOAT::raytracing::Raytrace_pure rp(S);
	rp.setNumReflex(0);
	rp.trace();
	 dp.save("H:\\data\\test.dat");
	// rp.trace("/home/weigel/data/rays.dat");
	// rp.trace("H:\\data\\rays.dat");


	return 0;
}
