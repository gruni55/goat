#include "raytrace.h"
#include "cone.h"


int main(int argc, char** argv)
{
	/*
	GOAT::maths::Vector<double> EllPos(0, 0, 0);
	GOAT::maths::Vector<double> EllR(10, 50, 50);
	GOAT::raytracing::Ellipsoid Ell(EllPos, EllR, 1.5);
	*/

	GOAT::maths::Vector<double> conePos(0, 0, 10);
	double coneHeight = 3;
	double coneRadius = 10;
	GOAT::raytracing::Cone ConeObj(conePos, coneRadius, coneHeight, 1.5);
	ConeObj.setConeAngle(5.0 / 180.0 * M_PI);

	GOAT::maths::Vector<double> DetBoxPos(60, 0, 0);
	GOAT::maths::Vector<double> DetBoxDim(100, 50, 50);
	GOAT::raytracing::Box DetBox(DetBoxPos, DetBoxDim, 1.0);

	
	GOAT::maths::Vector<double> LSPos(0, 0, -15);
	double LSD =10.0;
	int LSnumRays = 500;
	double LSwvl = 1.0;
	// GOAT::raytracing::LightSrcPlane_mc LS(LSPos, LSnumRays, LSwvl, LSD);
	GOAT::raytracing::LightSrcRing_mc LS(LSPos, LSnumRays, LSwvl, 5.0, 6.0);
	LS.setk(GOAT::maths::Vector<double>(0.0, 0.0, 1.0));
	GOAT::raytracing::Scene S;
	S.addObject(&ConeObj);

	//S.addObject(&Ell);	
	// S.addObject(&DetBox);
	//S.setRaytype(LIGHTSRC_RAYTYPE_RAY);
	
	S.addLightSource(&LS);
	S.setr0(50.0);
	S.setnS(1.0);
	
	GOAT::raytracing::Raytrace_Path rp(S);
	rp.setNumReflex(0);
	// rp.trace("/home/weigel/data/rays.dat");
	rp.trace("H:\\data\\rays.dat");


	return 0;
}
