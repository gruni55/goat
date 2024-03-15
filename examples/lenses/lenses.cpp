#include "raytrace.h"
#include "raytrace_inel.h"
#include "lens.h"
#include "sphericLens.h"

int main(int argc, char** argv)
{
	int nRays = 20000;
	GOAT::raytracing::Scene S;

	GOAT::maths::Vector<double> LSPos(0, 0, -50);
	GOAT::raytracing::LightSrcRing LS(LSPos, nRays, 1.0, 0, 14.5);
	// GOAT::raytracing::LightSrcPlane LS(LSPos, nRays, 1.0,20.0);
	LS.setk(GOAT::maths::ez);

	S.setr0(100);
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
	S.addObject(&Lens);

	GOAT::maths::Vector<double> boxPos(0, 0, 110);
	GOAT::maths::Vector<double> boxDim(1, 150, 150);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);
	box.setActive(true);

	GOAT::raytracing::Raytrace_Inel rt(S, 1000);
	rt.setExcitationFieldOnly();
	GOAT::raytracing::RRTParms parms;
	
	rt.trace(parms);
	rt.exportExcitation("h:\\data\\test.dat");

	/*GOAT::raytracing::Raytrace_Path rt(S);
	rt.setNumReflex(0);
	rt.trace("H:\\data\\path.dat");
	*/

	return 0;
}