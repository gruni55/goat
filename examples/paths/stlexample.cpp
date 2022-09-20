/*
*  This example shows you how to import an object from a STL-file in a simple ray path log. 
*  A gaussian light source is used
*/
#include "raytrace.h"

int main()
{
	GOAT::raytracing::surface O(GOAT::maths::dzero, 1.5);
	// O.importBinSTL("rotor_30Grad.stl");
	O.importBinSTL("utahteapot.stl");
	O.setn(1.5);
    // O.scale(20.0);
	
	// LightSrcGauss LS(-200*ez, 20, 1.0, 1.0, dzero);
	GOAT::raytracing::LightSrcPlane LS = GOAT::raytracing::LightSrcPlane(-15 * GOAT::maths::ex+2.5* GOAT::maths::ez, 5, 1.0, 3);
	GOAT::raytracing::Scene S;
	S.setr0(20.0);
	LS.setk(GOAT::maths::ex);
	// S.setnS(1.33);
	S.addObject(&O);
	S.addLightSource(&LS);
	

//	LS.setPos(-20 * ez - 10 * ex);
//	LS.setFocuspos(-3 * ex);
//	LS.setNA(1.2);

	GOAT::raytracing::Raytrace_Path rp;
	rp.setShowOutgoingRays(true);
	rp.setNumReflex(0);
	rp.setScene(S);
	rp.trace("test.dat");
}
