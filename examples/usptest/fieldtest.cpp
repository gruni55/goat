#include "raytrace_field.h"

int main (int argc, char **argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(1E+3);
	S.setnS(1.0);

	int numRays = 1000000;
	GOAT::raytracing::LightSrcPlane_mc LS(-100 * GOAT::maths::ex, numRays, 1.0,100.0);
	LS.setk(GOAT::maths::ex);
	S.addLightSource(&LS);

	GOAT::maths::Vector<double> ellPos;
	GOAT::maths::Vector<double> ellDim(50, 50, 50);
	GOAT::raytracing::Ellipsoid ell(ellPos, ellDim, 1.5);
	S.addObject(&ell);

	GOAT::maths::Vector<double> boxPos;
	GOAT::maths::Vector<double> boxDim(120, 120, 120);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);

	GOAT::raytracing::Raytrace_Field rf(S);
	rf.addBoxDetector(&box);
	rf.setResolution(0.1);
	rf.trace();

	GOAT::raytracing::saveFullE(rf.SE, "C:\\users\\weigetz9\\data\\field.dat");
	return 0;
}
