#include "raytrace.h"
#include "asphericLens.h"

int main(int argc, char** argv)
{
	int nRays = 100;
	GOAT::raytracing::Scene S;

	GOAT::maths::Vector<double> LSPos(0, 0, -500);
	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, nRays, 1.0,5.0);
	LS.setk(GOAT::maths::ez);

	S.setr0(1E+5);
	S.addLightSource(&LS);

	GOAT::raytracing::asphericLensParms lensParms;
	/*
	lensParms.left.A = {0,0,0,0,-8.168899E-3,0,-1.995690E-4,0,2.200933E-4,0,-3.982084E-5,0,2.656612E-6};
	lensParms.left.k = -0.431729;
	lensParms.left.R = 3.879161;

	lensParms.right.A = { 0,0,0,0,0.008168899,0.000199569,-0.0002200933, 3.98208400000000e-05, -2.656612e-06 };
	lensParms.right.k= -0.431729;
	lensParms.right.R = -3.879161;
	*/

	lensParms.left.A = { 0,0,0,0,-8.168899E-3,0,-1.995690E-4,0,2.200933E-4,0,-3.982084E-5,0,2.656612E-6 };
	lensParms.left.k = -0.431729;
	lensParms.left.R = 3.879161;

	lensParms.right.A = { 0,0,0,0,0.008168899,0.000199569,-0.0002200933, 3.98208400000000e-05, -2.656612e-06 };
	lensParms.right.k = -0.431729;
	lensParms.right.R = -3.879161;

	lensParms.radius = 2.25;
	lensParms.offset = 1.13;

	GOAT::maths::Vector<double> lensPos;
	GOAT::raytracing::asphericLens aLens(lensPos,1.5,lensParms);
	S.addObject(&aLens);
	GOAT::raytracing::Raytrace_Path rt(S);
	rt.trace("C:\\Users\\weigetz9\\data\\path.dat");

	return 0;
}