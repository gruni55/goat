#include "pulsecalculation_kirchhoff.h"
#include "raytrace.h"
#include "lightsrc_mc.h"

int main(int argc, char** argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(200);
	GOAT::raytracing::Box kBox(GOAT::maths::Vector<double>(10, 10, 10), GOAT::maths::Vector<double>(10,10,10),1.0);
	GOAT::raytracing::DetectorPlane det(GOAT::maths::Vector<double>(0, 0, 0), GOAT::maths::Vector<double>(200, 0, 0), GOAT::maths::Vector<double>(0, 200, 0), 200, 200);
	GOAT::raytracing::Kirchhoff3D k3d(&kBox, 200);
	GOAT::raytracing::LightSrcRing_mc ls(GOAT::maths::Vector<double>(0, 0, 0), 10000, 1.0, 0, 10);
	S.addLightSource(&ls);
		k3d.addDetector(&det);
	S.addKirchhoff3D(&k3d);
	GOAT::raytracing::pulseCalculationKirchhoff pc(S);
	GOAT::computeSettings settings;
	settings.numThreads = 20;
	pc.setComputeSettings(settings);
	pc.field(0.0);
	GOAT::raytracing::saveFullE(pc.S.k3D[0]->field3D, "c:\\tmp\\testKirchhoff.dat", 0);

	return 0;
}