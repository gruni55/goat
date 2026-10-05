#include "pulsecalculation_kirchhoff.h"
#include "raytrace.h"
#include "lightsrc_mc.h"

int main(int argc, char** argv)
{
	GOAT::raytracing::Scene S;
	S.setr0(5000);
	GOAT::raytracing::Box kBox(GOAT::maths::Vector<double>(0, 0, 25), GOAT::maths::Vector<double>(10,10,500),1.0);
	GOAT::raytracing::DetectorPlane det(GOAT::maths::Vector<double>(0, 0, 0), GOAT::maths::Vector<double>(0, 0, -1),10, 2);
	GOAT::raytracing::Kirchhoff3D k3d(&kBox, 1);	
	GOAT::raytracing::LightSrcRing_mc ls(GOAT::maths::Vector<double>(0, 0, 0),1000000, 1.0, 0, 10);	
	S.setr0(1000);
	S.addLightSource(&ls);
	S.addDetector(&det);
	k3d.addDetector(&det);	
	S.addKirchhoff3D(&k3d);	
	GOAT::raytracing::pulseCalculationKirchhoff pc(S);
	pc.setSpectralRanges(1);
	GOAT::computeSettings settings;
	settings.numThreads = 20;
	std::cout << "test\tr0=" << pc.S.r0 << std::endl;
	pc.setSpatialResolution(1);
	std::cout << "test\tr0=" << pc.S.r0 << std::endl;
	pc.setComputeSettings(settings);
	pc.field(0.0);
	GOAT::raytracing::saveFullE(pc.S.k3D[0]->field3D, "c:\\tmp\\blubb.dat", 0);
	return 0;
}