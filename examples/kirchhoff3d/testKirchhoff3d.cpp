#include "tinyxml2.h"
#include <iostream>
#include "xml.h"
#include "raytrace.h"
#include "kirchhoff.h"

int main(int argc, char** argv)
{
	GOAT::XML::xmlReader xmlr;	
	xmlr.readXML(argv[1]);
	GOAT::raytracing::Raytrace_pure rp(xmlr.S);
	rp.trace();
	GOAT::raytracing::Kirchhoff3D k3d((GOAT::raytracing::Box *)xmlr.S.Obj[0],xmlr.S.NumCellsPerDir);
	for (auto det : xmlr.S.Det)
		k3d.addDetector((GOAT::raytracing::DetectorPlane*)det);
	k3d.calc(1.0,10);
	auto SA=k3d.field(); 
	//SA.write("c:\\tmp\\testKirchhoff3DField.dat");
	GOAT::raytracing::saveFullE(k3d.field3D, "c:\\tmp\\testKirchhoff3DField.dat", 0);
	return 0;
}
