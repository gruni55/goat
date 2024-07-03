#include "tinyxml2.h"
#include <iostream>
#include "xml.h"
#include "raytrace.h"





int main(int argc, char** argv)
{
	
	
	GOAT::XML::xmlReader xmlr;
	GOAT::raytracing::Scene S;
        xmlr.readXML (argv[1]);
// 	xmlr.readXML("H:\\github\\goat\\examplesXML\\full.xml");
/*	GOAT::raytracing::Raytrace_Path rp(xmlr.S);
	rp.trace("H:\\data\\path.dat"); */
	return 0;
}
