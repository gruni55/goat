#include "tinyxml2.h"
#include <iostream>
#include "xml.h"





int main(int argc, char** argv)
{
	
	
	GOAT::XML::xmlReader xmlr;
	GOAT::raytracing::Scene S;
	xmlr.readXML("C:\\Users\\thomas\\source\\repos\\goat\\examplesXML\\full.xml");
	
	return 0;
}
