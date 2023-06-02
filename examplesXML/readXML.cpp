#include "tinyxml2.h"
#include <iostream>
#include "xml.h"





int main(int argc, char** argv)
{
	/*ticpp::Document doc("C:\\Users\\weigetz9\\source\\repos\\goat\\examplesXML\\full.xml");
	doc.LoadFile();
	ticpp::Element* rootElement = doc.FirstChildElement();
	ticpp::Element* ell = rootElement->FirstChildElement();	
	ticpp::Element* ells = ell->FirstChildElement();
	ells = ells->NextSiblingElement();
	std::cout << "The Name of ells:" << ells->Value() << std::endl;
	doc.Clear();
	*/

	
	GOAT::XML::xmlReader xmlr;
	GOAT::raytracing::Scene S;
	xmlr.readXML("C:\\Users\\weigetz9\\source\\repos\\goat\\examplesXML\\full.xml");
	
/*
	GOAT::maths::Vector<double> Pos;
	std::complex<double> n;
	GOAT::raytracing::surface* S1 = new GOAT::raytracing::surface(Pos, n);
	GOAT::raytracing::surface* S2 = new GOAT::raytracing::surface(Pos, n);*/
	return 0;
}
