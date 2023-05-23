#include "tinyxml2.h"
#include "raytrace.h"
#include <string>

/*****************************************************************//**
 * \file   xml.cpp
 * \brief  Here are a collection of functions to read and write XML-files for GOAT
 *
 * \author Thomas Weigel
 * \date   May 2023
 *********************************************************************/
#include <string>
namespace GOAT
{
	namespace XML
	{
		#define numXMLRootElements   3
		#define XML_NONE         -1 
		#define XML_SCENE_LIGHTSOURCES 0
		#define XML_SCENE_OBJECTS      1
		#define XML_SCENE_DETECTORS    2

		#define numXML_LS_TYPES   4
		#define XML_LS_TYPE_PLANE    0
		#define XML_LS_TYPE_GAUSS    1
		#define XML_LS_TYPE_PLANE_MC 2
		#define XML_LS_TYPE_GAUSS_MC 3


		const std::string sceneXMLElements[] = { "r0","nS","LightSources","Objects","Detectors" };
		const std::string LSXMLAttributes[] = { "Type","Size","Wavelength","NumRays" };
		const std::string LSXMLTYPES[] = {"Plane","Gaussian","Plane_mc","Gaussian_mc"};
	


		class xmlReader
		{
			public:
				void readXML(std::string fname, GOAT::raytracing::Scene& S);

			private:				
	/*			int readLightSrc(ticpp::Element* ell);
				void readLightSrcDefaults(ticpp::Element* ell, double& size, int& numRays, double& wvl, GOAT::maths::Vector<double> &Pos);
				GOAT::maths::Vector<double> readDVector(ticpp::Element* ell);
				int checkSceneElements(std::string);
				int checkLSElements(std::string);

				ticpp::Element* rootElement;
				ticpp::Element* sceneElement;*/
				tinyxml2::XMLElement* rootElement;
				tinyxml2::XMLElement* sceneElement;
				GOAT::raytracing::Scene S;
		};
	}
}
