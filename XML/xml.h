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
	
		/**
		 * @brief This class provides functionality to read a XML-file 
		 * This function reads the XML file . It can also make calculations.
		 */

		class xmlReader
		{
			public:
				void readXML(std::string fname); ///< function to read the file
				GOAT::raytracing::Scene S; ///< The scene that was read from the file is saved here

			private:				
				GOAT::maths::Vector<double> readVector(tinyxml2::XMLElement* ell, double x = 0, double y = 0, double z = 0);
				std::complex<double> readCmplx(tinyxml2::XMLElement* ell, double defre=0.0, double defim=0.0 );
				GOAT::maths::Vector<double> readVector(tinyxml2::XMLElement* ell, int& xmlError);				
				tinyxml2::XMLElement* rootElement;
				tinyxml2::XMLElement* sceneElement;				
		};
	}
}
