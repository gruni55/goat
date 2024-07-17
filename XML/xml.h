#pragma once
#include "tinyxml2.h"
#include "raytrace.h"
#include <string>
#include "xmltoken.h"
/*****************************************************************//**
 * \file   xml.cpp
 * \brief  Here are a collection of functions to read and write XML-files for GOAT
 *
 * \author Thomas Weigel
 * \date   May 2023
 *********************************************************************/
#include <string>
#include <vector>
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
		 * This function reads the XML file . It can also perform calculations. 
		 * Beware: The parser is case sensitive. Element names starting with capital letters (except of "real" and "imag" for complex numbers), whereas values are always lower case
		 * The structure of the XML-File is as follows (the order of the subdivision may be changed): 
		 * Signs: o: optional, n: necessary, i: independend parameters, i.e. this parameter is valid for all types
		 * Values: float: floating point number
		 *         int: integer number
		 *         bool: possible values: true or false
		 *             
		 * <table>
		 *  <caption id="xml_structure">Structure  of a XML-File</caption>
		 * <tr> <th colspan="5"> Root <th> Value
		 * <tr> <td colspan="5"> Scene 
		 * <tr> <td> <td colspan="4"> Objects 
		 * <tr> <td> <td> <td colspan="3"> Object
		 * <tr> <td> ni <td> <td> <td> Type <td> <td> \"ellipsoid\",\"surface\",\"cone\",<br>\"aspheric_lens\" or \"spheric_lens\" 
		 * <tr> <td> oi <td> <td> <td> Alpha <td> <td> float 
		 * <tr> <td> oi <td> <td> <td> Beta <td> <td> float
		 * <tr> <td> oi <td> <td> <td> Gamma <td> <td> float
		 * <tr> <td> oi <td> <td> <td> IsActive <td> <td> bool
		 * <tr> <td> oi <td> <td> <td> n  
		 * <tr> <td> oin* <td> <td> <td>  <td> real <td> double 
		 * <tr> <td> oin* <td> <td> <td>  <td> imag <td> double
		 * <tr> <td> oi <td> <td> <td> Position
		 * <tr> <td> oin* <td> <td> <td>  <td> x <td> double 
		 * <tr> <td> oin* <td> <td> <td>  <td> y <td> double 
		 * <tr> <td> oin* <td> <td> <td>  <td> z <td> double 
		 * </table>
		*/

		class xmlReader
		{
			public:
				void readXML(const char* fname); ///< function to read the file and store 
				GOAT::raytracing::Scene S; ///< The scene that was read from the file is saved here				

			private:				
                void readScene();
				void readLightSources();
				void readCommands();
				void readObjects();
				void readDetectors();
				void doCalculations();
				GOAT::maths::Vector<double> readVector(tinyxml2::XMLElement* ell, double x = 0, double y = 0, double z = 0);
                GOAT::maths::Vector<double> readVector(tinyxml2::XMLElement* ell, int& xmlError);
                GOAT::maths::Vector<std::complex<double> > readCmplxVector (tinyxml2::XMLElement* ell, int& xmlError);
				std::complex<double> readCmplx(tinyxml2::XMLElement* ell, double defre=0.0, double defim=0.0 );
                std::complex<double> readCmplx(tinyxml2::XMLElement* ell, int& xmlError);
				tinyxml2::XMLNode* rootElement;
				tinyxml2::XMLElement* sceneElement;				
				tinyxml2::XMLElement* calculationElement;
				std::vector<GOAT::raytracing::ObjectShape*> Obj;
				GOAT::raytracing::LightSrc** LS=0;
				int numObj = 0;
				int numLS = 0;

		};
	}
}
