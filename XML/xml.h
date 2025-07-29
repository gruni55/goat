#pragma once
#include "tinyxml2.h"
#include "raytrace.h"
#include <string>
#include "xmltoken.h"
/*****************************************************************//**
 * \file   xml.cpp
 * \brief  Here are a collection of functions to read and write XML-files for GOAT
 * This part provides a collection of classes, which enables to describe the scene and the calculations.
 * The xmlreader class enables to read the XML file.
 * \author Thomas Weigel
 * \date   May 2023
 *********************************************************************/
#include <string>
#include <vector>
#include <sstream>
#include <locale>
#include <iomanip>
namespace GOAT
{
	namespace XML
	{
		#define numXMLRootElements   3
		#define XML_NONE         -1 
		#define	XML_SCENE_R0		   0
		#define	XML_SCENE_NS 1
		#define XML_SCENE_NCELLS 2
		#define XML_SCENE_LIGHTSOURCES 3
		#define XML_SCENE_OBJECTS      4
		#define XML_SCENE_DETECTORS    5

		#define numXML_LS_TYPES   4
		#define XML_LS_TYPE_PLANE    0
		#define XML_LS_TYPE_GAUSS    1
		#define XML_LS_TYPE_PLANE_MC 2
		#define XML_LS_TYPE_GAUSS_MC 3




		const std::string sceneXMLElements[] = { "r0","ns","CellsPerDir","lightsources","objects","detectors"};
		const std::string LSXMLAttributes[] = { "type","size","wavelength","numrays","numraysRT"};
		const std::string LSXMLTYPES[] = {"plane","gaussian","plane_mc","gaussian_mc"};
		const std::string LSTYPES[] = { "plane","gaussian","ring","tophat","line","point","plane_mc","gaussian_mc","ring_mc","line_mc","point_mc" };
	
	   /**
		* @brief check if the fname has a given extension
		* This function checks if the fname ends with a given extension
		* \param[in] fname Filename
		* \param[in] extension Extension to search (without "."!)
		* \return True if fname has the given extension
	    */
        bool findExtension (std::string fname, std::string extension); 
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
			 /**
			  * @brief read the XML-file
			  * This methode reads the XML-file 
			  * @param fname the filename
			  * @param path path of the XML-file. If path is empty, the path will be extracted from the filename. With this, the library is able to 
			  * handle relative paths for the files used e.g. for the detectors
			  */
				void readXML(std::string fname, bool calc_enabled=true, std::string path = "");
				GOAT::raytracing::Scene S; ///< The scene that was read from the file is saved here				
				void setEnableCalculation(bool enable) { calculation_enabled=enable;}
				bool isCalculationEnabled() {return calculation_enabled;}
			private:				
                void readScene(); ///< read the entire Scene
				void readLightSources(); ///< (used in readScene) read the light sources from the file
				void readCommands(); ///< (used in readXML) read and execute the commands for calculation
				void readObjects(); ///< (used in readScene) read the objects from the file
				void readDetectors(); ///< (used in readScene) read the detectors from the file (deprecated ?)
				void doCalculations(); ///< (used in readScene) read and execute the commands for calculation
                void doPulseCalculation(tinyxml2::XMLElement* objEll); ///< (used in doCalculations) to pulsed Calculation (rt + integral)
				void doPulseCalculation_rt(tinyxml2::XMLElement* objEll); ///< do pulse Calculation with raytracing only
				/**
				 * @brief read double vector from file
				 * This method reads a 3D-Vector from the XML-file represented by the corresponding XMLElement 
				 * @param ell XML-Element to read from
				 * @param x default parameter for x-component (if not given in the file)
				 * @param y default parameter for y-component (if not given in the file)
				 * @param z default parameter for z-component (if not given in the file)
				 * @return double 3D vector
				 */
				GOAT::maths::Vector<double> readVector(tinyxml2::XMLElement* ell, double x = 0, double y = 0, double z = 0); 

				/**
				 * @brief read double vector from file
				 * This method reads a 3D-Vector from the XML-file represented by the corresponding XMLElement. The default value for 
				 * not given attributes is 0 for all components
				 * @param[in] ell XML-Element to read from
				 * @param[out] xmlError error number from tinyxml
				 */
                GOAT::maths::Vector<double> readVector(tinyxml2::XMLElement* ell, int& xmlError);

				/**
				 * @brief read double vector from file
				 * This method reads a complex valued 3D-Vector from the XML-file represented by the corresponding XMLElement. The default value for 
				 * not given attributes is 0 for all components
				 * @param[in] ell XML-Element to read from
				 * @param[out] xmlError error number from tinyxml
				 * @return complex 3D vector
				 */
                GOAT::maths::Vector<std::complex<double> > readCmplxVector (tinyxml2::XMLElement* ell, int& xmlError);
				
				/**
				 * @brief read complex value from file
				 * This method reads a complex value from file with given default values for the real and the imaginary part
				 * @param ell XML-Element to read from
				 * @param defre default value for the real part
				 * @param defim default value for the imaginary part 
				 * @return complex value
				 */
				std::complex<double> readCmplx(tinyxml2::XMLElement* ell, double defre=0.0, double defim=0.0 );

				/**
				 * @brief read complex value from file
				 * This method reads a complex value from file with given default values for the real and the imaginary part
				 * @param[in] ell XML-Element to read from
				 * @param[out] xmlError error number from tinyxml
				 * @return complex value
				 */
                std::complex<double> readCmplx(tinyxml2::XMLElement* ell, int& xmlError);

				tinyxml2::XMLNode* rootElement; ///< pointer to root element of the XML
				tinyxml2::XMLElement* sceneElement;	 ///< pointer to the scene element of the XML			 
				tinyxml2::XMLElement* calculationElement; ///< pointer to the calculation elemeent of the XML
				std::vector<GOAT::raytracing::ObjectShape*> Obj; ///< vector, which carries all objects (as pointers)
				std::vector<GOAT::raytracing::Detector*> Det; ///< vector, which carries all detectors (as pointers)
				std::vector< GOAT::raytracing::LightSrc*> LS; ///< vector, which carries all light sources (as pointers)
				int numObj = 0; ///< number of objects
				int numLS = 0; ///< number of light sources
				int numDet = 0; ///< number of detectors
				std::string path; ///< path of the XML-File
				bool calculation_enabled=true; 
		};




        class xmlWriter
        {
			public:
                xmlWriter(const GOAT::raytracing::Scene &S);
                void write (std::string fname);
			private:
			inline std::string formatDouble(double val, int precision = 17) 
				{
 					   std::ostringstream oss;
    					oss.imbue(std::locale::classic()); // Erzwingt "." statt "," als Dezimaltrennzeichen
    					oss << std::fixed << std::setprecision(precision) << val;
    					return oss.str();
				}
				void writeLightSrc(int i); ///< write the i-th light source to the file
				void writeObject(int i); ///< write the i-th object to the file
				void writeDetector(int i); ///< write the i-th detector to the file
				

				/**
				 * @brief write double vector to file
				 * This method writes a double 3D vector to XML-file
				 * @param name Element name of the vector
				 * @param v double vector
				 * @return pointer to the corresponding XMLElemment
				 */
				tinyxml2::XMLElement* writeVectorD(std::string name, maths::Vector<double> v);
				
				/**
				 * @brief write complex vector to file
				 * This method writes a complex 3D vector to XML-file
				 * @param name Element name of the vector
				 * @param v double vector
				 * @return pointer to the corresponding XMLElemment
				 */
				tinyxml2::XMLElement* writeVectorC(std::string name, maths::Vector<std::complex<double>> v);

				/**
				 * @brief write complex number to file
				 * This method writes a complex number to the XML-file
				 * @param  name Element name of the complex number
				 * @param z complex number 
				 */
				tinyxml2::XMLElement* writeComplex(std::string name, std::complex<double> z);
                const GOAT::raytracing::Scene &S; ///< the scene
                tinyxml2::XMLDocument doc; ///< the xml document
				tinyxml2::XMLElement* root; ///< root XML Element
				tinyxml2::XMLElement* scene; ///< XML Element to the Scene section
				tinyxml2::XMLElement* lightSrcs; ///< XML Element to the LightSources section
				tinyxml2::XMLElement* objects; ///< XML Element to the Objects section
				tinyxml2::XMLElement* detectors; ///< XML Element to the Detectors section 
				
        };
	}
}
