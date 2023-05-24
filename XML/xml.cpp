
#include "tinyxml2.h"
#include "xml.h"


namespace GOAT
{
	namespace XML
	{
		void xmlReader::readXML(std::string fname)
		{
			tinyxml2::XMLElement* ell;
			tinyxml2::XMLDocument doc;
			double dv;
			doc.LoadFile(fname.c_str());
			rootElement = doc.RootElement();
			if (rootElement != NULL)
			{
				/*  Looking for the scene definition */
				sceneElement= rootElement->FirstChildElement("Scene");
				if (sceneElement != NULL)
				{
					dv = sceneElement->DoubleAttribute("r0",1000.0);    // reading the radius of the calculation sphere (if not: default 1000)
					S.setr0(dv);                                     

					ell = sceneElement->FirstChildElement("nS");        // reading refractive index of the surrounding medium 
					if (ell != NULL)
					{
						
						double re = ell->DoubleAttribute("real",1.0);
						double im = ell->DoubleAttribute("imag", 0.0);
						S.setnS(std::complex<double>(re, im));
					}

					/* Looking for light sources */
					ell = sceneElement->FirstChildElement("LightSources");
					if (ell != NULL)
					{
						// tinyxml2::XMLElement* lsEll = ell->FirstChildElement("LightSource");
						GOAT::maths::Vector<double> Pos;
						int numRays;
						double wavelength;
						double size;
						GOAT::raytracing::LightSrc* ls = NULL;

						for (tinyxml2::XMLElement *lsEll= ell->FirstChildElement("LightSource"); lsEll!=NULL; lsEll=lsEll->NextSiblingElement("LightSource"))
						{	
							std::string typeStr;
							typeStr = lsEll->Attribute("Type");							
							Pos = readVector(lsEll->FirstChildElement("Position"));
							numRays = lsEll->IntAttribute("NumRays", 100);
							wavelength = lsEll->DoubleAttribute("Wavelength", 1.0);
							size = lsEll->DoubleAttribute("Size", 10.0);
							
							if (typeStr.compare("plane") == 0)
							{
								ls = new GOAT::raytracing::LightSrcPlane(Pos, numRays, size, wavelength);
								GOAT::maths::Vector<double> k=readVector(lsEll->FirstChildElement("Direction"));
								ls->setk(k);								
							}

							if (typeStr.compare("plane_mc") == 0)
							{
								ls = new GOAT::raytracing::LightSrcPlane_mc(Pos, numRays, size, wavelength);
								GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
								ls->setk(k);								
							}


							if (typeStr.compare("gaussian") == 0)
							{
								double w0;
								double NA;
								tinyxml2::XMLError err;
								w0 = lsEll->DoubleAttribute("w0", 1.0);
								GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
								ls = new GOAT::raytracing::LightSrcGauss(Pos, numRays, wavelength,w0,focusPos);
								err = lsEll->QueryDoubleAttribute("NA", &NA);
								if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss *)ls)->setNA(NA);
							}

							if (typeStr.compare("gaussian_mc") == 0)
							{
								double w0;
								double NA=0.9;
								tinyxml2::XMLError err;
								w0 = lsEll->DoubleAttribute("w0", 1.0);
								GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
								ls = new GOAT::raytracing::LightSrcGauss_mc(Pos, numRays, wavelength, w0, focusPos);
								err = lsEll->QueryDoubleAttribute("NA" ,&NA);
								if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss_mc*)ls)->setNA(NA);
							}
							
							if (ls != NULL)
							{
								S.addLightSource(ls);
							}

						} // while loop						
					}
					/* End of light sources */

					/* Now, let's look for objects */
					ell = sceneElement->FirstChildElement("Objects");
					if (ell != NULL)
					{
						GOAT::maths::Vector<double> Pos;
						std::string typeStr;
						std::complex<double> n;
						bool isActive;
						double alpha = 0;
						double beta = 0;
						double gamma = 0;

						for (tinyxml2::XMLElement* objEll = ell->FirstChildElement("Object"); objEll != NULL; objEll = objEll->NextSiblingElement("Object"))
						{							
							Pos = readVector(objEll->FirstChildElement("Position"));
							typeStr = objEll->Attribute("Type");
							alpha = objEll->DoubleAttribute("Alpha", 0.0);
							beta = objEll->DoubleAttribute("Beta", 0.0);
							gamma = objEll->DoubleAttribute("Gamma", 0.0);
							isActive = objEll->BoolAttribute("IsActive", false);
							GOAT::raytracing::ObjectShape* obj = NULL;
							n = readCmplx(objEll->FirstChildElement("n"), 1.0);
							if (typeStr.compare("ellipsoid") == 0) // is an elliptic particle
							{
								GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10, 10, 10);
								
								obj = new GOAT::raytracing::Ellipsoid(Pos, Dimensions,n);								
							}

							if (typeStr.compare("box") == 0) // is a box
							{
								GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10, 10, 10);
								GOAT::raytracing::Box *obj= new GOAT::raytracing::Box(Pos, Dimensions, n);								
							}

							if (obj != NULL)
							{
								obj->setMatrix(alpha, beta, gamma);
								obj->setActive(isActive);
								S.addObject(obj);
							}

						} // while loop
					}
					/* End of objects */
				}
			}
		}

		GOAT::maths::Vector<double> xmlReader::readVector(tinyxml2::XMLElement* ell,  double x, double y, double z)
		{
			GOAT::maths::Vector<double> P;
			if (ell != NULL)
			{			
				P[0] = ell->DoubleAttribute("x", x);
				P[1] = ell->DoubleAttribute("y", y);
				P[2] = ell->DoubleAttribute("z", z);
			}
			return P;
		}

		GOAT::maths::Vector<double> xmlReader::readVector(tinyxml2::XMLElement* ell, int &xmlError)
		{
			double x=0, y=0, z=0;
			if (ell != NULL)
			{
				xmlError = ell->QueryDoubleAttribute("x", &x);
				xmlError = xmlError | ell->QueryDoubleAttribute("y", &y);
				xmlError = xmlError | ell->QueryDoubleAttribute("z", &z);
			}
			return GOAT::maths::Vector<double>(x,y,z);
		}

		std::complex<double> xmlReader::readCmplx(tinyxml2::XMLElement* ell, double defre, double defim)
		{
			double im=defim, re=defre;
			if (ell != NULL)
			{
				re = ell->DoubleAttribute("real", defre);
				im = ell->DoubleAttribute("imag", defim);
			}
			return std::complex<double>(re, im);
		}
	}
}