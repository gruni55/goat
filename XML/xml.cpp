
#include "tinyxml2.h"
#include "xml.h"


namespace GOAT
{
	namespace XML
	{
		void xmlReader::readXML(std::string fname, GOAT::raytracing::Scene& S)
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
						tinyxml2::XMLElement* lsEll = ell->FirstChildElement("LightSource");
						GOAT::maths::Vector<double> Pos;
						int numRays;
						double wavelength;
						double size;
						std::string typeStr;

						while (lsEll)
						{
							Pos = readVector(lsEll->FirstChildElement("Position"));
							numRays = lsEll->IntAttribute("NumRays", 100);
							wavelength = lsEll->DoubleAttribute("Wavelength", 1.0);
							size = lsEll->DoubleAttribute("Size", 10.0);
							typeStr = lsEll->Attribute("Type", "plane");
							if (typeStr.compare("plane") == 0)
							{
								GOAT::raytracing::LightSrcPlane* ls = new GOAT::raytracing::LightSrcPlane(Pos, numRays, size, wavelength);
								GOAT::maths::Vector<double> k=readVector(lsEll->FirstChildElement("Direction"));
								ls->setk(k);
								S.addLightSource(ls);
							}

							if (typeStr.compare("plane_mc") == 0)
							{
								GOAT::raytracing::LightSrcPlane_mc* ls = new GOAT::raytracing::LightSrcPlane_mc(Pos, numRays, size, wavelength);
								GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
								ls->setk(k);
								S.addLightSource(ls);
							}


							if (typeStr.compare("gaussian") == 0)
							{
								double w0;
								double NA;
								tinyxml2::XMLError err;
								w0 = lsEll->DoubleAttribute("w0", 1.0);
								GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
								GOAT::raytracing::LightSrcGauss* ls = new GOAT::raytracing::LightSrcGauss(Pos, numRays, wavelength,w0,focusPos);
								err = lsEll->QueryDoubleAttribute("NA", &NA);
								if (err == tinyxml2::XML_SUCCESS) ls->setNA(NA);
							}

							if (typeStr.compare("gaussian_mc") == 0)
							{
								double w0;
								double NA;
								tinyxml2::XMLError err;
								w0 = lsEll->DoubleAttribute("w0", 1.0);
								GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
								GOAT::raytracing::LightSrcGauss_mc* ls = new GOAT::raytracing::LightSrcGauss_mc(Pos, numRays, wavelength, w0, focusPos);
								err = lsEll->QueryDoubleAttribute("NA", &NA);
								if (err == tinyxml2::XML_SUCCESS) ls->setNA(NA);
							}


						};
					}
				}
			}
		}

		GOAT::maths::Vector<double> xmlReader::readVector(tinyxml2::XMLElement* ell)
		{
			GOAT::maths::Vector<double> P;
			P[0] = ell->DoubleAttribute("x", 0.0);
			P[1] = ell->DoubleAttribute("y", 0.0);
			P[2] = ell->DoubleAttribute("z", 0.0);
			return P;
		}

		GOAT::maths::Vector<double> xmlReader::readVector(tinyxml2::XMLElement* ell, int &xmlError)
		{
			double x, y, z;
			xmlError = ell->QueryDoubleAttribute("x",&x);
			xmlError = xmlError | ell->QueryDoubleAttribute("y", &y);
			xmlError = xmlError | ell->QueryDoubleAttribute("z", &z);			
			return GOAT::maths::Vector<double>(x,y,z);
		}
	}
}