
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

						while (lsEll)
						{

                          
						};
					}
				}
			}
		}
	}
}