
#include "tinyxml2.h"
#include "xml.h"

namespace GOAT
{
	namespace XML
	{
		void xmlReader::readXML(const char* fname)
		{
			
			tinyxml2::XMLDocument doc;
			tinyxml2::XMLError eResult=doc.LoadFile(fname);
						
			if (eResult == tinyxml2::XML_SUCCESS)
			{
				rootElement = doc.RootElement();
				readScene();
			}
			else
				std::cerr << "Could not read XML-File:" << fname << std::endl;
			
		}

		void xmlReader::readScene()
		{
			sceneElement = rootElement->FirstChildElement("Scene");
			if (sceneElement != NULL)
			{
				double dv;

				dv = sceneElement->DoubleAttribute("r0", 1000.0);    // reading the radius of the calculation sphere (if not: default 1000)
				S.setr0(dv);
				tinyxml2::XMLElement* ell;
				ell = sceneElement->FirstChildElement("nS");        // reading refractive index of the surrounding medium 
				if (ell != NULL)
				{

					double re = ell->DoubleAttribute("real", 1.0);
					double im = ell->DoubleAttribute("imag", 0.0);
					S.setnS(std::complex<double>(re, im));
				}

				/* Looking for light sources */
			    readLightSources();
				

				/* Now, let's look for objects */
				readObjects();

			}		
		}
		void xmlReader::readLightSources()
		{
			tinyxml2::XMLElement* ell;
			ell = sceneElement->FirstChildElement("LightSources");
			if (ell != NULL)
			{
			   GOAT::maths::Vector<double> Pos;
				int numRays;
				double wavelength;
				double size;
				GOAT::raytracing::LightSrc* ls = NULL;
				for (tinyxml2::XMLElement* lsEll = ell->FirstChildElement("LightSource"); lsEll != NULL; lsEll = lsEll->NextSiblingElement("LightSource"))
				{
					std::string typeStr;
					typeStr = lsEll->Attribute("Type");
					Pos = readVector(lsEll->FirstChildElement("Position"));
					numRays = lsEll->IntAttribute("NumRays", 100);
					wavelength = lsEll->DoubleAttribute("Wavelength", 1.0);
					size = lsEll->DoubleAttribute("Size", 10.0);
					LS = (GOAT::raytracing::LightSrc**)realloc(LS, numLS + 1);

					if (typeStr.compare("plane") == 0)
					{
						LS[numLS] = new GOAT::raytracing::LightSrcPlane(Pos, numRays, wavelength, size);
						GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
						LS[numLS]->setk(k);
					}

					if (typeStr.compare("plane_mc") == 0)
					{
						LS[numLS] = new GOAT::raytracing::LightSrcPlane_mc(Pos, numRays, wavelength, size);
						GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
						LS[numLS]->setk(k);
					}


					if (typeStr.compare("gaussian") == 0)
					{
						double w0;
						double NA;
						tinyxml2::XMLError err;
						w0 = lsEll->DoubleAttribute("w0", 1.0);
						GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
						LS[numLS] = new GOAT::raytracing::LightSrcGauss(Pos, numRays, wavelength, w0, focusPos);
						err = lsEll->QueryDoubleAttribute("NA", &NA);
						if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss*)ls)->setNA(NA);
					}

					if (typeStr.compare("gaussian_mc") == 0)
					{
						double w0;
						double NA = 0.9;
						tinyxml2::XMLError err;
						w0 = lsEll->DoubleAttribute("w0", 1.0);
						GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
						LS[numLS] = new GOAT::raytracing::LightSrcGauss_mc(Pos, numRays, wavelength, w0, focusPos);
						err = lsEll->QueryDoubleAttribute("NA", &NA);
						if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss_mc*)ls)->setNA(NA);
					}

					numLS++;

					
				} // while loop			
				S.addLightSourceList(numLS, LS);
			}
		}

		void xmlReader::readCommands()
		{
			int numReflex;
			std::string typeStr, fname;
			tinyxml2::XMLElement* ell;
			ell = rootElement->FirstChildElement("Calculations");
			if (ell != NULL)
			{
				for (tinyxml2::XMLElement* calcEll = ell->FirstChildElement("Calculation"); calcEll != NULL; calcEll = calcEll->NextSiblingElement("Calculation"))
				{
					typeStr = calcEll->Attribute("Type");
					fname = calcEll->Attribute("Filename");
					numReflex = calcEll->IntAttribute("numReflex", 0);
					if (typeStr.compare("path") && (!fname.empty()) )
					{
						bool outgoingRays= calcEll->BoolAttribute("outgoing rays", false);
						GOAT::raytracing::Raytrace_Path rp(S);
						rp.setShowOutgoingRays(outgoingRays);
						rp.setNumReflex(numReflex);
						rp.trace(fname);
					}
					

				}
			}

		}

		void xmlReader::readObjects()
		{
			tinyxml2::XMLElement* ell;
			GOAT::raytracing::surface objS;
			GOAT::maths::Vector<double> Pos;
			std::string typeStr;
			std::string fileTypeStr;
			std::string fileName;
			std::complex<double> n;
			bool isActive;
			double alpha = 0;
			double beta = 0;
			double gamma = 0;

			ell = sceneElement->FirstChildElement("Objects");
			if (ell != NULL)
			{
				

				for (tinyxml2::XMLElement* objEll = ell->FirstChildElement("Object"); objEll != NULL; objEll = objEll->NextSiblingElement("Object"))
				{		
					/*if (numObj == 0)
						Obj = (GOAT::raytracing::ObjectShape**) malloc(sizeof(GOAT::raytracing::ObjectShape*));
					else 
						Obj = (GOAT::raytracing::ObjectShape**) realloc(Obj, numObj + 1);*/
					Pos = readVector(objEll->FirstChildElement("Position"));
					typeStr = objEll->Attribute("Type");
					alpha = objEll->DoubleAttribute("Alpha", 0.0);
					beta = objEll->DoubleAttribute("Beta", 0.0);
					gamma = objEll->DoubleAttribute("Gamma", 0.0);
					isActive = objEll->BoolAttribute("IsActive", false);
					// GOAT::raytracing::ObjectShape* obj = NULL;
					n = readCmplx(objEll->FirstChildElement("n"), 1.0);
					if (typeStr.compare("ellipsoid") == 0) // is an elliptic particle
					{
						GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10.0, 10.0, 10.0);
						// Obj[numObj] = new GOAT::raytracing::Ellipsoid(Pos, Dimensions, n);
						Obj.push_back(new GOAT::raytracing::Ellipsoid(Pos, Dimensions, n));
						S.addObject(Obj[numObj]);
					}

					if (typeStr.compare("box") == 0) // is a box
					{
						GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10, 10, 10);
						Obj[numObj] = new GOAT::raytracing::Box(Pos, Dimensions, n);
					}


					//GOAT::raytracing::ObjectShape* obj;
				/*	if (typeStr.compare("surface") == 0)
					{
						Obj[numObj] = new GOAT::raytracing::surface(Pos, n);
						fileTypeStr = objEll->Attribute("Filetype");

						if (fileTypeStr.compare("srf") == 0)
						{
							fileName = objEll->Attribute("Filename");
							if (!fileName.empty()) objS.createsurface(fileName);
						}

						if (fileTypeStr.compare("stl") == 0)
						{
							fileName = objEll->Attribute("Filename");
							if (!fileName.empty())
								((GOAT::raytracing::surface *)Obj[numObj])->importBinSTL(fileName);
						}

					}*/
					numObj++;
				} // while loop
				  // S.addObjectList(numObj, Obj);

			}
			/* End of objects */
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
			else
				P = GOAT::maths::Vector<double>(x, y, z);
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