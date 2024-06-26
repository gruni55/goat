
#include "tinyxml2.h"
#include "xml.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation.h"

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
					
					int type = mapString2LightSourceToken(typeStr);
					switch (type)
					{
					case TOKEN_LIGHTSOURCE_PLANE: {
													LS[numLS] = new GOAT::raytracing::LightSrcPlane(Pos, numRays, wavelength, size);
													GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
													LS[numLS]->setk(k);
												   }
												   break;

					case TOKEN_LIGHTSOURCE_PLANE_MC: {
														LS[numLS] = new GOAT::raytracing::LightSrcPlane_mc(Pos, numRays, wavelength, size);
														GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
														LS[numLS]->setk(k);
													 }
														  break;

					case TOKEN_LIGHTSOURCE_GAUSSIAN: {
														double w0;
														double NA;
														tinyxml2::XMLError err;
														w0 = lsEll->DoubleAttribute("w0", 1.0);
														GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
														LS[numLS] = new GOAT::raytracing::LightSrcGauss(Pos, numRays, wavelength, w0, focusPos);
														err = lsEll->QueryDoubleAttribute("NA", &NA);
														if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss*)ls)->setNA(NA);
													  }
													 	  break;

					case TOKEN_LIGHTSOURCE_GAUSSIAN_MC: {
															double w0;
															double NA = 0.9;
															tinyxml2::XMLError err;
															w0 = lsEll->DoubleAttribute("w0", 1.0);
															GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
															LS[numLS] = new GOAT::raytracing::LightSrcGauss_mc(Pos, numRays, wavelength, w0, focusPos);
															err = lsEll->QueryDoubleAttribute("NA", &NA);
															if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss_mc*)ls)->setNA(NA);
														}
													  break;
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
			// GOAT::raytracing::surface objS;
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
					
					Pos = readVector(objEll->FirstChildElement("Position"));
					typeStr = objEll->Attribute("Type");
					alpha = objEll->DoubleAttribute("Alpha", 0.0);
					beta = objEll->DoubleAttribute("Beta", 0.0);
					gamma = objEll->DoubleAttribute("Gamma", 0.0);
					isActive = objEll->BoolAttribute("IsActive", false);
					// GOAT::raytracing::ObjectShape* obj = NULL;
					n = readCmplx(objEll->FirstChildElement("n"), 1.0);
					int type = mapString2ObjectToken(typeStr);
					switch (type)
					{
					case TOKEN_OBJECT_ELLIPSOID: {
													GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10.0, 10.0, 10.0);
													Obj.push_back(new GOAT::raytracing::Ellipsoid(Pos, Dimensions, n));
													Obj[numObj]->setMatrix(alpha, beta, gamma);
													Obj[numObj]->setActive(isActive);
													S.addObject(Obj[numObj]);
												 }
											   break;
					case TOKEN_OBJECT_BOX: {
													GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10, 10, 10);
													Obj.push_back(new GOAT::raytracing::Box(Pos, Dimensions, n));
													Obj[numObj]->setMatrix(alpha, beta, gamma);
													Obj[numObj]->setActive(isActive);
													S.addObject(Obj[numObj]);
										    }
										 break;
					case TOKEN_OBJECT_SURFACE: 
										   {
											Obj.push_back(new GOAT::raytracing::surface(Pos, n));
											fileTypeStr = objEll->Attribute("Filetype");

											if (fileTypeStr.compare("srf") == 0)
											{
												fileName = objEll->Attribute("Filename");
												if (!fileName.empty()) ((GOAT::raytracing::surface*)Obj[numObj])->createsurface(fileName);
											}

											if (fileTypeStr.compare("stl") == 0)
											{
												fileName = objEll->Attribute("Filename");
												if (!fileName.empty())
												((GOAT::raytracing::surface*)Obj[numObj])->importBinSTL(fileName);
											}
										   }
										   Obj[numObj]->setMatrix(alpha, beta, gamma);
										   Obj[numObj]->setActive(isActive);
										   S.addObject(Obj[numObj]);
										   break;

					case TOKEN_OBJECT_SPHERIC_LENS:
											{
												GOAT::raytracing::lensParms lensparms;												
												std::string leftCurvatureStr, rightCurvatureStr;
												
												// left Side
												tinyxml2::XMLElement* leftEll = objEll->FirstChildElement("left");
												leftCurvatureStr = leftEll->Attribute("Curvature");
												if (!leftCurvatureStr.empty())
												{
													if (leftCurvatureStr.compare("flat") == 0) lensparms.left.curvature = GOAT::raytracing::flat;
													if (leftCurvatureStr.compare("concave") == 0) lensparms.left.curvature = GOAT::raytracing::concave;
													if (leftCurvatureStr.compare("convex") == 0) lensparms.left.curvature = GOAT::raytracing::convex;
												}
												else lensparms.left.curvature = GOAT::raytracing::flat;
												lensparms.left.R=leftEll->DoubleAttribute("R",0.0);
												
												tinyxml2::XMLElement* rightEll = objEll->FirstChildElement("right");
												rightCurvatureStr = rightEll->Attribute("Curvature");
												if (!rightCurvatureStr.empty())
												{
												if (rightCurvatureStr.compare("flat") == 0) lensparms.right.curvature = GOAT::raytracing::flat;
												if (rightCurvatureStr.compare("concave") == 0) lensparms.right.curvature = GOAT::raytracing::concave;
												if (rightCurvatureStr.compare("convex") == 0) lensparms.right.curvature = GOAT::raytracing::convex;
												}
												else lensparms.right.curvature = GOAT::raytracing::flat;
												lensparms.right.R = rightEll->DoubleAttribute("R", 0.0);
												lensparms.offset = objEll->DoubleAttribute("Offset", 0.0);
												lensparms.radius = objEll->DoubleAttribute("Radius", 0.0);



												Obj.push_back(new GOAT::raytracing::sphericLens(Pos,n,lensparms));
												Obj[numObj]->setMatrix(alpha, beta, gamma);
												Obj[numObj]->setActive(isActive);
												S.addObject(Obj[numObj]);							
											}
					}
									
					numObj++;
				} // while loop
				doCalculations();
				  // S.addObjectList(numObj, Obj);

			}
			/* End of objects */
		}

		void xmlReader::doCalculations()
		{
			int type;
			std::string typeStr;
			tinyxml2::XMLElement* ell;
			ell = rootElement->FirstChildElement("Calculations");
			if (ell != NULL)
			{
				for (tinyxml2::XMLElement* objEll = ell->FirstChildElement("Calculation"); objEll != NULL; objEll = objEll->NextSiblingElement("Calculation"))
				{
					std::string inactiveStr;
					const char* hStr;
						hStr=objEll->Attribute("Inactive");
						if (hStr != NULL) inactiveStr = hStr;
						else inactiveStr = "false";
						if (inactiveStr.compare("true")!=0)
					{
						
						typeStr = objEll->Attribute("Type");
						if (!typeStr.empty())
						{
							type = mapString2CalculationToken(typeStr);
							switch (type)
							{
							case TOKEN_CALCULATION_PATH:
							{
								std::string fname = objEll->Attribute("Filename");
								int numRays;
								int numReflex;
								if (!fname.empty())
								{
									GOAT::raytracing::Raytrace_Path rt(S);
									// optionally, the number of rays of all sources can be set to numRays
									numRays = objEll->IntAttribute("numRays", 0);									
									std::vector<int> numRays_old;
									if (numRays > 0)
									{
										// store the old values 									 
										for (int i = 0; i < S.nLS; i++)
										{
											numRays_old.push_back(S.LS[i]->getNumRays());
											S.LS[i]->setNumRays(numRays);
										}
									}

									numReflex = objEll->IntAttribute("numReflex", 0);									
									rt.setNumReflex(numReflex);

									rt.trace(fname);

									if (numRays > 0)
									{
										// restore the old values  
										for (int i = 0; i < S.nLS; i++)
											S.LS[i]->setNumRays(numRays_old[i]);
									} 
								}
								else
									std::cerr << "Path calculation: You forgot to give an appropriate file name for the output!!" << std::endl;
								break;
							} // case path calculation
							case TOKEN_CALCULATION_PULSE:
							{
								std::string fname = objEll->Attribute("Filename");
								if (!fname.empty())
								{
									GOAT::raytracing::pulseCalculation pc(S);
									GOAT::raytracing::TrafoParms trafoparms;
									trafoparms = pc.getTrafoParms();
									pc.setCenterWavelength(objEll->DoubleAttribute("Wavelength", trafoparms.wvl));
									pc.setNumReflex(objEll->IntAttribute("NumReflexions", trafoparms.nR));
									pc.setNumWavelengthsPerRange(objEll->IntAttribute("NumWavelengthsPerRange", trafoparms.nS));
									pc.setSpectralRanges(objEll->IntAttribute("NumSpectralRanges", trafoparms.nI));
									pc.setReferenceTime(objEll->IntAttribute("Reference_time", pc.getReferenceTime()));
									double repRate = objEll->DoubleAttribute("Repetition_rate", -1);
									if (repRate > 0) pc.setRepetitionRate(repRate);
									double dx = 2.0 * S.r0 / (double)pc.getNumCellsPerDirection();
									pc.setSpatialResolution(objEll->DoubleAttribute("Spatial_resolution", dx));
								

									char cs[3];
									std::vector< std::function< std::complex< double >(double) > > nList;
									std::string refFuncName;
									bool failed = false;

									tinyxml2::XMLElement* refEll = objEll->FirstChildElement("RefractiveIndexList");
									if (refEll == NULL)
									{
										std::cerr << "Pulse calculation: Refractive index function list is missing! Stopped calculation" << std::endl;
										break;
									}
								
									std::string refStr;
									int refIndexToken;
									for (int i=0; (i<S.nObj) && (!failed); i++)
									{
										sprintf(cs, "n%i", i);
										hStr = refEll->Attribute(cs);
										if (hStr == NULL)
										{
											
											std::cerr << "Pulse calculation: Refractive index function for object " << i << " is missing! Stopped calculation" << std::endl;
											failed = true;											
										}
										else
										{
											refStr = hStr;
											refIndexToken = mapString2RefractiveIndexToken(refStr);
											if (refIndexToken == TOKEN_NOT_FOUND)
											{
												std::cerr << "Pulse calculation: Wrong refractive index function name (" << refStr << ") !Calculation stopped!" << std::endl;
												failed = true;
											}
											addFunction2IndexList(nList, refIndexToken);
										}
									}

									if (failed) break;


									hStr = refEll->Attribute("nS");
									if (hStr == NULL)
									{
										std::cerr << "Pulse calculation: Refractive index function for surrounding medium is missing! Stopped calculation" << std::endl;
										failed = true;
									}
									if (failed) break;
									refStr = hStr;
									refIndexToken = mapString2RefractiveIndexToken(refStr);
									if (refIndexToken == TOKEN_NOT_FOUND)
									{
										std::cerr << "Pulse calculation: Wrong refractive index function name (" << refStr << ") !Calculation stopped!" << std::endl;
										failed = true;
									}

									if (failed) break;
									addFunction2IndexList(nList, refIndexToken);
									pc.setRefractiveIndexFunctions(nList);

									double time = objEll->DoubleAttribute("Time", -1);
									if (time < 0)
									{
										double offset = objEll->DoubleAttribute("TimeOffset", 0);
										int objEstimate = objEll->IntAttribute("EstimateTimeForObject", 0);
										time = pc.findHitTime(objEstimate) + offset;
									}

									
									pc.field(time);

									std::string fullfname;

									for (int i = 0; i < S.nObj; i++)
									{
										if (S.Obj[i]->isActive())
										{
											fullfname = fname + std::to_string(i) + ".dat";
											GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
										}
									}
								}
								else
									std::cerr << "Path calculation: You forgot to give an appropriate file name for the output!!" << std::endl;
								break;
							}
							} // switch 
						} // is type given ?
					} // is inactive ?
				} // for loop
			} // if no Calculations
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
