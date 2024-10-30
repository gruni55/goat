

#include "tinyxml2.h"
#include "xml.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation.h"
#include "pulsecalculation_rt.h"
#include "pulsecalculation_field.h"
#include "raytrace_inel.h"
#include <chrono>

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

				/* look for the detectors */
				readDetectors();

				/* Looking for light sources */
			    readLightSources();
				

                /* Now, let's look for objects */
				readObjects();

                doCalculations();
			}		
		}

        void xmlReader::readDetectors()
		{
			tinyxml2::XMLElement* ell;
			ell = sceneElement->FirstChildElement("Detectors");

			if (ell != NULL)
			{
				int n1, n2;
				
				for (tinyxml2::XMLElement* detEll = ell->FirstChildElement("Detector"); detEll != NULL; detEll = detEll->NextSiblingElement("Detector"))
				{
					maths::Vector<double> Pos = readVector(detEll->FirstChildElement("Position"));
					maths::Vector<double> Dir = readVector(detEll->FirstChildElement("Direction"));
					std::string typeStr;
					typeStr = detEll->Attribute("Type");                    
					std::string filename;
					filename = detEll->Attribute("Filename");
					int type = mapString2DetectorToken(typeStr);
					switch (type)
					{
					  case TOKEN_DETECTOR_PLANE : 
													{
                                                     double d = detEll->DoubleAttribute("d", 1);
													 int n = detEll->IntAttribute("n", 1);
                                                     Det.push_back(new raytracing::DetectorPlane(Pos, Dir, d,n));
													 Det[numDet]->fname = filename;
													 S.addDetector(Det[numDet]);
                                                    numDet++;
													}
													
					}
				}
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
                    case TOKEN_LIGHTSOURCE_LINE: {
                                                        
                                                        GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
                                                        GOAT::maths::Vector<double> D = readVector(lsEll->FirstChildElement("lateral_direction"));
                                                        LS[numLS] = new GOAT::raytracing::LightSrcLine(Pos, numRays, wavelength, size, k, D);
                                                        }
                                               break;
                    case TOKEN_LIGHTSOURCE_LINE_MC: {
                                                        GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
                                                        GOAT::maths::Vector<double> D = readVector(lsEll->FirstChildElement("lateral_direction"));
                                                        LS[numLS] = new GOAT::raytracing::LightSrcLine_mc(Pos, numRays, wavelength, size, k, D);
                                                    }
                                               break;

					case TOKEN_LIGHTSOURCE_GAUSSIAN: {
														double w0;
														double NA=1.0;
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
															double NA = 1.0;
															tinyxml2::XMLError err;
															w0 = lsEll->DoubleAttribute("w0", 1.0);
															GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
															LS[numLS] = new GOAT::raytracing::LightSrcGauss_mc(Pos, numRays, wavelength, w0, focusPos);
															err = lsEll->QueryDoubleAttribute("NA", &NA);
															if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss_mc*)ls)->setNA(NA);
														}
													  break;
                    case TOKEN_LIGHTSOURCE_RING:        {
                                                         double rmin, rmax;
                                                         rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                         rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                         LS[numLS]=new GOAT::raytracing::LightSrcRing(Pos, numRays, wavelength, rmin,rmax);
														 GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
														 LS[numLS]->setk(k);
                                                        break;
                                                        }
                    case TOKEN_LIGHTSOURCE_RING_MC:
                                                        {
                                                          double rmin, rmax;
                                                          rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                          rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                           LS[numLS]=new GOAT::raytracing::LightSrcRing_mc(Pos, numRays, wavelength, rmin,rmax);
														   GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
														   LS[numLS]->setk(k);
														   break;
                                                        }
                    case TOKEN_LIGHTSOURCE_GAUSSIAN_RING_MC :
                                                        {
                                                         double rmin, rmax;
                                                         double width;
                                                         rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                         rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                         width=lsEll->DoubleAttribute("width",rmax);
                                                          LS[numLS]=new GOAT::raytracing::LightSrcRingGauss_mc(Pos, numRays, wavelength, rmin, rmax);
                                                         ((GOAT::raytracing::LightSrcRingGauss_mc *)LS[numLS])->setFWHM(width);
														 GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
														 LS[numLS]->setk(k);
                                                         break;
                                                        }


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
                                                break;						
											}
                    case TOKEN_OBJECT_CONE: 
                                            {
                                                double height, radius;
                                                height=objEll->DoubleAttribute("height",100);
                                                radius=objEll->DoubleAttribute("radius",100);
                                                Obj.push_back(new GOAT::raytracing::Cone(Pos,radius,height,n));
                                                Obj[numObj]->setMatrix(alpha, beta, gamma);
												Obj[numObj]->setActive(isActive);
												S.addObject(Obj[numObj]);	
                                                break;
                                            }

					}
									
					numObj++;
				} // while loop

				  // S.addObjectList(numObj, Obj);

			}
			/* End of objects */
		}

		void xmlReader::doCalculations()
		{
           int type;
			std::string typeStr;
			tinyxml2::XMLElement* ell;
            std::cout << "doCalculations" << std::endl;
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
                        if (inactiveStr.compare("false")==0)
						{						
						typeStr = objEll->Attribute("Type");
                        std::cout << "typeStr=" << typeStr << std::endl;
						if (!typeStr.empty())
						{							
							type = mapString2CalculationToken(typeStr);
                            switch (type)
							{
                            case TOKEN_CALCULATION_PURE:
                            {
                             GOAT::raytracing::Raytrace_pure rt(S);
                             rt.trace();
                             break;
                            }
							case TOKEN_CALCULATION_PATH:
							{
                                std::cout << "do path calculation" << std::endl;
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
								std::string methodStr;
								const char* hStr;
								hStr=objEll->Attribute("Method");
								if (hStr != NULL) methodStr = hStr;
								else methodStr = "mixed";
								if (methodStr.compare("rtonly")==0)
										doPulseCalculation_rt(objEll); 
								else 
										doPulseCalculation(objEll);
								break;
							}

                            case TOKEN_CALCULATION_PULSE_FIELD:
                            {
                                std::string fname = objEll->Attribute("Filename");
                                if (!fname.empty())
                                {
                                    GOAT::raytracing::pulseCalculation_Field pc(S);
                                    GOAT::raytracing::TrafoParms trafoparms;
                                    trafoparms = pc.getTrafoParms();
                                    pc.setCenterWavelength(objEll->DoubleAttribute("Wavelength", trafoparms.wvl));
                                    pc.setNumReflex(objEll->IntAttribute("NumReflexions", trafoparms.nR));
                                    pc.setNumWavelengthsPerRange(objEll->IntAttribute("NumWavelengthsPerRange", trafoparms.nS));
                                    pc.setPulseWidth(objEll->DoubleAttribute("Pulse_width",trafoparms.dt));
                                    pc.setSpectralRanges(objEll->IntAttribute("NumSpectralRanges", trafoparms.nI));
                                    pc.setReferenceTime(objEll->IntAttribute("Reference_time", pc.getReferenceTime()));
                                    double repRate = objEll->DoubleAttribute("Repetition_rate", -1);
                                    if (repRate > 0) pc.setRepetitionRate(repRate);
                                    double dx = 2.0 * S.r0 / (double)pc.getNumCellsPerDirection();
                                    pc.setSpatialResolution(objEll->DoubleAttribute("Spatial_resolution", dx));
                                    double D=objEll->DoubleAttribute("D",-1.0);
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
                                        double offset = objEll->DoubleAttribute("Time_offset", 0);
                                        int objEstimate = objEll->IntAttribute("EstimateTimeForObject", 0);
                                        time = pc.findHitTime(objEstimate);
                                            time+= offset;
                                    }

                                    std::string fullfname;
                                    double d;
                                    if (D>0)
                                    {
                                        const char* hStr;
                                        std::string corrFilename;
                                        std::ofstream corrOS;
                                        hStr=objEll->Attribute("CorrelationFilename");
                                        if (hStr != NULL)
                                        {
                                            corrOS.open(hStr);
                                        }




//										do
                                        {
                                          pc.field(time);

                                          for (int i = 0; i < S.nObj; i++)
                                          {
                                            if (S.Obj[i]->isActive())
                                            {
                                                fullfname = fname + std::to_string(i) + ".dat";
                                                GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
                                            }
                                          }
                                          if (hStr != NULL) corrOS << d << std::endl;
                                        }// while (d>D);
                                      if (hStr != NULL) corrOS.close();
                                    }
                                    else
                                    {
                                        pc.field(time);
                                        for (int i = 0; i < S.nObj; i++)
                                          {
                                            if (S.Obj[i]->isActive())
                                            {
                                                fullfname = fname + std::to_string(i) + ".dat";
                                                GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
                                            }
                                          }
                                    }
                                }
                                else
                                    std::cerr << "Path calculation: You forgot to give an appropriate file name for the output!!" << std::endl;
                                break;
                            }
                            case TOKEN_CALCULATION_INELASTIC:
							{
								std::string fname = objEll->Attribute("Filename");
								if (fname.empty())
								{
									fname = "dummy";
								}
								int n = objEll->IntAttribute("n",500);
									GOAT::raytracing::Raytrace_Inel rt(S,n);									
									bool fieldonly=true;
								
									const char *str = objEll->Attribute("FieldOnly");
									if (str !=NULL)
									{
										fieldonly = (strcmp(str,"true") == 0);
									}
									
									
									if (fieldonly) rt.setExcitationFieldOnly();
									GOAT::raytracing::RRTParms rrtparms;
									
									rt.trace(rrtparms);
									std::string fullfname;
									rt.exportExcitation(fname, GOAT::raytracing::INEL_EXPORT_EXCITATION_FIELD_VECTOR);
									/*for (int i = 0; i < S.nObj; i++)
									{
										if (S.Obj[i]->Active)
										{
											fullfname = fname + std::to_string(i) + ".dat";
											GOAT::raytracing::saveFullE(rt.SGE[0], fullfname, i);
										}
									}*/
								
							}
							} // switch
							for (int i = 0; i < numDet; i++)
								S.Det[i]->save(Det[i]->fname.c_str());
						} // is type given ?
					} // is inactive ?
				} // for loop
			} // if no Calculations
		}

void xmlReader::doPulseCalculation(tinyxml2::XMLElement* objEll)
        {			
            std::cout << "------------------ DO PULSED CALCULATION  (mixed) -----------------" << std::endl;
            const char* hStr;
            std::string fname = objEll->Attribute("Filename");
            if (!fname.empty())
            {
                GOAT::raytracing::pulseCalculation pc(S);
                GOAT::raytracing::TrafoParms trafoparms;
                trafoparms = pc.getTrafoParms();
                pc.setCenterWavelength(objEll->DoubleAttribute("Wavelength", trafoparms.wvl));
                pc.setNumReflex(objEll->IntAttribute("NumReflexions", trafoparms.nR));
                pc.setNumWavelengthsPerRange(objEll->IntAttribute("NumWavelengthsPerRange", trafoparms.nS));
                pc.setPulseWidth(objEll->DoubleAttribute("Pulse_width",trafoparms.dt));
                pc.setSpectralRanges(objEll->IntAttribute("NumSpectralRanges", trafoparms.nI));
                pc.setReferenceTime(objEll->IntAttribute("Reference_time", pc.getReferenceTime()));
                pc.setNumberOfThreads(objEll->IntAttribute("NumberOfThreads",pc.getNumberOfThreads()));
                double repRate = objEll->DoubleAttribute("Repetition_rate", -1);
                if (repRate > 0) pc.setRepetitionRate(repRate);
                double dx = 2.0 * S.r0 / (double)pc.getNumCellsPerDirection();
				
                pc.setSpatialResolution(objEll->DoubleAttribute("Spatial_resolution", dx));
                double D=objEll->DoubleAttribute("D",-1.0);
                char cs[3];
                std::vector< std::function< std::complex< double >(double) > > nList;
                std::string refFuncName;
                bool failed = false;

                tinyxml2::XMLElement* refEll = objEll->FirstChildElement("RefractiveIndexList");
                if (refEll == NULL)
                {
                    std::cerr << "Pulse calculation: Refractive index function list is missing! Stopped calculation" << std::endl;
                    return;
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

                if (failed) return;


                hStr = refEll->Attribute("nS");
                if (hStr == NULL)
                {
                    std::cerr << "Pulse calculation: Refractive index function for surrounding medium is missing! Stopped calculation" << std::endl;
                    failed = true;
                }
                if (failed) return;
                refStr = hStr;
                refIndexToken = mapString2RefractiveIndexToken(refStr);
                if (refIndexToken == TOKEN_NOT_FOUND)
                {
                    std::cerr << "Pulse calculation: Wrong refractive index function name (" << refStr << ") !Calculation stopped!" << std::endl;
                    failed = true;
                }

                if (failed) return;
                addFunction2IndexList(nList, refIndexToken);
                pc.setRefractiveIndexFunctions(nList);

                double time = objEll->DoubleAttribute("Time", -1);
				std::cout << "time:" << time << std::endl;
               if (time < 0)
                {
                    double offset = objEll->DoubleAttribute("Time_offset", 0);
                    int objEstimate = objEll->IntAttribute("EstimateTimeForObject", 0);                    
                    time = pc.findHitTime(objEstimate);                    
                    std::cout << "estimated time: " << time << std::endl << std::flush;
                    time+= offset;
                }


                std::string fullfname;
                double d;
                if (D>0)
                {
                    const char* hStr;
                    std::string corrFilename;
                    std::ofstream corrOS;
                    hStr=objEll->Attribute("CorrelationFilename");
                    if (hStr != NULL)
                    {
                        corrOS.open(hStr);
                    }

                    int loopno=0;
                    do
                    {
                      d=pc.field(time,GOAT::raytracing::PULSECALCULATION_NOT_CLEAR_RESULT);								      
				      for (int i = 0; i < S.nObj; i++)
                      {
                        if (S.Obj[i]->isActive())
                        {
                            fullfname = fname + std::to_string(i) + ".dat";
                            GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);				
				        }
                      }
                      if (hStr != NULL) corrOS << d << std::endl;
                      loopno++;
                      std::cout << "loopno=" << loopno << std::endl;
                    } while (true); // while ( (d>D) || (loopno<2));
                  if (hStr != NULL) corrOS.close();
                }

                else
                {                    
                    pc.field(time);
                    for (int i = 0; i < S.nObj; i++)
                      {
                        if (S.Obj[i]->isActive())
                        {
                            fullfname = fname + std::to_string(i) + ".dat";
                            GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
                        }
                      }
                }
            }
            else
                std::cerr << "Path calculation: You forgot to give an appropriate file name for the output!!" << std::endl;
            return;
        }



        void xmlReader::doPulseCalculation_rt(tinyxml2::XMLElement* objEll)
        {			
            std::cout << "------------------ DO PULSED CALCULATION -----------------" << std::endl;
            const char* hStr;
            std::string fname = objEll->Attribute("Filename");
            if (!fname.empty())
            {
                int numLoops = objEll->IntAttribute("NumLoops", -1);
                GOAT::raytracing::pulseCalculation_rt pc(S);
                GOAT::raytracing::TrafoParms trafoparms;
                //trafoparms = pc.getTrafoParms();
                pc.setCenterWavelength(objEll->DoubleAttribute("Wavelength", trafoparms.wvl));
                pc.setNumReflex(objEll->IntAttribute("NumReflexions", trafoparms.nR));
                //pc.setNumWavelengthsPerRange(objEll->IntAttribute("NumWavelengthsPerRange", trafoparms.nS));
                pc.setPulseWidth(objEll->DoubleAttribute("Pulse_width",trafoparms.dt));
                pc.setSpectralRanges(objEll->IntAttribute("NumSpectralRanges", trafoparms.nI));
                //pc.setReferenceTime(objEll->IntAttribute("Reference_time", pc.getReferenceTime()));
                // pc.setNumberOfThreads(objEll->IntAttribute("NumberOfThreads",pc.getNumberOfThreads()));
                double repRate = objEll->DoubleAttribute("Repetition_rate", -1);
                if (repRate > 0) pc.setRepetitionRate(repRate);
                double dx = 2.0 * S.r0 / (double)pc.getNumCellsPerDirection();
				
                pc.setSpatialResolution(objEll->DoubleAttribute("Spatial_resolution", dx));
                double D=objEll->DoubleAttribute("D",-1.0);
                char cs[3];
                std::vector< std::function< std::complex< double >(double) > > nList;
                std::string refFuncName;
                bool failed = false;

                tinyxml2::XMLElement* refEll = objEll->FirstChildElement("RefractiveIndexList");
                if (refEll == NULL)
                {
                    std::cerr << "Pulse calculation: Refractive index function list is missing! Stopped calculation" << std::endl;
                    return;
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

                if (failed) return;


                hStr = refEll->Attribute("nS");
                if (hStr == NULL)
                {
                    std::cerr << "Pulse calculation: Refractive index function for surrounding medium is missing! Stopped calculation" << std::endl;
                    failed = true;
                }
                if (failed) return;
                refStr = hStr;
                refIndexToken = mapString2RefractiveIndexToken(refStr);
                if (refIndexToken == TOKEN_NOT_FOUND)
                {
                    std::cerr << "Pulse calculation: Wrong refractive index function name (" << refStr << ") !Calculation stopped!" << std::endl;
                    failed = true;
                }

                if (failed) return;
                addFunction2IndexList(nList, refIndexToken);
                pc.setRefractiveIndexFunctions(nList);

                double time = objEll->DoubleAttribute("Time", -1);
				std::cout << "time:" << time << std::endl;
               if (time < 0)
                {
                    double offset = objEll->DoubleAttribute("Time_offset", 0);
                    int objEstimate = objEll->IntAttribute("EstimateTimeForObject", 0);                    
//                    time = pc.findHitTime(objEstimate);                    
                    std::cout << "estimated time: " << time << std::endl << std::flush;
                    time+= offset;
                }


                std::string fullfname;
                double d;
               // if (D>0)
                {
                    const char* hStr;
                    std::string corrFilename;
                    std::ofstream corrOS;
                    hStr=objEll->Attribute("CorrelationFilename");
                    if (hStr != NULL)
                    {
                        corrOS.open(hStr);
                    }

                    int loopno=0;
                    bool cancel=false;
                    do
                    {
                //      d=pc.field(time,GOAT::raytracing::PULSECALCULATION_NOT_CLEAR_RESULT);								      
						pc.field(time);
				
                      for (int i = 0; i < S.nObj; i++)
                      {
                        if (S.Obj[i]->isActive())
                        {
                            fullfname = fname + std::to_string(i) + ".dat";
                            // GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
							
		        	GOAT::raytracing::saveFullE(pc.rt.SA[0], fullfname, i);
					
					d=sumabs2(pc.rt.SA[0],i);
				        }
                      }
                      if (hStr != NULL) corrOS << d << std::endl;
                      loopno++;
                      cancel = (loopno >= numLoops) && (numLoops >= 0);
                      std::cout << "loopno=" << loopno << std::endl;
                    } while (!cancel); // while ( (d>D) || (loopno<2));
                  if (hStr != NULL) corrOS.close();
                }
/*
                else
                {                    
                    pc.field(time);
                    for (int i = 0; i < S.nObj; i++)
                      {
                        if (S.Obj[i]->isActive())
                        {
                            fullfname = fname + std::to_string(i) + ".dat";
                            GOAT::raytracing::saveFullE(pc.rt.SA[0], fullfname, i);
                        }
                      }
                }
*/
            }
            else
                std::cerr << "Path calculation: You forgot to give an appropriate file name for the output!!" << std::endl;
            return;
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

        std::complex<double> xmlReader::readCmplx(tinyxml2::XMLElement* ell, int &xmlError)
        {
            double re,im;
            if (ell != NULL)
            {
                xmlError=ell->QueryDoubleAttribute("real",&re);
                xmlError = xmlError | ell->QueryDoubleAttribute("imag", &im);
            }
            return std::complex<double> (re,im);
        }

        GOAT::maths::Vector<std::complex<double> > xmlReader::readCmplxVector (tinyxml2::XMLElement* ell, int& xmlError)
        {
          std::complex<double> x,y,z;
          int errorX, errorY, errorZ;
          if (ell != NULL)
          {
              x=readCmplx(ell,errorX);
              y=readCmplx(ell,errorY);
              z=readCmplx(ell,errorZ);
              xmlError = errorX | errorY | errorZ;
          }
          return GOAT::maths::Vector<std::complex<double> > (x,y,z);
        }

        /*------------------------------- XML-Writer Implementation ----------------------------------------- */
        xmlWriter::xmlWriter(GOAT::raytracing::Scene S)
        {
            this->S=S;
        }

        void xmlWriter::write (std::string fname)
        {
          root=doc.NewElement("Root");
          doc.InsertFirstChild(root);
          scene=doc.NewElement("Scene");


        }
	}
}
