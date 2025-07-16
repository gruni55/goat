

#include "tinyxml2.h"
#include "xml.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation.h"
#include "pulsecalculation_rt.h"
#include "pulsecalculation_field.h"
#include "raytrace_inel.h"
#include <chrono>
#include <goodies.h>
#include <filesystem>


#define tl(s) GOAT::maths::tl(s)	

namespace GOAT
{
	namespace XML
	{

        bool findExtension (std::string fname, std::string extension)
        {
            std::size_t ifound = fname.find_last_of (".");
            std::string ext=fname.substr(ifound+1,std::string::npos);

            bool found = ext.compare (extension)==0;
            return found;
        }

		// void xmlReader::readXML(const char* fname, char* path)
        void xmlReader::readXML(std::string fname, bool calc_enabled, std::string path)
		{
            setEnableCalculation(calc_enabled);
            std::ofstream os;
            os.open("test.log");
            os << "Filename :" << fname << std::flush << std::endl;
            os.close();
            this->path = path;
			setlocale(LC_NUMERIC, "C");
			tinyxml2::XMLDocument doc;
            std::cout << "fname=" << fname << "\tpath=" << path << std::endl;
            // check, if path is given separatly 
            if (path.size() >0)
            {
                std::string fstr = std::string(path) + "/" + std::string(fname);
                fname = fstr.c_str();
            }
            else // path is not given separatly => try to extract it from the filename
            {
                
                std::filesystem::path p(fname);
                std::filesystem::path dir=p.parent_path();
                std::filesystem::path filename=p.filename();
                if (p.is_absolute())
                {
                    dir = ""; 
                    std::cout << "path is absolute" << std::endl;
                }
                path=dir.string();
               // fname=filename.string();
            }

			tinyxml2::XMLError eResult=doc.LoadFile(fname.c_str());                        
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

                if (calculation_enabled) doCalculations();
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
					typeStr = detEll->Attribute("type");                    
					std::string filename;
					filename = detEll->Attribute("filename");

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
                                                    Det[numDet]->load(filename.c_str());
                                                    numDet++;
                                                    std::cout << "Detector filename=" << filename << std::endl;
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
               GOAT::maths::Vector<std::complex<double> > Pol;                
               GOAT::maths::Vector<double> Pold;
			   GOAT::maths::Vector<double> Pos;
				int numRays;
				double wavelength;
				double size;
                double power;
				GOAT::raytracing::LightSrc* ls = NULL;
				for (tinyxml2::XMLElement* lsEll = ell->FirstChildElement("LightSource"); lsEll != NULL; lsEll = lsEll->NextSiblingElement("LightSource"))
				{
			        	std::string typeStr;
					typeStr = lsEll->Attribute("type");
					Pos = readVector(lsEll->FirstChildElement("Position"));
					numRays = lsEll->IntAttribute("numRays", 100);
					wavelength = lsEll->DoubleAttribute("wavelength", 1.0);
                    std::cout << "read: Wavelength:" << wavelength << std::endl;
					size = lsEll->DoubleAttribute("size", 10.0);
                 
					Pold=readVector(lsEll->FirstChildElement("Polarisation"),1,0,0);
                 
                    power = lsEll->DoubleAttribute("power", 1.0);
                    Pol[0]=Pold[0];
                    Pol[1]=Pold[1];
                    Pol[2]=Pold[2];
					int type = mapString2LightSourceToken(typeStr);
                    

                   switch (type)
					{
					case TOKEN_LIGHTSOURCE_PLANE: {
													ls = new GOAT::raytracing::LightSrcPlane(Pos, numRays, wavelength, size, Pol);
													GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"),0,0,1);
													ls->setk(k);
                                                    LS.push_back(ls);
												   }
												   break;

					case TOKEN_LIGHTSOURCE_PLANE_MC: {                                                              
														ls = new GOAT::raytracing::LightSrcPlane_mc(Pos, numRays, wavelength, size, Pol);
														GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
														ls->setk(k);
                                                        LS.push_back(ls);
													 }
                                                      break;
                    case TOKEN_LIGHTSOURCE_LINE: {
                                                        
                                                        GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
                                                        GOAT::maths::Vector<double> D = readVector(lsEll->FirstChildElement("lateral_direction"));
                                                        ls = new GOAT::raytracing::LightSrcLine(Pos, numRays, wavelength, size, k);
                                                         ls->setPol(Pol);
                                                         LS.push_back(ls);
                                                        }
                                               break;
                    case TOKEN_LIGHTSOURCE_LINE_MC: {
                                                        GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
                                                        GOAT::maths::Vector<double> D = readVector(lsEll->FirstChildElement("lateral_direction"));
                                                        ls = new GOAT::raytracing::LightSrcLine_mc(Pos, numRays, wavelength, size, k, D);
                                                        ls->setPol(Pol);
                                                        LS.push_back(ls);
                                                    }
                                               break;

                    case TOKEN_LIGHTSOURCE_POINT_MC: {
                                                        ls = new GOAT::raytracing::LightSrcPoint_mc(Pos, numRays, wavelength);
                                                        LS.push_back(ls);
                                                     }
                                                   break;

					case TOKEN_LIGHTSOURCE_GAUSSIAN: {
														double w0;
														double NA=1.0;
														tinyxml2::XMLError err;
														w0 = lsEll->DoubleAttribute("w0", 1.0);
                                                        size=lsEll->DoubleAttribute("size",1.0);
														GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
														ls = new GOAT::raytracing::LightSrcGauss(Pos, numRays, wavelength, w0, focusPos,size,Pol);
														err = lsEll->QueryDoubleAttribute("NA", &NA);
														if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss*)ls)->setNA(NA);
                                                        LS.push_back(ls);
													  }
													 	  break;

					case TOKEN_LIGHTSOURCE_GAUSSIAN_MC: {
															double w0;
															double NA = 1.0;
															tinyxml2::XMLError err;
															w0 = lsEll->DoubleAttribute("w0", 1.0);
                                                            size=lsEll->DoubleAttribute("size",1.0);
															GOAT::maths::Vector<double> focusPos = readVector(lsEll->FirstChildElement("FocusPosition"));
															ls = new GOAT::raytracing::LightSrcGauss_mc(Pos, numRays, wavelength, w0, focusPos,size,Pol);
															err = lsEll->QueryDoubleAttribute("NA", &NA);
															if (err == tinyxml2::XML_SUCCESS) ((GOAT::raytracing::LightSrcGauss_mc*)ls)->setNA(NA);
                                                            LS.push_back(ls);
														}
													  break;
                    case TOKEN_LIGHTSOURCE_RING:        {
                                                         double rmin, rmax;
                                                         rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                         rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                         ls=new GOAT::raytracing::LightSrcRing(Pos, numRays, wavelength, rmin,rmax,Pol);
														 GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
														 ls->setk(k);
                                                         LS.push_back(ls);
                                                        break;
                                                        }
                    case TOKEN_LIGHTSOURCE_RING_MC:
                                                        {
                                                          double rmin, rmax;
                                                          rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                          rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                           ls=new GOAT::raytracing::LightSrcRing_mc(Pos, numRays, wavelength, rmin,rmax,Pol);
														   GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"));
                                                           std::cout << "1 Pold=" << ls->Pol << "\t" << ls->initPol << std::endl;
														   ls->setk(k);
                                                           std::cout << "2 Pold=" << ls->Pol << "\t" << ls->initPol << std::endl;
                                                           LS.push_back(ls);
														   break;
                                                        }
                    case TOKEN_LIGHTSOURCE_GAUSSIAN_RING_MC :
                                                        {
                                                            std::cout << "light source gaussian ring mc" << std::endl;
                                                         double rmin, rmax;
                                                         double width;
                                                         rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                         rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                         width=lsEll->DoubleAttribute("width",rmax);
                                                           ls=new GOAT::raytracing::LightSrcRingGauss_mc(Pos, numRays, wavelength, rmin, rmax,Pol);
                                                         ((GOAT::raytracing::LightSrcRingGauss_mc *)ls)->setFWHM(width);
														 GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"),0,0,1);
														 ls->setk(k);
                                                         LS.push_back(ls);
                                                         break;
                                                        }


					}
                   LS[numLS]->P0 = power;
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
					typeStr = calcEll->Attribute("type");
					fname = calcEll->Attribute("filename");
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
					typeStr = objEll->Attribute("type");
                    alpha = objEll->DoubleAttribute("alpha", 0.0) / 180.0 * M_PI;
					beta = objEll->DoubleAttribute("beta", 0.0) / 180.0 * M_PI;
					gamma = objEll->DoubleAttribute("gamma", 0.0) / 180.0 * M_PI;
					isActive = objEll->BoolAttribute("isactive", false);
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
											fileTypeStr = objEll->Attribute("filetype");

											if (fileTypeStr.compare(".srf") == 0)
											{
												fileName = objEll->Attribute("filename");
                                          
                                                std::cout << "path.size()=" << path.size() << std::endl;
                                                if (path.size()>0)
                                                {
                                                    std::string sep = "/";
                                                    std::filesystem::path p(fileName);
                                                    if (p.is_relative())
                                                    fileName = path + sep + fileName;
                                                    std::cout << "fileName:" << fileName << std::endl;
                                                }

                                                

												if (!fileName.empty()) ((GOAT::raytracing::surface*)Obj[numObj])->createsurface(fileName);
											}

											if (fileTypeStr.compare(".stl") == 0)
											{
												fileName = objEll->Attribute("filename");
                                                std::cout << "path.size()=" << path.size() << std::endl;
                                                if (path.size() > 0)
                                                {
                                                    std::filesystem::path p(fileName);
                                                    if (p.is_relative())
                                                    {
                                                        std::string sep = "/";
                                                        fileName = path + sep + fileName;
                                                    }
                                                    std::cout << "fileName:" << fileName << std::endl;
                                                }
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
												else lensparms.left.curvature = GOAT::raytracing::flat; // default value is "flat"
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
                                            	lensparms.offset = objEll->DoubleAttribute("offset", 0.0);
												lensparms.radius = objEll->DoubleAttribute("radius", 0.0);                                                


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

                    case TOKEN_OBJECT_CYLINDER:
                                            {                                                
                                                double height, radius;
                                                height = objEll->DoubleAttribute("height", 1);
                                                radius = objEll->DoubleAttribute("radius", 1);                                                
                                                Obj.push_back(new GOAT::raytracing::Cylinder(Pos, radius, height, n));
                                                Obj[numObj]->setMatrix(alpha, beta, gamma);
                                                Obj[numObj]->setActive(isActive);
                                                S.addObject(Obj[numObj]);
                                                break;
                                            }
                    case TOKEN_OBJECT_VORTEX_PLATE:
                                            {                                                
                                                double height, radius, dh;
                                                int m;
                                                height = objEll->DoubleAttribute("height", 1);
                                                radius = objEll->DoubleAttribute("radius", 1);
                                                m = objEll->IntAttribute("m", 1);
                                                dh = objEll->DoubleAttribute("dh", 1);
                                                Obj.push_back(new GOAT::raytracing::VortexPlate(Pos, radius, height, dh, m, n));
                                                Obj[numObj]->setMatrix(alpha, beta, gamma);
                                                Obj[numObj]->setActive(isActive);
                                                S.addObject(Obj[numObj]);
                                                break;
                            
                                            }

					}
					double sf=objEll->DoubleAttribute("scaling",1);
                    if ((sf!=1) && (sf>0)) Obj[numObj]->scale(sf); 
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
						hStr=objEll->Attribute("inactive");
						if (hStr != NULL) inactiveStr = hStr;
						else inactiveStr = "false";
                        if (inactiveStr.compare("false")==0)
						{						
						typeStr = objEll->Attribute("type");
                        std::cout << "typeStr=" << typeStr << std::endl;

                        // change number of rays, if given
                        int numRays;
                        std::vector<int> numRays_old;
                        numRays = objEll->IntAttribute("numRays", 0);
                        bool numRaysChanged=false;
                        if (numRays > 0)
									{
										// store the old values 									 
                                        numRaysChanged=true;
										for (int i = 0; i < S.nLS; i++)
										{
											numRays_old.push_back(S.LS[i]->getNumRays());
											S.LS[i]->setNumRays(numRays);
										}
									}

                        int numReflex;
                        numReflex = objEll->IntAttribute("numReflex", 0);

						if (!typeStr.empty())
						{							
							type = mapString2CalculationToken(typeStr);
                            switch (type)
							{
                            case TOKEN_CALCULATION_PURE:
                            {
                             GOAT::raytracing::Raytrace_pure rt(S);      
                             rt.setNumReflex(numReflex);                       
                             rt.trace();
                             break;
                            }
							case TOKEN_CALCULATION_PATH:
							{
                                std::cout << "do path calculation" << std::endl;
                                std::string fname = objEll->Attribute("filename");								
								int numDet=S.nDet;
                                // S.nDet=0;
								if (!fname.empty())
								{
									GOAT::raytracing::Raytrace_Path rt(S);
									rt.setNumReflex(numReflex);
									rt.trace(fname);
								}
								else
									std::cerr << "Path calculation: You forgot to give an appropriate file name for the output!!" << std::endl;
                                // S.nDet=numDet;
								break;
							} // case path calculation

                            case TOKEN_CALCULATION_PULSE: 
							{
								std::string methodStr;
								const char* hStr;
								hStr=objEll->Attribute("method");
								if (hStr != NULL) methodStr = hStr;
								else methodStr = "rtonly";
								if (methodStr.compare("rtonly")==0)
										doPulseCalculation_rt(objEll); 
								else 
										doPulseCalculation(objEll);
								break;
							}

                            case TOKEN_CALCULATION_PULSE_FIELD:
                            {
                                std::string fname = objEll->Attribute("filename");
                                if (!fname.empty())
                                {
                                    GOAT::raytracing::pulseCalculation_Field pc(S);
                                    GOAT::raytracing::TrafoParms trafoparms;
                                    trafoparms = pc.getTrafoParms();
                                    pc.setCenterWavelength(objEll->DoubleAttribute("wavelength", trafoparms.wvl));
                                    pc.setNumReflex(numReflex);
                                    trafoparms.nR=numReflex;
                                    // pc.setNumReflex(objEll->IntAttribute("NumReflex", trafoparms.nR));
                                    pc.setNumWavelengthsPerRange(objEll->IntAttribute("NumWavelengthsPerRange", trafoparms.nS));
                                    pc.setPulseWidth(objEll->DoubleAttribute("pulseWidth",trafoparms.dt));
                                    pc.setSpectralRanges(objEll->IntAttribute("numSpectralRanges", trafoparms.nI));
                                    pc.setReferenceTime(objEll->IntAttribute("referenceTime", pc.getReferenceTime()));                                    
                                    double repRate = objEll->DoubleAttribute("repetitionRate", -1);
                                    if (repRate > 0) pc.setRepetitionRate(repRate);
                                    double dx = 2.0 * S.r0 / (double)pc.getNumCellsPerDirection();
                                    pc.setSpatialResolution(objEll->DoubleAttribute("spatialResolution", dx));
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

                                    double time = objEll->DoubleAttribute("time", -1);
                                    if (time < 0)
                                    {
                                        double offset = objEll->DoubleAttribute("timeOffset", 0);
                                        int objEstimate = objEll->IntAttribute("estimateTimeForObject", 0);
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
                                        hStr=objEll->Attribute("correlationFilename");
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
								std::string fname = objEll->Attribute("filename");
								if (fname.empty())
								{
									fname = "dummy";
								}
								int n = objEll->IntAttribute("n",500);
                                S.setNumberOfCellsPerDirection(n);
									GOAT::raytracing::Raytrace_Inel rt(S);									
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
                            std::cout << "number of detectors:" << numDet << std::endl;
                            double normfac = 0;
                            double Iall = 0;
                            for (int i = 0; i < numLS; i++)
                            {
                                normfac += LS[i]->getNumRays();
                                // normfac += LS[i]->P0 * LS[i]->getNumRays() * raytracing::mu0 * raytracing::c_light / (LS[i]->area() * 1E-12 * LS[i]->getIsum1());  // factor 1E-12 to convert �m^2 into m^2 
                                Iall += LS[i]->getIsum1();
                                std::cout.precision (6);
                                std::cout << "area=" << LS[i]->area() << "\tIsum1=" << LS[i]->getIsum1() << std::endl;
                            }
                            
                            // normfac /= Iall;
                            normfac =  sqrt(normfac) ;                            
                            std::cout << "Iall=" << Iall << "\tnormfac=" << normfac << std::endl;
                            for (int i = 0; i < numDet; i++)
                            {
                                S.Det[i]->mult(normfac);
                                S.Det[i]->save(Det[i]->fname.c_str());
                            }
                            if (numRaysChanged)
									{
										// restore the old values  
										for (int i = 0; i < S.nLS; i++)
											S.LS[i]->setNumRays(numRays_old[i]);
									} 
						} // is type given ?
					} // is inactive ?
				} // for loop
			} // if no Calculations
		}

void xmlReader::doPulseCalculation(tinyxml2::XMLElement* objEll)
        {			
            std::cout << "------------------ DO PULSED CALCULATION  (mixed) -----------------" << std::endl;
            const char* hStr;
            std::string fname = objEll->Attribute("filename");
            if (!fname.empty())
            {
                GOAT::raytracing::pulseCalculation pc(S);
                GOAT::raytracing::TrafoParms trafoparms;
                trafoparms = pc.getTrafoParms();
                pc.setCenterWavelength(objEll->DoubleAttribute("wavelength", trafoparms.wvl));
                pc.setNumReflex(objEll->IntAttribute("numReflex", trafoparms.nR));
                pc.setNumWavelengthsPerRange(objEll->IntAttribute("NumWavelengthsPerRange", trafoparms.nS));
                pc.setPulseWidth(objEll->DoubleAttribute("pulseWidth",trafoparms.dt));
                pc.setSpectralRanges(objEll->IntAttribute("numSpectralRanges", trafoparms.nI));
                pc.setReferenceTime(objEll->IntAttribute("referenceTime", pc.getReferenceTime()));
                pc.setNumberOfThreads(objEll->IntAttribute("NumberOfThreads",pc.getNumberOfThreads()));
                double repRate = objEll->DoubleAttribute("repetitionRate", -1);
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
            std::string fname = objEll->Attribute("filename");
            if (!fname.empty())
            {
                int numLoops = objEll->IntAttribute("numLoops", -1);
                GOAT::raytracing::pulseCalculation_rt pc(S);
                GOAT::raytracing::TrafoParms trafoparms;
                //trafoparms = pc.getTrafoParms();
                pc.setCenterWavelength(objEll->DoubleAttribute("wavelength", trafoparms.wvl));
                pc.setNumReflex(objEll->IntAttribute("numReflex", trafoparms.nR));
                //pc.setNumWavelengthsPerRange(objEll->IntAttribute("numWavelengthsPerRange", trafoparms.nS));
                pc.setPulseWidth(objEll->DoubleAttribute("pulseWidth",trafoparms.dt));
                pc.setSpectralRanges(objEll->IntAttribute("numSpectralRanges", trafoparms.nI));
                //pc.setReferenceTime(objEll->IntAttribute("Reference_time", pc.getReferenceTime()));
                // pc.setNumberOfThreads(objEll->IntAttribute("NumberOfThreads",pc.getNumberOfThreads()));
                double repRate = objEll->DoubleAttribute("repetitionRate", -1);
                if (repRate > 0) pc.setRepetitionRate(repRate);
                double dx = 2.0 * S.r0 / (double)pc.getNumCellsPerDirection();
				
                pc.setSpatialResolution(objEll->DoubleAttribute("spatialResolution", dx));
                 std::cout << "dx=" << dx << std::endl;
                                   
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

                double time = objEll->DoubleAttribute("time", -1);
				std::cout << "time:" << time << std::endl;
               if (time < 0)
                {
                    double offset = objEll->DoubleAttribute("timeOffset", 0);
                    int objEstimate = objEll->IntAttribute("estimateTimeForObject", 0);                    
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
                    hStr=objEll->Attribute("correlationFilename");
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
                       //     GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
							
		        	 GOAT::raytracing::saveFullE(pc.rt.SA[0], fullfname, i);
					
					           d=sumabs2(pc.rt.SA[0],i);
                             std::cout << "d=" << d << std::endl;
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
            bool cartesian=false;
            int error;
			if (ell != NULL)
			{
				xmlError = ell->QueryDoubleAttribute("x", &x);
                cartesian = xmlError==tinyxml2::XML_SUCCESS;
                error = ell->QueryDoubleAttribute("y", &y);
                cartesian = cartesian | (error==tinyxml2::XML_SUCCESS);
				xmlError = xmlError | error;
				error = ell->QueryDoubleAttribute("z", &z);
                cartesian = cartesian | (error==tinyxml2::XML_SUCCESS);
				xmlError = xmlError | error;
                // use caresian coordinates, if x,y or z is given 
                // otherwise try spherical coordinates                 
                if (!cartesian) 
                {
                    double r,theta,phi;
                   xmlError = ell->QueryDoubleAttribute("r", &r); 
                   xmlError = xmlError | ell->QueryDoubleAttribute("theta", &theta);
                   xmlError = xmlError | ell->QueryDoubleAttribute("phi", &phi);
                   if (xmlError == tinyxml2::XML_SUCCESS) 
                   {
                    x=r * cos(phi) * sin(theta);
                    y=r * sin(phi) * sin(theta);
                    z=r * cos(theta);
                   }
                }
			}
			return GOAT::maths::Vector<double>(x,y,z);
		}


		std::complex<double> xmlReader::readCmplx(tinyxml2::XMLElement* ell, double defre, double defim)
		{
			double im=defim, re=defre;
			if (ell != NULL)
			{       
                /*         
				re = ell->DoubleAttribute("real", defre);
				im = ell->DoubleAttribute("imag", defim);
                */              
                ell->QueryDoubleAttribute("real", &re);
				ell->QueryDoubleAttribute("imag", &im);
			}
			return std::complex<double>(re, im);
		}

        std::complex<double> xmlReader::readCmplx(tinyxml2::XMLElement* ell, int &xmlError)
        {
            double re=0,im=0;
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
            tinyxml2::XMLDeclaration* decl = doc.NewDeclaration(R"(xml version="1.0" encoding="utf-8")");
            doc.InsertFirstChild(decl);
            std::cout << "write :" << fname << std::endl;
          root=doc.NewElement("Root");
          doc.InsertEndChild(root);
          scene=doc.NewElement("Scene");
          scene->SetAttribute("r0", formatDouble(S.r0).c_str());
          root->InsertEndChild(scene);
          std::cout << "no. of light sources: "<< S.nLS << std::endl;
          if (S.nLS > 0)
          {
            std::cout << "write light sources" << std::endl;
              tinyxml2::XMLElement* lightSrc;
              lightSrcs = doc.NewElement("LightSources");
               for (int i = 0; i < S.nLS; i++)
                   writeLightSrc(i);              
               scene->InsertEndChild(lightSrcs);
          }

          if (S.nObj > 0)
          {
               objects=doc.NewElement("Objects");
               for (int i=0; i<S.nObj; i++)
                    writeObject(i);
                scene->InsertEndChild(objects);
          }

          if (S.nDet > 0)
          {
            detectors=doc.NewElement("Detectors");
            for (int i=0; i<S.nDet; i++)
                writeDetector(i);
            scene->InsertEndChild(detectors);
          }
          
          tinyxml2::XMLError e = doc.SaveFile(fname.c_str());


        }
        void xmlWriter::writeLightSrc(int i)
        {
            
            auto lightSrc = doc.NewElement("LightSource");
            int type = S.LS[i]->type;
            int typeh = type < 10 ? type-1 : type - 5;
            std::cout << "typeh=" << typeh << "\ttype=" << type << std::endl;
            lightSrc->SetAttribute("type", LSTYPES[typeh].c_str());
            lightSrc->SetAttribute("numRays", S.LS[i]->getNumRays());
            lightSrc->SetAttribute("wavelength", formatDouble(S.LS[i]->getWavelength()).c_str());

            lightSrc->InsertEndChild(writeVectorD("Position", S.LS[i]->Pos));
            lightSrc->InsertEndChild(writeVectorC("Polarisation", S.LS[i]->Pol));
            
            switch (type)
            {
              case raytracing::LIGHTSRC_SRCTYPE_PLANE:
              case raytracing::LIGHTSRC_SRCTYPE_PLANE_MC:
                  {
                  raytracing::LightSrcPlane* ls = (raytracing::LightSrcPlane*)S.LS[i];
                   lightSrc->InsertEndChild(writeVectorD("Direction", ls->getk()));
                   lightSrc->SetAttribute("size", formatDouble(ls->D).c_str());                   
                  }
                  break;

              case raytracing::LIGHTSRC_SRCTYPE_GAUSS :
              case raytracing::LIGHTSRC_SRCTYPE_GAUSS_MC:
              {
                  raytracing::LightSrcGauss* ls = (raytracing::LightSrcGauss*)S.LS[i];
                  lightSrc->SetAttribute("size", formatDouble(ls->D).c_str());
                  lightSrc->SetAttribute("w0", formatDouble(ls->w0).c_str());
              }
              break;

              case raytracing::LIGHTSRC_SRCTYPE_RING :
              case raytracing::LIGHTSRC_SRCTYPE_RING_MC:
              {
                  raytracing::LightSrcRing* ls = (raytracing::LightSrcRing*)S.LS[i];
                  lightSrc->SetAttribute("rmin", formatDouble(ls->getRmin()).c_str());
                  lightSrc->SetAttribute("rmax", formatDouble(ls->getRmax()).c_str());
                  lightSrc->InsertEndChild(writeVectorD("Direction", ls->getk()));
              }
            }
            lightSrcs->InsertEndChild(lightSrc);
        }

        void xmlWriter::writeObject(int i)
        {
            auto object = doc.NewElement("Object");
            int typeh = S.Obj[i]->type-10000;
            int type = S.Obj[i]->type;


            // ---------------- global parameters ----------------
            object->SetAttribute("type",objectToken[typeh].c_str());                        
            object->InsertEndChild(writeVectorD("Position", S.Obj[i]->P));
            object->SetAttribute("alpha",formatDouble(S.Obj[i]->Ealpha/M_PI*180.0).c_str());
            object->SetAttribute("beta",formatDouble(S.Obj[i]->Ebeta/M_PI*180.0).c_str());
            object->SetAttribute("gamma",formatDouble(S.Obj[i]->Egamma/M_PI*180.0).c_str());
            object->SetAttribute("isactive",S.Obj[i]->isActive());
            object->InsertEndChild(writeComplex("n",S.Obj[i]->n));            
            object->SetAttribute("scaling",formatDouble(S.Obj[i]->sf).c_str());

            // --------------- special parameters ----------------
            switch (type) 
            {
                case OBJECTSHAPE_ELLIPSOID : 
                    {
                    auto obj=(raytracing::Ellipsoid *) S.Obj[i];
                        object->InsertEndChild(writeVectorD("Dimension",obj->r));
                    }
                    break;
                case OBJECTSHAPE_SURFACE : 
                    {
                    auto obj=(raytracing::surface *) S.Obj[i];
                        std::string ft;
                         switch (obj->filetype)
                         {
                            case OBJECTSHAPE_SURFACE_FILETYPE_STL : ft=".stl"; break;
                            case OBJECTSHAPE_SURFACE_FILETYPE_SRF : ft=".srf"; break;                            
                         }

                        if (obj->filetype!=OBJECTSHAPE_SURFACE_FILETYPE_NONE) object->SetAttribute("filetype",ft.c_str());
                        object->SetAttribute("filename",obj->getFilename().c_str());
                    }
                    break;
                case OBJECTSHAPE_CONE : break;
                case OBJECTSHAPE_ASPHERIC_LENS : break;
                case OBJECTSHAPE_SPHERIC_LENS : 
                    {
                        auto obj=(raytracing::sphericLens *) S.Obj[i];
                        raytracing::lensParms lensparms = obj->getParms();        
                        std::cout << "-> radius=" << lensparms.radius << std::endl;                                        
                        object->SetAttribute("radius",formatDouble(lensparms.radius).c_str());
                        object->SetAttribute("offset",formatDouble(lensparms.offset).c_str());

                        auto left= doc.NewElement("left");
                        switch (lensparms.left.curvature)
                        {
                            case raytracing::convex : left->SetAttribute("Curvature","convex"); break;                           
                            case raytracing::concave : left->SetAttribute("Curvature","concave"); break; 
                            case raytracing::flat : left->SetAttribute("Curvature","flat"); break;
                        }
                        left->SetAttribute("R",formatDouble(lensparms.left.R).c_str());
                        object->InsertEndChild(left);

                        auto right= doc.NewElement("right");
                        switch (lensparms.right.curvature)
                        {
                            case raytracing::convex : right->SetAttribute("Curvature","convex"); break;                           
                            case raytracing::concave : right->SetAttribute("Curvature","concave"); break; 
                            case raytracing::flat : right->SetAttribute("Curvature","flat"); break;
                        }
                        right->SetAttribute("R",formatDouble(lensparms.right.R).c_str());
                        object->InsertEndChild(right);                        
                    }
                    break;

                case OBJECTSHAPE_BOX: 
                    {
                    auto obj=(raytracing::Box *) S.Obj[i];
                        object->InsertEndChild(writeVectorD("Dimension",obj->d));
                    }
                    break;
                case OBJECTSHAPE_CYLINDER: 
                    {
                    auto obj=(raytracing::Cylinder *) S.Obj[i];
                        object->SetAttribute("height",formatDouble(obj->height()).c_str());
                        object->SetAttribute("radius",formatDouble(obj->radius()).c_str());
                    }
                    break;                
                case OBJECTSHAPE_VORTEX_PLATE: 
                    {
                        auto obj=(raytracing::VortexPlate *) S.Obj[i];
                        object->SetAttribute("height",formatDouble(obj->height()).c_str());
                        object->SetAttribute("radius",formatDouble(obj->radius()).c_str());
                        object->SetAttribute("m",obj->order());
                        object->SetAttribute("dh",formatDouble(obj->vortexHeight()).c_str());
                    }
                break;
            }
            objects->InsertEndChild(object);
            
        }

        void xmlWriter::writeDetector(int i)
        {
            auto detector = doc.NewElement("Detector");
            int type=S.Det[i]->Type();
            int typeh=type-20000;
            detector->SetAttribute("type",detectorToken[typeh].c_str());
            detector->InsertEndChild(writeVectorD("Position",S.Det[i]->position()));
            detector->InsertEndChild(writeVectorD("Direction",S.Det[i]->norm()));
            detector->SetAttribute("filename",S.Det[i]->fname.c_str());
            S.Det[i]->save(S.Det[i]->fname.c_str());
            switch (type)
            {
                case DETECTOR_PLANE : 
                    {
                     auto det=(raytracing::DetectorPlane *) S.Det[i];
                     detector->SetAttribute("d",formatDouble(det->D1()).c_str()); // we assume, that d1=d2 
                     detector->SetAttribute ("n", det->N1()); // we also assume that n1=n2                            
                    }
            }
            detectors->InsertEndChild(detector);
        }

        tinyxml2::XMLElement* xmlWriter::writeVectorD(std::string name, maths::Vector<double> v)
        {
            auto vell = doc.NewElement(name.c_str());
            vell->SetAttribute("x", formatDouble(v[0]).c_str());
            vell->SetAttribute("y", formatDouble(v[1]).c_str());
            vell->SetAttribute("z", formatDouble(v[2]).c_str());
            return vell;
        }

        tinyxml2::XMLElement* xmlWriter::writeVectorC(std::string name, maths::Vector<std::complex<double> > v)
        {
            auto vell = doc.NewElement(name.c_str());
            vell->InsertEndChild(writeComplex("x", v[0]));
            vell->InsertEndChild(writeComplex("y", v[1]));
            vell->InsertEndChild(writeComplex("z", v[2]));

            return vell;
        }
        tinyxml2::XMLElement* xmlWriter::writeComplex(std::string name, std::complex<double> z)
        {
            auto cell = doc.NewElement(name.c_str());
            cell->SetAttribute("real", formatDouble(real(z)).c_str());
            cell->SetAttribute("imag", formatDouble(imag(z)).c_str());
            return cell;
        }


	}
}
