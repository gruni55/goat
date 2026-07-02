

#include "tinyxml2.h"
#include "xml.h"
#include "lens.h"
#include "sphericLens.h"
#include "pulsecalculation.h"
#include "pulsecalculation_rt.h"
#include "pulsecalculation_field.h"
#include "raytrace_inel.h"
#include "kirchhoff.h"
#include "angularSpectrum.h"
#include "detector.h"
#include "goat_defines.h"
#include "roughObject.h"    
#include <chrono>
#include <goodies.h>
#include <filesystem>
#include "refractive_index_functions.h"
#include <sstream>
#include <iomanip>
#include <locale>
#include <cmath>

#define tl(s) GOAT::maths::tl(s)	

namespace GOAT
{
	namespace XML
	{
        inline std::string formatDouble(double value, int precision = 17)
        {
            // Optional: NaN / Inf explizit behandeln
            if (std::isnan(value)) return "nan";
            if (std::isinf(value)) return (value > 0) ? "inf" : "-inf";

            std::ostringstream oss;
            oss.imbue(std::locale::classic());   // erzwingt '.' als Dezimaltrennzeichen
            oss << std::setprecision(precision) << value;
            return oss.str();
        }

        void createXMLElementWithParam(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* parent, const calculationParam& param)
        {
            tinyxml2::XMLElement* paramElement = doc.NewElement("Param");
            paramElement->SetAttribute("name", param.name.c_str());
            if (std::holds_alternative<int>(param.value))
            {
                paramElement->SetAttribute("type", "int");
                paramElement->SetAttribute("value", std::get<int>(param.value));
            }
            else if (std::holds_alternative<long long>(param.value))
            {
                paramElement->SetAttribute("type", "longlong");
		auto v = std::get<long long>(param.value);
		paramElement->SetAttribute("value", static_cast<uint64_t>(v));
                //paramElement->SetAttribute("value", std::get<long long>(param.value));
            }
            else if (std::holds_alternative<double>(param.value))
            {
                paramElement->SetAttribute("type", "double");
                paramElement->SetAttribute("value", formatDouble(std::get<double>(param.value)).c_str());
            }
            else if (std::holds_alternative<bool>(param.value))
            {
                paramElement->SetAttribute("type", "bool");
                paramElement->SetAttribute("value", std::get<bool>(param.value) ? "true" : "false");
            }
            else if (std::holds_alternative<std::string>(param.value))
            {
                paramElement->SetAttribute("type", "string");
                paramElement->SetAttribute("value", std::get<std::string>(param.value).c_str());
            }
			parent->InsertEndChild(paramElement);
        }


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
                    }
                path=dir.string();
               // fname=filename.string();
            }

			tinyxml2::XMLError eResult=doc.LoadFile(fname.c_str());                        
			if (eResult == tinyxml2::XML_SUCCESS)
			{
				rootElement = doc.RootElement();
				readScene();
				// if (calculation_enabled) doCalculations();
			}
			else
				std::cerr << "Could not read XML-File:" << fname << std::endl;
			
		}

        bool xmlReader::readRequest(std::string& request)
        {
         //   tinyxml2::XMLDocument doc;
			calculation_enabled = false;
            if (doc.Parse(request.c_str(), request.size()) != tinyxml2::XML_SUCCESS)
                return false;
			rootElement = doc.RootElement();
			readScene();
            readJobs();
            return true;
        }

        void xmlReader::readScene()
		{
			sceneElement = rootElement->FirstChildElement("Scene");
			if (sceneElement != NULL)
			{
				double dv;
                int iv;

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

                iv = sceneElement->IntAttribute("nCellsPerDir", 1000);
#ifdef WITH_NP
                S.setNumberOfCellsPerDirection(iv);
#endif
                S.nS = readCmplx(sceneElement->FirstChildElement("nS"),1.0);
                S.setNumReflex(sceneElement->IntAttribute("nReflex", 0));
				/* look for the detectors */
				readDetectors();

				/* Looking for light sources */
			    readLightSources();
				

                /* Now, let's look for objects */
				readObjects();

                if (calculation_enabled) doCalculations();
			}		
		}

         bool xmlReader::readParam(tinyxml2::XMLElement* paramEll, calculationParam &p)
        {     
            p.name=paramEll->Attribute("name");
            std::string typeStr=paramEll->Attribute("type");
            if (typeStr.compare("int")==0)
            {
                int val=paramEll->IntAttribute("value",0);
                p.value=val;
            }
            else if (typeStr.compare("longlong")==0)
            {
                long long val=paramEll->Int64Attribute("value",0);
                p.value=val;
            }
            else if (typeStr.compare("double")==0)
            {
                double val=paramEll->DoubleAttribute("value",0.0);
                p.value=val;
            }
            else if (typeStr.compare("bool")==0)
            {
                bool val=paramEll->BoolAttribute("value",false);
                p.value=val;
            }
            else if (typeStr.compare("string")==0)
            {
                std::string val=paramEll->Attribute("value");
                p.value=val;
            }
            else 
				return false;
			return true;
        }

         void xmlReader::readJobs()
         {
             tinyxml2::XMLElement* ell = rootElement->FirstChildElement("Calculations");
             if (ell != NULL)
             {
                 for (auto calcEll = ell->FirstChildElement("Calculation"); calcEll != NULL; calcEll = calcEll->NextSiblingElement("Calculation"))
                 {
                     calculationJob job;
                     job.id = "";
                     job.type = mapString2CalculationToken(calcEll->Attribute("type"));
                     switch (job.type)
                     {
                     case TOKEN_CALCULATION_PULSE:
                     {
                         pulseJobParms parms;
                         parms.trafo.wvl = calcEll->DoubleAttribute("wavelength", 1.0);
                         parms.trafo.nR = calcEll->IntAttribute("numReflex", raytracing::INEL_MAX_NREFLEX);
                         parms.trafo.dt = calcEll->DoubleAttribute("pulseWidth", 100.0);
                         parms.trafo.nI = calcEll->IntAttribute("numSpectralRanges", 20);
						 parms.trafo.nS = calcEll->IntAttribute("numWavelengthsPerRange", 10);
                         parms.trafo.repetitionTime = calcEll->DoubleAttribute("repetitionTime", 1000.0);
						 parms.trafo.spatialResolution = calcEll->DoubleAttribute("spatialResolution", 1.0);
                         parms.trafo.number_of_threads = calcEll->IntAttribute("numThreads", 5);
                         parms.time = calcEll->DoubleAttribute("time", 0.0);
                         parms.offsetTime = calcEll->DoubleAttribute("offsetTime", 0.0);
                         int i = 0;
						 auto refractiveIndexListEll = calcEll->FirstChildElement("RefractiveIndexList");
                         for (auto obj = S.Obj.begin(); obj != S.Obj.end(); ++obj)
                         {
                          //   if ((*obj)->isActive())
                             {
                                 std::string name = "n" + std::to_string(i);         
                                 //std::string funcname = refractiveIndexListEll->Attribute("n0");

								  std::string funcname = refractiveIndexListEll->Attribute(name.c_str());
                                 raytracing::nFnPtr fp = GOAT::raytracing::keyToN.at(funcname);  // fp ist cplx(*)(double)								
								 parms.trafo.nList.push_back(fp);
                             }
                             i++;
                         }
						 job.parms = parms;
                         jobs.push_back(job);
                         break;
                     }

                     
                     }
                 }
             }
         }

         struct detectorLink
         {
             raytracing::Propagator* det;
             std::vector<std::string> linkIDs;
		 };

        void xmlReader::readDetectors()
        {
			std::vector<detectorLink> pendingLinks; // to store the links for Kirchhoff detectors until all detectors are read
			tinyxml2::XMLElement* ell;
			ell = sceneElement->FirstChildElement("Detectors");
            std::vector<std::vector<std::string> > linkList;
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
                    std::string ID = detEll->Attribute("ID");

                    int type = mapString2DetectorToken(typeStr);
                    switch (type)
                    {
                    case TOKEN_DETECTOR_PLANE:
                    {
                        std::vector<std::string> dummy;
                        linkList.push_back(dummy);
                        double d = detEll->DoubleAttribute("d", 1);
                        int n = detEll->IntAttribute("n", 1);
						int n1 = detEll->IntAttribute("n1", -1);
                        if (n1 == -1) n1 = n;
                        int n2 = detEll->IntAttribute("n2", -1);
                        if (n2 == -1) n2 = n;
						double d1 = detEll->DoubleAttribute("d1", -1);
						if (d1 < 0) d1 = d;
						double d2 = detEll->DoubleAttribute("d2", -1);
						if (d2 < 0) d2 = d;
                        Det.push_back(new raytracing::DetectorPlane(Pos, Dir, d1, d2, n1, n2));
                        Det[numDet]->fname = filename;
                        S.addDetector(Det[numDet]);
                        Det[numDet]->load(filename.c_str());
                        Det[numDet]->setID(ID);
                        numDet++;
                       }
                    break;

                    case TOKEN_DETECTOR_KIRCHHOFF:
                    {						
                        double d = detEll->DoubleAttribute("d", 1);
                        int n = detEll->IntAttribute("n", 1);
                        std::vector<raytracing::DetectorPlane *> sources;
                        bool cancel = false;
                        Det.push_back(new raytracing::Kirchhoff(1.0,Pos, Dir, d, n));
						double wvl = detEll->DoubleAttribute("wavelength", 1.0);
                        detectorLink links;
						links.det = (raytracing::Kirchhoff*)Det[numDet];
                        for (tinyxml2::XMLElement* link=detEll->FirstChildElement ("Link"); link != NULL; link = link->NextSiblingElement("Link"))
                        {
                            std::string linkID = link->Attribute("ID");
                            raytracing::Detector* det = S.getDetector(linkID);
                            cancel = det == NULL; // check, if detector with ID exists
							if (!cancel)
                            {
								links.linkIDs.push_back(linkID); 
                            }
                        }
						pendingLinks.push_back(links);
                        Det[numDet]->setID(ID);
                        S.addDetector(Det[numDet]);
                        numDet++;				
                    }
                    break;

                    case TOKEN_DETECTOR_ANGULAR_SPECTRUM:
                        {
                            double d = detEll->DoubleAttribute("d", -1);
                            double d1, d2;
                            if (d == -1)
                            {
                                d1 = detEll->DoubleAttribute("d1", 1);
                                d2 = detEll->DoubleAttribute("d2", 1);
                            }
                            else
                            {
                                d1 = d;
                                d2 = d;
                            }

                            int n = detEll->IntAttribute("n", -1);
                            int n1, n2;
                            if (n == -1)
                            {
                                n1 = detEll->IntAttribute("n1", 1);
                                n2 = detEll->IntAttribute("n2", 1);
                            }
                            else
                            {
                                n1 = n;
                                n2 = n;
                            }
                            double wvl = detEll->DoubleAttribute("wavelength", 1.0);
                            maths::Vector<double> e1, e2;
                            e1 = readVector(detEll->FirstChildElement("e1"));
                            e2 = readVector(detEll->FirstChildElement("e2"));
							Det.push_back(new raytracing::AngularSpectrum(wvl, Pos, e1, e2, n1, n2));
                            bool cancel = false;
                            detectorLink links;
                            for (tinyxml2::XMLElement* link = detEll->FirstChildElement("Link"); link != NULL; link = link->NextSiblingElement("Link"))
                            {
                                std::string linkID = link->Attribute("ID");
                                raytracing::Detector* det = S.getDetector(linkID);
                                cancel = det == NULL; // check, if detector with ID exists
                                if (!cancel)
                                {
                                    links.linkIDs.push_back(linkID);
									links.det =(GOAT::raytracing::AngularSpectrum*) Det[numDet];
                                }
                            }
                            pendingLinks.push_back(links);
                            Det[numDet]->fname = filename;  
                            Det[numDet]->load(filename.c_str());
                            Det[numDet]->setID(ID);
                            S.addDetector(Det[numDet]);
                            numDet++;
                        }
                        break;
                    } // switch(type)
				} // for...
              
                for (auto& links : pendingLinks)
                {
                    for (auto linkID : links.linkIDs)
                        links.det->addDetector((raytracing::DetectorPlane *)(S.getDetector(linkID)));

                }

	   } // if (ell != NULL)

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
				GOAT::raycount_t numRays;
                int numRaysRT;
				double wavelength;
				double size;
                double power;
				GOAT::raytracing::LightSrc* ls = NULL;
				for (tinyxml2::XMLElement* lsEll = ell->FirstChildElement("LightSource"); lsEll != NULL; lsEll = lsEll->NextSiblingElement("LightSource"))
				{
			        	std::string typeStr;
					typeStr = lsEll->Attribute("type");
					Pos = readVector(lsEll->FirstChildElement("Position"));
					 double nrays = lsEll->DoubleAttribute("numRays", 100);
					 if (nrays < 0) nrays = 100;
					 numRays = static_cast<GOAT::raycount_t>(nrays);
                    numRaysRT = lsEll->IntAttribute("numRaysRT", 10);
					wavelength = lsEll->DoubleAttribute("wavelength", 1.0);
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
                                                        double thetaMax = lsEll->DoubleAttribute("thetaMax", M_PI);
                                                        ((GOAT::raytracing::LightSrcPoint_mc *)ls)->setThetamax(thetaMax);
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
                                                           ls->setk(k);
                                                           LS.push_back(ls);
														   break;
                                                        }
                    case TOKEN_LIGHTSOURCE_GAUSSIAN_RING_MC :
                                                        {
                                                           double rmin, rmax;
                                                         double width;
                                                         double FWHM;
                                                         rmin=lsEll->DoubleAttribute("rmin",0.0);
                                                         rmax=lsEll->DoubleAttribute("rmax",100.0);
                                                        // width=lsEll->DoubleAttribute("width",rmax);
                                                         FWHM = lsEll->DoubleAttribute("FWHM", rmax);
                                                           ls=new GOAT::raytracing::LightSrcRingGauss_mc(Pos, numRays, wavelength, rmin, rmax,Pol);
                                                         ((GOAT::raytracing::LightSrcRingGauss_mc *)ls)->setFWHM(FWHM);
														 GOAT::maths::Vector<double> k = readVector(lsEll->FirstChildElement("Direction"),0,0,1);
														 ls->setk(k);
                                                                                     
                                                         LS.push_back(ls);
                                                         break;
                                                        }


					}
                    LS[numLS]->setNumRays(numRays);
                    LS[numLS]->setNumRaysRT(numRaysRT);

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
            std::string ID;
			std::complex<double> n;
            bool isRough;
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
					isRough = objEll->BoolAttribute("isrough", false);
                    auto text=objEll->Attribute("ID");
                    if (text == nullptr) ID = "new_Object";
                    else ID = text;
					// GOAT::raytracing::ObjectShape* obj = NULL;
					n = readCmplx(objEll->FirstChildElement("n"), 1.0);
					int type = mapString2ObjectToken(typeStr);
					raytracing::ObjectShape* obj = nullptr;
					switch (type)
					{ 

					case TOKEN_OBJECT_ELLIPSOID: {
													GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10.0, 10.0, 10.0);
                                                    obj = new GOAT::raytracing::Ellipsoid(Pos, Dimensions, n);                                                    
													obj->setMatrix(alpha, beta, gamma);
                                                    obj->setActive(isActive);
												 }
											   break;
					case TOKEN_OBJECT_BOX: {
													GOAT::maths::Vector<double> Dimensions = readVector(objEll->FirstChildElement("Dimension"), 10, 10, 10);
													obj = new GOAT::raytracing::Box(Pos, Dimensions, n);
													obj->setMatrix(alpha, beta, gamma);
													obj->setActive(isActive);													
										    }
										 break;
					case TOKEN_OBJECT_SURFACE: 
										   {
											obj=new GOAT::raytracing::surface(Pos, n);
											fileTypeStr = objEll->Attribute("filetype");

											if (fileTypeStr.compare(".srf") == 0)
											{
												fileName = objEll->Attribute("filename");
                                          
                                                if (path.size()>0)
                                                {
                                                    std::string sep = "/";
                                                    std::filesystem::path p(fileName);
                                                    if (p.is_relative())
                                                    fileName = path + sep + fileName;
                                                }

                                                

												if (!fileName.empty()) ((GOAT::raytracing::surface*)obj)->createsurface(fileName);
											}

											if (fileTypeStr.compare(".stl") == 0)
											{
												fileName = objEll->Attribute("filename");
                                                if (path.size() > 0)
                                                {
                                                    std::filesystem::path p(fileName);
                                                    if (p.is_relative())
                                                    {
                                                        std::string sep = "/";
                                                        fileName = path + sep + fileName;
                                                    }
                                                }
												if (!fileName.empty())
												((GOAT::raytracing::surface*)obj)->importBinSTL(fileName);
											}
										   }
										   obj->setMatrix(alpha, beta, gamma);
										   obj->setActive(isActive);
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


												obj=new GOAT::raytracing::sphericLens(Pos,n,lensparms);
												obj->setMatrix(alpha, beta, gamma);
												obj->setActive(isActive);
                                                break;						
											}
                    case TOKEN_OBJECT_CONE: 
                                            {
                                                double height, radius;
                                                height=objEll->DoubleAttribute("height",100);
                                                radius=objEll->DoubleAttribute("radius",100);
                                                obj=new GOAT::raytracing::Cone(Pos,radius,height,n);
                                                obj->setMatrix(alpha, beta, gamma);
												obj->setActive(isActive);
												break;
                                            }

                    case TOKEN_OBJECT_CYLINDER:
                                            {                                                
                                                double height, radius;
                                                height = objEll->DoubleAttribute("height", 1);
                                                radius = objEll->DoubleAttribute("radius", 1);                                                
                                                obj=new GOAT::raytracing::Cylinder(Pos, radius, height, n);
                                                obj->setMatrix(alpha, beta, gamma);
                                                obj->setActive(isActive);
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
                                                obj=new GOAT::raytracing::VortexPlate(Pos, radius, height, dh, m, n);
                                                obj->setMatrix(alpha, beta, gamma);
                                                obj->setActive(isActive);
                                                break;
                            
                                            }

					}
                    if (obj != nullptr)
                    {
                        if (isRough)
                        {
							double sigma = objEll->DoubleAttribute("sigma", 0.0);
                            obj = makeRough(obj,sigma);
                        }
                        Obj.push_back(obj);
                        S.addObject(obj);
                        double sf = objEll->DoubleAttribute("scaling", 1);
                        if ((sf != 1) && (sf > 0)) Obj[numObj]->scale(sf);
                        Obj[numObj]->nFunc() = GOAT::raytracing::n_Vacuum;
                        Obj[numObj]->setPos(Pos);

                        std::string objID = "object_" + std::to_string(numObj);
                        Obj[numObj]->setID(objID);
                        numObj++;
                    }
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
						std::cout << "Calculation type: " << typeStr << std::endl;
                        // change number of rays, if given
                        int numRays;
                        std::vector<int> numRays_old;
                        numRays = objEll->IntAttribute("numRays", 0);
                        bool numRaysChanged=false;
                        if (numRays > 0)
									{
										// store the old values 									 
                                        numRaysChanged=true;
										for (int i = 0; i < S.getNumberOfLightSources(); i++)
										{
											numRays_old.push_back(S.LS[i]->getNumRays());
											S.LS[i]->setNumRays(numRays);
										}
									}

                        int numReflex;
                        numReflex = objEll->IntAttribute("numReflex", 0);
						int numThreads = objEll->IntAttribute("numThreads", 1);
                        S.setNumThreads(numThreads);
						if (!typeStr.empty())
						{							
							type = mapString2CalculationToken(typeStr);							
                            switch (type)
							{
                            case TOKEN_CALCULATION_PURE:
                            {
								std::cout << "do pure raytracing calculation" << std::endl;
                             GOAT::raytracing::Raytrace_pure rt(S);      
                             rt.setNumReflex(S.getNumReflex());                       
                             rt.trace();
							 std::cout << S.getNumberOfDetectors() << " detectors in scene" << std::endl;
                             for (auto det : S.Det)
                             {
                                 std::cout << "Detector: " << det->getID() << " with type" << det->Type() << std::endl;
                                 std::cout << "total intensity(rt): " << std::endl;
                                 std::cout << "total intensity: " << det->getTotalIntensity() << std::endl;
                                 int type = det->Type();
                                 if (type == TOKEN_DETECTOR_KIRCHHOFF || type == raytracing::DETECTOR_KIRCHHOFF)
                                 {
                                     GOAT::raytracing::Kirchhoff* K = (GOAT::raytracing::Kirchhoff*)det;
                                     std::cout << "Kirchhoff detector: " << K->getID() << "with " << K->numberOfSources() << " sources" << std::endl;
                                     K->setNumberOfThreads(numThreads);
                                     K->calc();
                                 }

                                 if (type == TOKEN_DETECTOR_ANGULAR_SPECTRUM || type == raytracing::DETECTOR_ANGULAR_SPECTRUM)
                                 {
                                     GOAT::raytracing::AngularSpectrum* AS = (GOAT::raytracing::AngularSpectrum*)det;
                                     std::cout << "Angular spectrum detector: " << AS->getID() << std::endl;
                                     AS->setNumberOfThreads(numThreads);
                                     AS->calc();
                                 }
                             }
                                 break;
                             }
							case TOKEN_CALCULATION_PATH:
							{
                                std::cout << "do path calculation" << std::endl;
                                std::string fname = objEll->Attribute("filename");								
								int numDet=S.getNumberOfDetectors();
                                // S.nDet=0;
								if (!fname.empty())
								{
                                    /*
									GOAT::raytracing::Raytrace_Inel rt(S);
                                    GOAT::raytracing::RRTParms rrtparms;
                                    
                                    rt.setExcitationFieldOnly();
                                    rt.setNumReflex(numReflex);
                                    rt.trace(rrtparms);
									*/
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
                                    for (int i=0; (i<S.getNumberOfObjects()) && (!failed); i++)
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

                                          for (int i = 0; i < S.getNumberOfObjects(); i++)
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
                                        for (int i = 0; i < S.getNumberOfObjects(); i++)
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
									/*for (int i = 0; i < S.getNumberOfObjects(); i++)
									{
										if (S.Obj[i]->Active)
										{
											fullfname = fname + std::to_string(i) + ".dat";
											GOAT::raytracing::saveFullE(rt.SGE[0], fullfname, i);
										}
									}*/
								
							}
							} // switch
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
										for (int i = 0; i < S.getNumberOfLightSources(); i++)
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
                for (int i=0; (i<S.getNumberOfObjects()) && (!failed); i++)
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
				      for (int i = 0; i < S.getNumberOfObjects(); i++)
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
                    for (int i = 0; i < S.getNumberOfObjects(); i++)
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
                pc.setNumWavelengthsPerRange(objEll->IntAttribute("numWavelengthsPerRange", trafoparms.nS));
                pc.setPulseWidth(objEll->DoubleAttribute("pulseWidth",trafoparms.dt));
                pc.setSpectralRanges(objEll->IntAttribute("numSpectralRanges", trafoparms.nI));
			    //pc.setReferenceTime(objEll->IntAttribute("Reference_time", pc.getReferenceTime()));
                // pc.setNumberOfThreads(objEll->IntAttribute("NumberOfThreads",pc.getNumberOfThreads()));
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
                    return;
                }

                std::string refStr;
                int refIndexToken;
                for (int i=0; (i<S.getNumberOfObjects()) && (!failed); i++)
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
				
                      for (int i = 0; i < S.getNumberOfObjects(); i++)
                      {
                        if (S.Obj[i]->isActive())
                        {
                            fullfname = fname + std::to_string(i) + ".dat";
                       //     GOAT::raytracing::saveFullE(pc.trafo.SAres, fullfname, i);
							
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
        xmlWriter::xmlWriter(const GOAT::raytracing::Scene &scene) : S(scene)
        {
            
        }

        void xmlWriter::write(std::string fname)
        {
			tinyxml2::XMLDocument doc;
            buildDOM(doc);
            tinyxml2::XMLError e = doc.SaveFile(fname.c_str());
        }

        std::string xmlWriter::prepareRequest(std::vector<calculationJob>& jobs)
        {
			tinyxml2::XMLDocument doc;
			

			buildDOM(doc);
            auto calculations = doc.NewElement("Calculations");
            auto results = doc.NewElement("Data");
            auto* root = doc.RootElement();
            root->InsertEndChild(calculations);
            root->InsertEndChild(results);
            for (const auto& job : jobs)
            {
                addCalculation2DOM(doc, calculations, job);
                addResult2DOM(doc, results, job);
            }
            
			tinyxml2::XMLPrinter printer;
            doc.Print(&printer);

            std::string request(printer.CStr(), printer.CStrSize() - 1);
            return request;
        }

        void xmlWriter::addCalculation2DOM(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* calculations, calculationJob job)
        {
			auto calculation = doc.NewElement("Calculation");
            calculation->SetAttribute("type", calculationToken[job.type-200].c_str());
            calculation->SetAttribute("numThreads", std::get<pulseJobParms>(job.parms).trafo.number_of_threads);
            calculation->SetAttribute("wavelength", formatDouble(std::get<pulseJobParms>(job.parms).trafo.wvl).c_str());
            switch (job.type)
            {
                case TOKEN_CALCULATION_PULSE:
					
                    calculation->SetAttribute("numReflex", std::get<pulseJobParms>(job.parms).trafo.nR);
                    calculation->SetAttribute("numWavelengthsPerRange", std::get<pulseJobParms>(job.parms).trafo.nS);
                    calculation->SetAttribute("pulseWidth", formatDouble(std::get<pulseJobParms>(job.parms).trafo.dt).c_str());
                    calculation->SetAttribute("numSpectralRanges", std::get<pulseJobParms>(job.parms).trafo.nI);
				     calculation->SetAttribute("spatialResolution", formatDouble(std::get<pulseJobParms>(job.parms).trafo.spatialResolution).c_str());
                    calculation->SetAttribute("repetitionTime", formatDouble(std::get<pulseJobParms>(job.parms).trafo.repetitionTime).c_str());
					calculation->SetAttribute("numLoops", std::get<pulseJobParms>(job.parms).numLoops);
                    int i = 0;
					auto refractiveIndexList = doc.NewElement("RefractiveIndexList");
                    for (auto obj=S.Obj.begin(); obj!=S.Obj.end(); ++obj)
                    {
                        if ((*obj)->isActive())
                        {
                            auto p = (*obj)->nFunc().target<raytracing::nFnPtr>();
                            std::string entry= GOAT::raytracing::nToKey.at(*p);
							std::string nStr = "n" + std::to_string(i);
                            refractiveIndexList->SetAttribute(nStr.c_str(), entry.c_str());
                        }
                        i++;
					}
					calculation->InsertEndChild(refractiveIndexList);
                    calculation->SetAttribute("time", std::get<pulseJobParms>(job.parms).time);
                    calculation->SetAttribute("offsetTime", std::get<pulseJobParms>(job.parms).offsetTime);
					break;

            }
            calculations->InsertEndChild(calculation);
		}

        void xmlWriter::addResult2DOM(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* results, calculationJob job)
        {
            auto result = doc.NewElement("HDF5");
                result->SetAttribute("filename","results.h5");
                results->InsertEndChild(result);
        }





        void xmlWriter::buildDOM (tinyxml2::XMLDocument &doc)
        {
       
            tinyxml2::XMLElement* root; ///< root XML Element
            tinyxml2::XMLElement* scene; ///< XML Element to the Scene section
            tinyxml2::XMLElement* lightSrcs; ///< XML Element to the LightSources section
            tinyxml2::XMLElement* objects; ///< XML Element to the Objects section
            tinyxml2::XMLElement* detectors; ///< XML Element to the Detectors section 
            tinyxml2::XMLElement* dataEntries;
            tinyxml2::XMLDeclaration* decl = doc.NewDeclaration(R"(xml version="1.0" encoding="utf-8")");
            doc.InsertFirstChild(decl);
          root=doc.NewElement("Root");
          doc.InsertEndChild(root);
          scene=doc.NewElement("Scene");
          scene->SetAttribute("r0", formatDouble(S.r0).c_str());
          scene->SetAttribute("nCellsPerDir", static_cast<int64_t> (S.getNumberOfCellsPerDirection()));
		  scene->SetAttribute("nReflex", S.getNumReflex());
		  scene->InsertEndChild(addComplex2DOM(doc, "nS", S.nS));
          root->InsertEndChild(scene);
          if (S.getNumberOfLightSources() > 0)
          {
              tinyxml2::XMLElement* lightSrc;
              lightSrcs = doc.NewElement("LightSources");
               for (int i = 0; i < S.getNumberOfLightSources(); i++)
                   addLightSrc2DOM(doc, lightSrcs, i);              
               scene->InsertEndChild(lightSrcs);
          }

          if (S.getNumberOfObjects() > 0)
          {
               objects=doc.NewElement("Objects");
               for (int i=0; i<S.getNumberOfObjects(); i++)
                    addObject2DOM(doc, objects, i);
                scene->InsertEndChild(objects);
          }

		  
          if (S.getNumberOfDetectors() > 0)
          {
            detectors=doc.NewElement("Detectors");
            for (int i = 0; i < S.getNumberOfDetectors(); i++)
            {
                addDetector2DOM(doc, detectors, i);
             }
            scene->InsertEndChild(detectors);
          }
          


        }

        void xmlWriter::addLightSrc2DOM(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* lightSrcs, int i)
        {
            
            auto lightSrc = doc.NewElement("LightSource");
            int type = S.LS[i]->type;
            int typeh = type < 10 ? type-1 : type - 5;
            lightSrc->SetAttribute("type", LSTYPES[typeh].c_str());
            lightSrc->SetAttribute("numRays", S.LS[i]->getNumRays());
			lightSrc->SetAttribute("numRays", std::to_string(S.LS[i]->getNumRays()).c_str());                     
            lightSrc->SetAttribute("numRaysRT", S.LS[i]->getNumRaysRT());
            lightSrc->SetAttribute("wavelength", formatDouble(S.LS[i]->getWavelength()).c_str());

            lightSrc->InsertEndChild(addVectorD2DOM(doc, "Position", S.LS[i]->Pos));
            lightSrc->InsertEndChild(addVectorC2DOM(doc, "Polarisation", S.LS[i]->Pol));
            
            switch (type)
            {
              case raytracing::LIGHTSRC_SRCTYPE_PLANE:
              case raytracing::LIGHTSRC_SRCTYPE_PLANE_MC:
                  {
                  raytracing::LightSrcPlane* ls = (raytracing::LightSrcPlane*)S.LS[i];
                   lightSrc->InsertEndChild(addVectorD2DOM(doc, "Direction", ls->getk()));
                   lightSrc->SetAttribute("size", formatDouble(ls->D).c_str());                   
                  }
                  break;

              case raytracing::LIGHTSRC_SRCTYPE_GAUSS:
              {
                  raytracing::LightSrcGauss* ls = (raytracing::LightSrcGauss*)S.LS[i];
                  lightSrc->InsertEndChild(addVectorD2DOM(doc, "FocusPosition", ls->getFocuspos()));
              }
              case raytracing::LIGHTSRC_SRCTYPE_GAUSS_MC:
              {
                  raytracing::LightSrcGauss* ls = (raytracing::LightSrcGauss*)S.LS[i];
                  lightSrc->SetAttribute("size", formatDouble(ls->D).c_str());
                  lightSrc->SetAttribute("w0", formatDouble(ls->w0).c_str());
              }
              break;

              case raytracing::LIGHTSRC_SRCTYPE_POINT:
              case raytracing::LIGHTSRC_SRCTYPE_POINT_MC:
              {
                  auto lsp = (raytracing::LightSrcPoint_mc*)S.LS[i];
                  lightSrc->SetAttribute("thetaMax", formatDouble(lsp->getThetamax()).c_str());
              }

              case raytracing::LIGHTSRC_SRCTYPE_RING_GAUSS_MC:
              {
                  auto lsg = (raytracing::LightSrcRingGauss_mc*)S.LS[i];
                  lightSrc->SetAttribute("FWHM", formatDouble(lsg->getFWHM()).c_str());
              }
              case raytracing::LIGHTSRC_SRCTYPE_RING :
              case raytracing::LIGHTSRC_SRCTYPE_RING_MC:
              {
                  raytracing::LightSrcRing* ls = (raytracing::LightSrcRing*)S.LS[i];
                  lightSrc->SetAttribute("rmin", formatDouble(ls->getRmin()).c_str());
                  lightSrc->SetAttribute("rmax", formatDouble(ls->getRmax()).c_str());
                  lightSrc->InsertEndChild(addVectorD2DOM(doc, "Direction", ls->getk()));
              }
            }
            lightSrcs->InsertEndChild(lightSrc);
        }

        void xmlWriter::addObject2DOM(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* objects, int i)
        {
            auto object = doc.NewElement("Object");
            int typeh = S.Obj[i]->Type()-10000;
            int type = S.Obj[i]->Type();


            // ---------------- global parameters ----------------
            object->SetAttribute("type",objectToken[typeh].c_str());                        
            object->InsertEndChild(addVectorD2DOM(doc,"Position", S.Obj[i]->pos()));
            object->SetAttribute("alpha",formatDouble(S.Obj[i]->getAlpha()/M_PI*180.0).c_str());
            object->SetAttribute("beta",formatDouble(S.Obj[i]->getBeta()/M_PI*180.0).c_str());
            object->SetAttribute("gamma",formatDouble(S.Obj[i]->getGamma()/M_PI*180.0).c_str());
            object->SetAttribute("isactive",S.Obj[i]->isActive());
            object->InsertEndChild(addComplex2DOM(doc, "n",S.Obj[i]->getn()));            
            object->SetAttribute("scaling",formatDouble(S.Obj[i]->getScale()).c_str());
			object->SetAttribute("ID", S.Obj[i]->getID().c_str());

            // --------------- special parameters ----------------
                     
            if (S.Obj[i]->isRough())
            {
                object->SetAttribute("isrough", true);
                raytracing::roughInterface* robj =dynamic_cast<raytracing::roughInterface *>(S.Obj[i]);
                if (robj) object->SetAttribute("sigma",formatDouble(robj->getSigma()).c_str());       
            }
            switch (type) 
            {
                case OBJECTSHAPE_ELLIPSOID : 
                    {
                    auto obj=(raytracing::Ellipsoid *) S.Obj[i];
                        object->InsertEndChild(addVectorD2DOM(doc, "Dimension",obj->r));
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
                case OBJECTSHAPE_CONE : 
                {
                    auto obj = (raytracing::Cone*)S.Obj[i];
                    object->SetAttribute("height", obj->getHeight());
                    object->SetAttribute("radius", obj->getRadius());
                }
                    break;
                case OBJECTSHAPE_ASPHERIC_LENS : break;
                case OBJECTSHAPE_SPHERIC_LENS : 
                    {
                        auto obj=(raytracing::sphericLens *) S.Obj[i];
                        raytracing::lensParms lensparms = obj->getParms();        
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
                        object->InsertEndChild(addVectorD2DOM(doc, "Dimension",obj->d));
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
			if (S.Obj[i]->isRough()) addRoughness2DOM(doc, object, S.Obj[i]);
            objects->InsertEndChild(object);
            
        }


        void xmlWriter::addRoughness2DOM(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* object, raytracing::ObjectShape *obj)
        {
			object->SetAttribute("isrough", true);
            raytracing::roughInterface* robj = obj->getRoughObj();
            auto rough = doc.NewElement("Roughness");
			switch (robj->getScatteringType())
			{
                case raytracing::ScatteringType::Gaussian:
                        rough->SetAttribute("type", "gaussian");
                        rough->SetAttribute("sigma", formatDouble(robj->getSigma()).c_str());
                        break;
                
                case raytracing::ScatteringType::CosineCone:
                        rough->SetAttribute("type", "cosine_cone");
                        break;

                case raytracing::ScatteringType::UniformCone :
                    rough->SetAttribute("type", "uniform_cone");
                    break;

				    
			}
            object->InsertEndChild(rough);
			
			
        }

        void xmlWriter::addDetector2DOM(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* detectors, int i)
        {
            auto detector = doc.NewElement("Detector");

            auto offCaller = (ptrdiff_t)((char*)&S.Det - (char*)&S);
            auto nGetter = S.getNumberOfDetectors();
            auto nDirect = S.Det.size();

         
			if (i >= S.getNumberOfDetectors())
            {
                std::cerr << "Error in addDetector2DOM: No such detector (i=" << i << ")! Skipped writing this detector!" << std::endl;
                return;
            }
            else
            {
                if (S.getNumberOfDetectors() > 0)
                {
                    int type = S.Det[i]->Type();
                    int typeh = type - 20000;
                    detector->SetAttribute("type", detectorToken[typeh].c_str());
                    detector->InsertEndChild(addVectorD2DOM(doc, "Position", S.Det[i]->position()));
                    detector->InsertEndChild(addVectorD2DOM(doc, "Direction", S.Det[i]->norm()));
                    detector->SetAttribute("filename", S.Det[i]->fname.c_str());
                    detector->SetAttribute("ID", S.Det[i]->getID().c_str());
                    S.Det[i]->save(S.Det[i]->fname.c_str());
                    switch (type)
                    {
                    case raytracing::DETECTOR_PLANE:
                    {
                        auto det = (raytracing::DetectorPlane*)S.Det[i];
                        detector->SetAttribute("d1", formatDouble(det->D1()).c_str());  
                        detector->SetAttribute("d2", formatDouble(det->D2()).c_str());
                        detector->SetAttribute("n1", det->N1());  
						detector->SetAttribute("n2", det->N2());
                    }
                    break;



#ifdef WITH_OPENMP
                    case raytracing::DETECTOR_KIRCHHOFF:
                    {
                        auto det = (raytracing::Kirchhoff*)S.Det[i];
                        detector->SetAttribute("d", formatDouble(det->D1()).c_str()); // we assume, that d1=d2 
                        detector->SetAttribute("n", det->N1()); // we also assume that n1=n2 
                        auto sources = det->getSources();
                        for (auto src : sources)
                        {
                            auto srcEll = doc.NewElement("Link");
                            detector->InsertEndChild(srcEll);
							srcEll->SetAttribute("ID", src->getID().c_str());
                        }
                    }
                    break;

                    case raytracing::DETECTOR_ANGULAR_SPECTRUM:
                        {
                            auto det = (raytracing::AngularSpectrum*)S.Det[i];
                            detector->SetAttribute("d1", formatDouble(det->D1()).c_str());
                            detector->SetAttribute("d2", formatDouble(det->D2()).c_str());
                            detector->SetAttribute("n1", det->N1());
                            detector->SetAttribute("n2", det->N2());
                            detector->InsertEndChild(addVectorD2DOM(doc, "e1", det->gete1() / abs(det->gete1()) * det->D1())); 
                            detector->InsertEndChild(addVectorD2DOM(doc, "e2", det->gete2() / abs(det->gete2()) * det->D2()));
							auto sources = det->getSources();
                            for (auto src : sources)
                            {
                                auto srcEll = doc.NewElement("Link");
                                detector->InsertEndChild(srcEll);
                                srcEll->SetAttribute("ID", src->getID().c_str());
                            }
                        }
                        break;
#endif
                    }
                detectors->InsertEndChild(detector);
                }
            }
        }

        

        tinyxml2::XMLElement* xmlWriter::addVectorD2DOM(tinyxml2::XMLDocument& doc, std::string name, maths::Vector<double> v)
        {
            auto vell = doc.NewElement(name.c_str());
            vell->SetAttribute("x", formatDouble(v[0]).c_str());
            vell->SetAttribute("y", formatDouble(v[1]).c_str());
            vell->SetAttribute("z", formatDouble(v[2]).c_str());
            return vell;
        }

        tinyxml2::XMLElement* xmlWriter::addVectorC2DOM(tinyxml2::XMLDocument& doc, std::string name, maths::Vector<std::complex<double> > v)
        {
            auto vell = doc.NewElement(name.c_str());
            vell->InsertEndChild(addComplex2DOM(doc,"x", v[0]));
            vell->InsertEndChild(addComplex2DOM(doc, "y", v[1]));
            vell->InsertEndChild(addComplex2DOM(doc, "z", v[2]));

            return vell;
        }
        tinyxml2::XMLElement* xmlWriter::addComplex2DOM(tinyxml2::XMLDocument& doc, std::string name, std::complex<double> z)
        {
            auto cell = doc.NewElement(name.c_str());
            cell->SetAttribute("real", formatDouble(real(z)).c_str());
            cell->SetAttribute("imag", formatDouble(imag(z)).c_str());
            return cell;
        }


	}
}
