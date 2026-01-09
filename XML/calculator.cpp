#include "calculator.h"
namespace GOAT
{
	namespace XML
	{ 
		Calculator::Calculator(raytracing::Scene &S, calculationJob &job)
		{
			this->job = job;
			this->S = S;
		}

		void Calculator::exec()
		{
			switch (job.type)
			{
			case TOKEN_CALCULATION_PULSE: pulseCalculation(); break;
 			}
		}

		void Calculator::pulseCalculation()
		{
			std::cout << "Starting pulse calculation..." << std::endl;
               // int numLoops = objEll->IntAttribute("numLoops", -1);
                GOAT::raytracing::pulseCalculation_rt pc(S);
                // GOAT::raytracing::TrafoParms trafoparms;
				// trafoparms = std::get<raytracing::TrafoParms> job.parms;
                auto& parms = std::get<pulseJobParms>(job.parms);
                pc.setCenterWavelength(parms.trafo.wvl);
                pc.setNumReflex(parms.trafo.nR);
                pc.setPulseWidth(parms.trafo.dt);
                pc.setSpectralRanges(parms.trafo.nI);
                pc.setNumWavelengthsPerRange(parms.trafo.nS);
                
                double repRate = parms.trafo.repetitionTime;
                if (repRate > 0) pc.setRepetitionRate(repRate);
                pc.setSpatialResolution(parms.trafo.spatialResolution);
                
                // double D = objEll->DoubleAttribute("D", -1.0);
                char cs[3];
                std::vector< std::function< std::complex< double >(double) > > nList;
                std::string refFuncName;
                bool failed = false;

                
                pc.setRefractiveIndexFunctions(parms.nFunc);

                double time = parms.time;
                std::cout << "time:" << time << std::endl;
                if (time < 0)
                {
                    double offset = parms.offsetTime;
                    // int objEstimate = objEll->IntAttribute("estimateTimeForObject", 0);
                    //std::cout << "estimated time: " << time << std::endl << std::flush;
                    time += offset;
                }


                std::string fullfname;
                double d;
                // if (D>0)
                {
                    const char* hStr;
                    std::string corrFilename;
                    std::ofstream corrOS;
                    /*hStr = objEll->Attribute("correlationFilename");
                    if (hStr != NULL)
                    {
                        corrOS.open(hStr);
                    }
                    */
                    int loopno = 0;
                    bool cancel = false;
                    do
                    {
                        pc.field(time);

                       /* for (int i = 0; i < S.nObj; i++)
                        {
                            if (S.Obj[i]->isActive())
                            {
                      //          fullfname = fname + std::to_string(i) + ".dat";
                
                                GOAT::raytracing::saveFullE(pc.rt.SA[0], fullfname, i);

                                d = sumabs2(pc.rt.SA[0], i);
                                std::cout << "d=" << d << std::endl;
                            }
                        }

                        
                        if (hStr != NULL) corrOS << d << std::endl;*/
                        loopno++;
                        cancel = (loopno >= parms.numLoops) && (parms.numLoops >= 0);
                        std::cout << "loopno=" << loopno << std::endl;
                    } while (!cancel); // while ( (d>D) || (loopno<2));
                    if (hStr != NULL) corrOS.close();
                }
           
		}
		
	}
}