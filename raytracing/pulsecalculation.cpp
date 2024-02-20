#include "pulsecalculation.h"
#include "fft.h"


namespace GOAT
{
	namespace raytracing
	{
		pulseCalculation::pulseCalculation(Scene S)
		{
			this->S = S;
			setDefaults();
			
			trafo = Trafo(trafoparms);			
		}

		double pulseCalculation::findHitTime(int ObjNo)
		{			
			int nI = trafoparms.nI;
			int nS = trafoparms.nS;
			setSpectralRanges(1);
			setNumWavelengthsPerRange(1);

			field(0);

			// find the first element which was hit by a ray
			bool found = false;
			INDEX_TYPE fix, fiy, fiz;
			for (INDEX_TYPE ix=0; (ix<=rt.SA[0].n[ObjNo][0]) && (!found); ix++)
				for (INDEX_TYPE iy = 0; (iy <= rt.SA[0].n[ObjNo][1]) && (!found); iy++)
					for (INDEX_TYPE iz = 0; (iz <= rt.SA[0].n[ObjNo][2]) && (!found); iz++)
					{
						found = !rt.SA[0].G[ObjNo][ix][iy][iz].empty();
						if (found) { fix = ix; fiy = iy; fiz = iz; }
					}
			
			
			
			double time = 0.0;
			if (found)
			for (auto se=rt.SA[0].G[ObjNo][fix][fiy][fiz][0].step.begin(); se!=rt.SA[0].G[ObjNo][fix][fiy][fiz][0].step.end(); se++)
			{
				time += se->l * real(trafo.nList[se->matIndex](trafoparms.wvl)) / C_LIGHT_MU_FS;
			}
			setSpectralRanges(nI);
			setNumWavelengthsPerRange(nS);
			return time;
		}

		void pulseCalculation::setReferenceTime(double tref)
		{
			this->tref = tref;
			trafo.setReferenceTime(tref);
		}
      
                void pulseCalculation::calcTrafoParms()
                {
 			// double Sigma= (2.0 * sqrt(M_LN2)) / trafoparms.dt; // Spectral sigma
					double Sigma = sqrt(2.0 * M_LN2) / trafoparms.dt;
					double Domega = 8.0 * Sigma;					
			        trafoparms.omegaStart = trafoparms.omega0 -  Domega/2.0;
			        trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
	                double lambdaStart=2.0*M_PI*C_LIGHT_MU_FS / trafoparms.omegaEnd; 
            		double lambdaEnd=2.0*M_PI*C_LIGHT_MU_FS / trafoparms.omegaStart; 
                    std::cout << "% Xwvl-range:" << lambdaStart << "\t" << lambdaEnd << std::endl;		
			trafo.setTrafoParms(trafoparms);
                }
				
		void pulseCalculation::setCenterWavelength(double wvl)
		{
			trafoparms.wvl = wvl;
			trafoparms.omega0 = C_LIGHT_MU_FS / wvl * 2.0 * M_PI;
                        calcTrafoParms();
		}

		void pulseCalculation::setBandwidth(double dWvl)
		{
			this->dWvl = dWvl;
			double Domega = 2.0 * M_PI * C_LIGHT_MU_FS * dWvl / (trafoparms.wvl * trafoparms.wvl);
			trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
			trafoparms.omegaStart = trafoparms.omega0 + Domega / 2.0;
		}


		void pulseCalculation::fieldCalculation()
		{			
			double Sigma = 2.3548 / trafoparms.dt;
			double Domega = 20.0 * Sigma;			
			std::cout << "% nI=" << trafoparms.nI << std::endl;
			double domega = Domega / (double)trafoparms.nI;
			double omega;
			double wvl;
			trafoparms.omegaStart = trafoparms.omega0 - Domega / 2.0;
			trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
			trafo.setTrafoParms(trafoparms);
			for (int i = 0; i < trafoparms.nI; i++)
			{
				omega = trafoparms.omegaStart + i * domega - domega / 2.0;
				wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega;
				for (int ls = 0; ls < S.nLS; ls++)
					S.LS[ls]->wvl = wvl;
				for (int lobj = 0; lobj < S.nObj; lobj++)
					S.Obj[lobj]->setn(trafoparms.nList[lobj](wvl));
				S.setnS(trafoparms.nList[S.nObj](wvl));
				
				rt.trace();				
//				save(rt.SA[1], "test.log");
	//			SA.push_back(rt.SA);
			}			
		}

		void pulseCalculation::fieldCalculation(double omega)
		{
			double wavelength = 2.0 * M_PI / omega * C_LIGHT_MU_FS;

			// clear the old raytracing results
		    // rt.clear();
			rt.setNumReflex(numReflex);

			// now, set the new refractive indices 
			for (int iObj = 0; iObj < S.nObj; iObj++)
				S.Obj[iObj]->setn(trafoparms.nList[iObj](wavelength));
			S.setnS(trafoparms.nList[S.nObj](wavelength));

			// set the wavelength for all light sources
			for (int ls = 0; ls < S.nLS; ls++)
				S.LS[ls]->wvl = wavelength;

			// do the raytracing
			rt.trace();
		}


		void pulseCalculation::setRepetitionRate(double rep)
		{
			double Domega = trafoparms.omegaEnd - trafoparms.omegaStart;
			trafoparms.nS = ceil(Domega / (rep * (double)trafoparms.nI));			
		}

		void pulseCalculation::field(double t)
		{
			double omega0 = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.wvl;
			double Domega = 8.0 * 4.0 * M_LN2 / trafoparms.dt;
std::cout << "Domega=" << Domega << std::endl;


//			double Domega = 8 * M_PI * C_LIGHT_MU_FS * dWvl / (4.0 * trafoparms.wvl * trafoparms.wvl - dWvl * dWvl);
//                        double Domega = 2.0 * M_PI * C_LIGHT_MU_FS dWvl / (trafoparms.wvl * trafoparms.wvl);
			double domega = Domega / (double)trafoparms.nI;
			double omegaStart = omega0 - Domega/2.0;
			double omega;
            double wvl1, wvl2;
			rt = Raytrace_usp(S, nn);
			double wvl;
		trafo.initResult(S.r0,rt.SA[0].nges[0], rt.SA[0].nges[1], rt.SA[0].nges[2],S.Obj,S.nObj);
		    // loop over the frequency ranges
			for (int iOmega = 0; iOmega < trafoparms.nI; iOmega++)
			{
				omega = omegaStart + (double)iOmega * domega;				
				wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega; // center wavelength of the current range

                // ------ for output only ------
                wvl1=  2.0 * M_PI * C_LIGHT_MU_FS / (omega-0.5*domega); 
                wvl2=  2.0 * M_PI * C_LIGHT_MU_FS / (omega+0.5*domega);
				std::cout << "%  " << iOmega << ":start FFT (" << wvl << "µm)" << "\twvl1=" << wvl1 << "\twvl2=" << wvl2 << "\tomega=" << omega << std::endl << std::flush;

				fieldCalculation(omega); // do the raytracing				
				trafo.calc(rt.SA, omega - domega * 0.5, omega + domega * 0.5, t); // do the Fourier transform
			}
		}

	/*	void pulseCalculation::field(double t)
		{
			
			if (trafoparms.nList.size() == S.nObj + 1) // process calculation only, if all necessary refractive index functions are given
			{
				if (!raytracingDone)
				{
					std::cout << "Start raytracing...";
					fieldCalculation(); // raytracing is necessary only once
					std::cout << "done." << std:: endl;
				}				
				trafo.calc(SA,t);
				raytracingDone = true;
			}
		}
		*/
		void pulseCalculation::reset()
		{
			// first, let's clear the array where the ray paths and the materials are stored
			for ( SuperArray<std::vector<gridEntry> >  SAElement : SA)
			{
				SAElement.clear();
			}
			SA.clear();
			raytracingDone = false;
		}


		void pulseCalculation::setPulseWidth(double dt)
		{
			trafoparms.dt = dt;
            calcTrafoParms();	 
		}

		void pulseCalculation::setDefaults()
		{
			trafoparms.dt = 100;
			trafoparms.wvl = 1.0;
			trafoparms.nI = 10;
			trafoparms.nR = 1;
			trafoparms.nS = 50;			
			setPulseWidth(trafoparms.dt);			
			setSpatialResolution(1.0);
		}

		void pulseCalculation::setSpectralRanges(int nI)
		{
			trafoparms.nI = nI;
			trafo.setTrafoParms(trafoparms);
		}

		void pulseCalculation::setNumWavelengthsPerRange(int nS)
		{
			trafoparms.nS = nS;
			trafo.setTrafoParms(trafoparms);
		}


		void pulseCalculation::setSpatialResolution(double dx)
		{
			nn = 2.0 * S.r0 / dx;
		}

		void pulseCalculation::setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList)
		{
			trafoparms.nList = nList;	
			trafo.setRefractiveIndexFunctions(nList);			
		}

		void pulseCalculation::setNumReflex(int numReflex)
		{
			this->numReflex = numReflex;
			rt.setNumReflex(numReflex);
		}
	}
}
