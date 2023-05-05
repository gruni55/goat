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
			setPulseWidth(40);
			setSpatialResolution(1.0);
			trafo = Trafo(trafoparms);			
		}

		void pulseCalculation::setReferenceTime(double tref)
		{
			this->tref = tref;
			trafo.setReferenceTime(tref);
		}


		void pulseCalculation::fieldCalculation()
		{			
			double Sigma = 2.3548 / trafoparms.dt;
			double Domega = 2.0 * Sigma;			
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
				rt = Raytrace_usp(S, nn);	
				rt.trace();				
				// save(rt.SA[0], "H:\\data\\data2.log");
				SA.push_back(rt.SA);
			}			
		}

		void pulseCalculation::field(double t)
		{
			
			if (trafoparms.nList.size() == S.nObj + 1) // process calculation only, if all necessary refractive index functions are given
			{
				if (!raytracingDone)
				{
				//	std::cout << "Start raytracing...";
					fieldCalculation(); // raytracing is necessary only once
				//	std::cout << "done." << std:: endl;
			    // save(rt.SA[0], "data.log");
				}				
				trafo.calc(SA,t);
				raytracingDone = true;
			}
		}

		void pulseCalculation::reset()
		{
			// first, let's clear the array where the ray paths and the materials are stored
			for (std::vector < SuperArray<std::vector<gridEntry> > > SAElement : SA)
			{
				for (SuperArray<std::vector<gridEntry> > SAGridElement : SAElement)
					SAGridElement.clear();
				SAElement.clear();
			}
			SA.clear();

			raytracingDone = false;
		}


		void pulseCalculation::setPulseWidth(double dt)
		{
			// dWvl = trafoparms.wvl * trafoparms.wvl * M_LN2 / (M_PI * M_PI * GOAT::raytracing::C_LIGHT_MU_FS * dt);	// FWHM
			trafoparms.dt = dt;
			trafoparms.omega0 = C_LIGHT_MU_FS * 2.0 * M_PI / trafoparms.wvl;			
			double Sigma= (2.0 * sqrt(2.0 * M_LN2)) / dt; // Spectral sigma
			trafoparms.omegaStart = trafoparms.omega0 -  Sigma;
			trafoparms.omegaEnd = trafoparms.omega0 + Sigma;
            double lambdaStart=2.0*M_PI*C_LIGHT_MU_FS / trafoparms.omegaEnd; 
            double lambdaEnd=2.0*M_PI*C_LIGHT_MU_FS / trafoparms.omegaStart; 

                        std::cout << "% wvl-range:" << lambdaStart << "\t" << lambdaEnd << std::endl;		
			trafo.setTrafoParms(trafoparms);
		}

		void pulseCalculation::setDefaults()
		{
			trafoparms.dt = 100;
			trafoparms.wvl = 1.0;
			trafoparms.nI = 1;
			trafoparms.nR = 1;
			trafoparms.nS = 15000;			
			setPulseWidth(trafoparms.dt);
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
	}
}
