#include "pulsecalculation.h"
#include "fft.h"


namespace GOAT
{
	namespace raytracing
	{
		pulseCalculation::pulseCalculation(Scene S)
		{
			this->S = S;
			setPulseWidth(10E-15);
			setSpatialResolution(1.0);
			trafo = Trafo(trafoparms);			
		}


		void pulseCalculation::fieldCalculation()
		{			
			for (int i = 0; i < trafoparms.nI; i++)
			{
				for (int ls = 0; ls < S.nLS; ls++)
					S.LS[ls]->wvl = trafoparms.lstart + i * dRWvl;
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
					std::cout << "Start raytracing...";
					fieldCalculation(); // raytracing is necessary only once
					std::cout << "done." << std:: endl;
				}
				save(rt.SA[0], "H:\\data\\data.log");
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
			dWvl = trafoparms.wvl * trafoparms.wvl * M_LN2 / (M_PI * M_PI * GOAT::raytracing::C_LIGHT_MU * dt);	// FWHM
			trafoparms.dt = dt;
			trafoparms.lstart = trafoparms.wvl - dWvl / 2.0;
			trafoparms.lstop = trafoparms.wvl + dWvl / 2.0;
			dRWvl = dWvl / trafoparms.nI;
		}

		void pulseCalculation::setDefaults()
		{
			trafoparms.dt = 10E-15;
			trafoparms.wvl = 1.0;
			trafoparms.nI = 4;
			trafoparms.nR = 2;
			trafoparms.nS = 1.0;
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
