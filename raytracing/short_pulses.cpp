#include "short_pulses.h"
#include "raytrace_usp.h"
#include "supergrid.h"

namespace GOAT
{
	namespace raytracing
	{
		shortPulse::shortPulse(std::vector<wavelengthRange> wvlRanges, Scene S, std::vector<std::function<std::complex<double>(double) > > nList)
		{
			this->numWvlRanges = numWvlRanges;
			this->wvlRanges = wvlRanges;
			this->S = S;
			this->nList = nList;
			if (S.nLS > 0)
				for (int i = 0; i < S.nLS; i++)
					LSwidth.push_back(standard_pulse_width);
		}

		void shortPulse::setPulseWidth(int i, double pulseWidth)
		{
			LSwidth[i] = pulseWidth;
		}
		 
		void shortPulse::startRaytracing()
		{
			double midWvl;
			double omega;
			double sigma2;

			if (wvlRanges.size() == S.nObj)
			{
				// loop over all wavelength ranges 
				for (wavelengthRange range : wvlRanges)
				{
					// -------  make raytracing for the mid wavelength ------
					midWvl = (range.wvlStart + range.wvlStop) / 2.0;
					omega = 2.0 * M_PI * C_LIGHT_MU / (midWvl*1E-6);

					// set all light sources to this wavelength and set the wavelength dependent power
					for (int i = 0; i < S.nLS; i++)
					{
						S.LS[i]->wvl = midWvl;						
						sigma2 = LSwidth[i] / (2.0 * sqrt(2.0 * log(2.0)));
						sigma2 *= sigma2;
						S.LS[i]->P0 = 1.0 / sigma2 * exp(-omega * omega / sigma2);
					}

					// set the refractive indices for this wavelength
					// first for the objects
					for (int i = 0; i < S.nObj; i++)
						S.LS[i]->Obj[0]->n = nList[i](midWvl);

					S.setnS(nList[S.nObj](midWvl)); // set refractive index of the surrounding medium 

					// Do raytracing
					GOAT::raytracing::Raytrace_usp ri(S, 2500);				    
				}
			}
		}
	}
}