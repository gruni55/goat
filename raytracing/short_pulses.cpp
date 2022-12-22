#include "short_pulses.h"
namespace GOAT
{
	namespace raytracing
	{
		shortPulse::shortPulse(int numWvlRanges, wavelengthRange* wvlRanges, Scene *S, std::vector<std::function<std::complex<double>(double) > > nList)
		{
			this->numWvlRanges = numWvlRanges;
			this->wvlRanges = wvlRanges;
			this->S = S;
			this->nList = nList;
		}
		 
		void shortPulse::startRaytracing()
		{
		}
	}
}