#pragma once
#include "raytrace.h"
#include <functional>
#include <vector>

namespace GOAT
{
	namespace raytracing
	{
		struct wavelengthRange
		{
			double wvlStart; 
			double wvlStop; 
			double dWvl; 
		};


	
		/**
		* @brief Class for the calculation of short pulses 
		* 
		* This class provides functions for the calculation of short pulses in the scene S. The short pulse is described by different 
		* wavelength calculations
		*/
		class shortPulse
		{
		  public:
			  shortPulse(int numWvlRanges, wavelengthRange* WvlRanges, Scene *S, std::vector<std::function<std::complex<double>(double) > > nList);

			  wavelengthRange* wvlRanges = 0; ///< Here, the wavelength ranges for the calculation are stored
			  int numWvlRanges = 0; ///< number of different wavelength ranges
			  Scene *S=0; ///< the scene - it contains all information about the light sources, objects etc.
			  std::vector<std::function<std::complex<double>(double) > > nList; ///< List of the refractive index function for each object
		};
	}
}
