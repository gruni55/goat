#pragma once
#include "raytrace.h"
namespace GOAT
{
	namespace raytracing
	typedef struct
	{
		double wvlStart; 
		double wvlStop; 
		double dWvl; 
	} wavelengthRange;


	{
		/**
		* @brief Class for the calculation of short pulses 
		* 
		* This class provides functions for the calculation of short pulses in the scene S. The short pulse is described by different 
		* wavelength calculations
		*/
		class shortPulse
		{
		  public:
			  shortPulse(int numWvlRanges, wavelengthRange* WvlRanges, Scene S);
			  wavelengthRange* WvlRanges = 0;
			  int numWvlRanges = 0;
			  Scene S;
		};
	}
}