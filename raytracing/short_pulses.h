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
		* This structure is a container for each step within the raytracing
		*/
		struct rayListEntry
		{			
			int refIndex; ///< index within the refractive index list (for further details refer to the description of class shortPulse)
			GOAT::maths::Vector<std::complex<double> > E; ///< electric field
		};

		/*
		*/
		struct indexList
		{
			int startIndex;
			int stopIndex;
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
			  void startRaytracing();

			  wavelengthRange* wvlRanges = 0; ///< Here, the wavelength ranges for the calculation are stored
			  int numWvlRanges = 0; ///< number of different wavelength ranges
			  Scene *S=0; ///< the scene - it contains all information about the light sources, objects etc.
			  std::vector<std::function<std::complex<double>(double) > > nList; ///< List of the refractive index function for each object
			  std::vector<rayListEntry> rayList; ///< List of all ray steps			  
		};
	}
}
