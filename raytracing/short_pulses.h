#pragma once
#include "raytrace.h"
#include <functional>
#include <vector>
#include "refractive_index_functions.h"
#include "superarray.h"

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

		constexpr int log2n=10; ///< 2^log2n support points for FFT
		constexpr double standard_pulse_width = 1E-9; ///< standard pulse width in seconds
	
		/**
		* @brief Class for the calculation of short pulses 
		* 
		* This class provides functions for the calculation of short pulses in the scene S. The short pulse is described by different 
		* wavelength calculations. Notice: Wavelengths are given here in um!
		*/
		class shortPulse
		{
		  public:
			  /**
			  * @brief Constructor
			  * @param WvlRanges: Defines the wavelength ranges for the calculation
			  * @param S: Scene, which holds all information about the light sources, objects etc.
			  * @param nList: holds the functions which describe the wavelength dependent refractive index (notice: nList[0...S.NObj-1] -> refractive index of the objects, nList[S.Nobj]: refractive index of the surrounding medium  
			  *                  
			  */
			  shortPulse(std::vector<wavelengthRange> WvlRanges, Scene S, std::vector<std::function<std::complex<double>(double) > > nList);
			  void setPulseWidth(int i, double pulseWidth); ///< set pulse width of the i-th light source
			  void startRaytracing(); ///< performs raytracing process
			  

			  std::vector<wavelengthRange> wvlRanges; ///< Here, the wavelength ranges for the calculation are stored
			  int numWvlRanges = 0; ///< number of different wavelength ranges
			  Scene S; ///< the scene - it contains all information about the light sources, objects etc.
			  std::vector<std::function<std::complex<double>(double) > > nList; ///< List of the refractive index function for each object
			  std::vector<rayListEntry> rayList; ///< List of all ray steps		
			  std::vector<double> LSwidth; ///< pulse widths of the light sources
			  SuperArray< maths::Vector<std::complex<double> > > field;
			  std::vector<unsigned int> bitLT; ///< bit reverse lookup table, used for FFT
		};
	}
}
