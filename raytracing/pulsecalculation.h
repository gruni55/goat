#pragma once

#include "superarray.h"
#include "fft.h"
#include "raytrace_usp.h"
#include <vector>
#include "superarray.h"



namespace GOAT
{
	namespace raytracing
	{
	    /*! \brief This class provides functionality to calculate field distributions for short pulses.  
		*   
		*  The class calculates the field distribution for a short pulsed light source. Also dispersion is considered, therefore a list of 
		*  functionsis required, which describe the wavelength dependence of all objects and the surrounding medium. 
		*  Since short pulses are considered, the light has a spectral width, which depends on the pulse width (Fourier transform). For the
		*  pulse a gaussian shape is assumed. All lengths and the wavelength is given in micro meters. The default wavelength is set to 1.0µm 
		*  and the pulse width 10fs. As spectral width, the full width at half maximum (FWHM) is used. The result is stored in a SuperArray SAres, 
		*  which holds the electric field at a certain time t which was given to the class by calling the function field
		*/
		class pulseCalculation
		{
			public:
				pulseCalculation(Scene S);
				void fieldCalculation(); ///< This function makes the raytracing (normally only used internally)
				void fieldCalculation(double omega); ///< This function makes one raytracing step at frequency omega
				void setPulseWidth(double dt); ///< Sets the spectral width according to the pulse width and adjusts the widht of the subdivisions				
				void setSpatialResolution(double dx); ///< sets the spatial resolution to a value near to dx
				void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList); ///< sets the list of functions, which describe the wavelength dependend refractive index (length must be number of objects + 1)
				void setNumReflex(int numReflex);
				void setSpectralRanges(int n); ///< Number of spectral ranges in which one raytracing calculation is made
				/**
				 * @brief Set number of wavelengths per range.
				 * For the calculation, the spectral range is subdivided into a number of wavelength ranges. This method set the
				 * number of wavelength steps per spectral range.
				 */
				void setNumWavelengthsPerRange(int nS);
				void setCenterWavelength(double wvl); ///< Set center wavelength of the pulse

				void field(double t); ///< This function calculates the fields at time. Keep in mind, that it works only if the class has the list with the refractive index functions
				void reset(); ///< Clears all arrays 		
				void setReferenceTime(double tref);
				Trafo trafo;
				SuperArray<GOAT::maths::Vector<std::complex<double> > > SAres;
					std::vector<SuperArray<std::vector<gridEntry> > >  SA;
					Raytrace_usp rt;

		
				

			private:	
				/* In this function the default values (trafoparms) for the calculations are set as follows:
				* dt  : 1E-14s
				* wvl : 1.0µm
				* nI  : 1
				* 
				*/
				void setDefaults();				
				// std::vector< std::vector<SuperArray<std::vector<gridEntry> > > > SA;

				double dWvl;  // spectral width of the light
				double dRWvl; // spectral width of one subdivision
				int nn;       // number of cells over the whole width of the calculation space (i.e. 2*r0)
				Scene S;
				
				bool raytracingDone = false; ///< If true, the raytracing part was done and the field calculation starts directly				

				TrafoParms trafoparms;
				double tref = 0.0; 
				int numReflex = INEL_MAX_NREFLEX;								
		};
	}
}
