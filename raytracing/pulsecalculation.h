#pragma once

#include "superarray.h"
#include "fft.h"
#include "raytrace_usp.h"
#include <vector>




namespace GOAT
{
	namespace raytracing
	{
		#define PULSECALCULATION_CLEAR_SA		  0
		#define PULSECALCULATION_NOT_CLEAR_SA     1  ///< Flag, set when the SA SuperArray should not be cleared
		constexpr int PULSECALCULATION_NOT_CLEAR_RESULT = 2;  ///< Flag, set when the results should not be cleared when calling field () in pulseCalculation

	    /*! \brief This class provides functionality to calculate field distributions for short pulses.  
		*   
        *  This class calculates the field distribution for a short pulsed light source. Also dispersion is considered, therefore a list of
        *  functions is required, which describe the wavelength dependence of all objects and the surrounding medium.
		*  Since short pulses are considered, the light has a spectral width, which depends on the pulse width (Fourier transform). For the
        *  pulse a gaussian shape is assumed. All lengths and the wavelength is given in micro meters. The default wavelength is set to 1.0&mu;m
        *  and the default pulse width 100fs. As spectral width, the full width at half maximum (FWHM) is used. The result is stored in a SuperArray trafo.SAres,
        *  which holds the electric field at a certain time t which was given to the class by calling the function field. For the calculation, the spectral range is
        *  divided into divisions. For each division, one raytracing is done with the central wavelength. Here, for all steps (one step is the path between the light source
        *  and the next surface or the distance between subsequent surfaces) the step length and an information about the material will be stored.
        *  With this informations, the electric field at each point of interest within this wavelength range is calculated.
        *  Within this calculation step, the wavelength range can be subdivided in a number of wavelength steps. Here it is considered, that
        *  even though the wavelength is changed from step to step, we consider that the path of the ray is not changing.
        *  At the end, all contributions for all wavelength divisions are summed up to get the overall electric field.
        *  Comment: The number of divisions and subdivision changes the repetition rate.
		*/
		class pulseCalculation
		{
			public:
				pulseCalculation(Scene S);
				/**
				* @brief Make an estimation, when the pulse hits the object the first time
				* Often it is a problem, to find the pulse in the time domain, especially for very short pulses. 
                * This function gives an estimation, by searching the first element of the array around the chosen object which was hit by
				* a ray. Then, the time the light needs to travel from the light source until this point will be calculated. 
				*/
				double findHitTime(int ObjNo); 
				void fieldCalculation(); ///< This function makes the raytracing (normally only used internally)
				void fieldCalculation(double omega); ///< This function makes one raytracing step at frequency omega
				void setPulseWidth(double dt); ///< Sets the spectral width according to the pulse width and adjusts the widht of the subdivisions				
				void setSpatialResolution(double dx); ///< sets the spatial resolution to a value near to dx
				void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList); ///< sets the list of functions, which describe the wavelength dependend refractive index (length must be number of objects + 1)
                void setNumReflex(int numReflex); ///< set the number of internal reflections to be considered
                void setSpectralRanges(int n); ///< set the number of spectral ranges in which one raytracing calculation is made
				/**
				 * @brief Set number of wavelengths per range.
				 * For the calculation, the spectral range is subdivided into a number of wavelength ranges. This method set the
				 * number of wavelength steps per spectral range.
				 */
                void setNumWavelengthsPerRange(int nS);
				void setCenterWavelength(double wvl); ///< Set center wavelength of the pulse
				/**
				* @brief Set Bandwith of the light source(s)
				* Approximately, \f[\frac{\Delta\lambda}{\lamda}=\frac{\Delta\omega}{\omega}\Rightarrow \Delta\omega=\Delta\lambda \cdot \frac{2\pi c}{\lambda^2}\f]
				*/
				void setBandwidth(double dWvl); ///< Set Bandwith of the light source(s)
				/**
				* @brief Set the repetition rate (in fs)
				* The repetition rate is the frequency resolution within the calculation. This function sets trafoparms::nS to
				* the next integer value. The used repetition rate used in the calculation is then given by bandwith df/(nI*nS)  
				*/
				void setRepetitionRate(double rep);
                void setNumberOfThreads(int n); ///< set the number of threads used for the calculation
                int getNumberOfThreads(); ///< get the number of threads used for the calculation

				double field(double t, int settings=PULSECALCULATION_CLEAR_SA); ///< This function calculates the fields at time. Keep in mind, that it works only if the class has the list with the refractive index functions
				void reset(); ///< Clears all arrays 		
				void setReferenceTime(double tref);
				Trafo trafo;
				SuperArray<GOAT::maths::Vector<std::complex<double> > > SAres;
					std::vector<SuperArray<std::vector<gridEntry> > >  SA; ///< Here, all infos are stored to calculate the pulse (step lengths, index of the medium etc.)
					Raytrace_usp rt;

		        double domega; ///< spectral resolution
				double dWvl=0.02;  ///< spectral width of the light (default 20nm)
				TrafoParms getTrafoParms() { return trafoparms; }
				double getReferenceTime() { return tref; }
                /*
                 * Returns the number of cells per direction, nn along each axis of the calculation space.
                 * The spatial resolution is then 2*r0/nn. This number is changed by calling #setSpatialResolution
                 */
				INDEX_TYPE getNumCellsPerDirection() { return nn; }

			private:	
				int settings=0;
				/* In this function the default values (trafoparms) for the calculations are set as follows:
				* dt  : 1E-14s
                * wvl : 1.0&mu;m
				* nI  : 1
				* 
				*/
				void setDefaults();				
                void calcTrafoParms();
				// std::vector< std::vector<SuperArray<std::vector<gridEntry> > > > SA;

				double dRWvl;      ///< spectral width of one subdivision
                INDEX_TYPE  nn;    ///< number of cells over the whole width of the calculation space (i.e. 2*r0).
				Scene S;
				
				bool raytracingDone = false; ///< If true, the raytracing part was done and the field calculation starts directly				

				TrafoParms trafoparms;
				double tref = 0.0; 
				int numReflex = INEL_MAX_NREFLEX;		
				int fieldCalls = 0; ///< how many times was field called (after last reset)
                int number_of_threads=5; ///< number of threads used for calculation
		};
	}
}
