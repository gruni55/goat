//#pragma once

#include "superarray.h"
#include "fft.h"
#include "raytrace_usp.h"
#include <vector>
#include "superarray.h"
#include "raytrace_field_usp.h"



namespace GOAT
{
	namespace raytracing
	{
	    /*! \brief This class provides functionality to calculate field distributions for short pulses.  
		*   
		*  The class calculates the field distribution for a short pulsed light source. Also dispersion is considered, therefore a list of 
		*  functionsis required, which describe the wavelength dependence of all objects and the surrounding medium. 
		*  Since short pulses are considered, the light has a spectral width, which depends on the pulse width (Fourier transform). For the
		*  pulse a gaussian shape is assumed. All lengths and the wavelength is given in micro meters. The default wavelength is set to 1.0�m 
		*  and the pulse width 10fs. As spectral width, the full width at half maximum (FWHM) is used. The result is stored in a SuperArray SAres, 
		*  which holds the electric field at a certain time t which was given to the class by calling the function field
		*/

template <class T>	class pulseCalculation
		{
			public:
				pulseCalculation();
				pulseCalculation(Scene S);
				/**
				* @brief Make an estimation, when the pulse hits the object the first time
				* Often it is a problem, to find the pulse in the time domain, especially for very short pulses. 
				* This function gives an estimation, by search�ng the first element of the array around the chosen object which was hit by
				* a ray. Then, the time the light needs to travel from the light source until this point will be calculated. 
				*/
				double findHitTime(int ObjNo); 
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

				void field(double t); ///< This function calculates the fields at time. Keep in mind, that it works only if the class has the list with the refractive index functions
				void reset(); ///< Clears all arrays 		
				void setReferenceTime(double tref);
				Trafo trafo;
				SuperArray<GOAT::maths::Vector<std::complex<double> > > SAres;
					std::vector<SuperArray<std::vector<gridEntry> > >  SA; ///< Here, all infos are stored to calculate the pulse (step lengths, index of the medium etc.)
					T rt;

		
				

			// protected:	
				/* In this function the default values (trafoparms) for the calculations are set as follows:
				* dt  : 1E-14s
				* wvl : 1.0�m
				* nI  : 1
				* 
				*/
				void setDefaults();				
                void calcTrafoParms();
				// std::vector< std::vector<SuperArray<std::vector<gridEntry> > > > SA;
				double domega; ///< spectral resolution
				double dWvl=0.02;  ///< spectral width of the light (default 20nm)
				double dRWvl;      ///< spectral width of one subdivision
				INDEX_TYPE  nn;    // number of cells over the whole width of the calculation space (i.e. 2*r0)
				Scene S;
				
				bool raytracingDone = false; ///< If true, the raytracing part was done and the field calculation starts directly				

				TrafoParms trafoparms;
				double tref = 0.0; 
				int numReflex = INEL_MAX_NREFLEX;					
		};

// ------- some specializations ------------
				class pulseCalculation_field : public pulseCalculation<Raytrace_field_usp>
				{
				public: 					
					pulseCalculation_field(Scene& S) : pulseCalculation<Raytrace_field_usp>(S) {}
					void addBoxDetector(Box* box) { rt.addBoxDetector(box); }; ///< add a box as detector
					void field(double t); 					
				};


				void pulseCalculation_field::field(double t)
				{
					double omega0 = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.wvl;
					double Domega = 8.0 * 4.0 * M_LN2 / trafoparms.dt;
					std::cout << "Domega=" << Domega << std::endl;


					//			double Domega = 8 * M_PI * C_LIGHT_MU_FS * dWvl / (4.0 * trafoparms.wvl * trafoparms.wvl - dWvl * dWvl);
					//                        double Domega = 2.0 * M_PI * C_LIGHT_MU_FS dWvl / (trafoparms.wvl * trafoparms.wvl);
					double domega = Domega / (double)trafoparms.nI;
					double omegaStart = omega0 - Domega / 2.0;
					double omega;
					double wvl1, wvl2;
					rt = Raytrace_field_usp(S);
					double wvl;
					trafo.initResult(S.r0, rt.SE.nges[0], rt.SE.nges[1], rt.SE.nges[2], S.Obj, S.nObj);
					// loop over the frequency ranges
					for (int iOmega = 0; iOmega < trafoparms.nI; iOmega++)
					{
						omega = omegaStart + (double)iOmega * domega;
						wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega; // center wavelength of the current range

						// ------ for output only ------
						wvl1 = 2.0 * M_PI * C_LIGHT_MU_FS / (omega - 0.5 * domega);
						wvl2 = 2.0 * M_PI * C_LIGHT_MU_FS / (omega + 0.5 * domega);
						std::cout << "%  " << iOmega << ":start FFT (" << wvl << "�m)" << "\twvl1=" << wvl1 << "\twvl2=" << wvl2 << "\tomega=" << omega << std::endl << std::flush;

						fieldCalculation(omega); // do the raytracing							
						trafo.calc(static_cast<std::vector<SuperArray<std::vector<gridEntry> > > >(rt.SA), omega - domega * 0.5, omega + domega * 0.5, t); // do the Fourier transform
					}
				}
				
	}
}

#include "pulsecalculation.hpp"