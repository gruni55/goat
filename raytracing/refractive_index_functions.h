#pragma once
#include <complex>
namespace GOAT
{
	namespace raytracing
	{
		/**
		* @brief refractive index of BK7 (glass, Schott) according to refractiveindex.info
		* @param wvl wavelength
		* @return refractive index
		*/
		std::complex<double> n_BK7(double wvl);

		/**
		* @brief refractive index of 1.5 
		* @param wvl wavelength (in um)
		* @return refractive index
		*/
		std::complex<double> n_Glass(double wvl);

		/**
		* @brief refractive index of fused silica according to refractiveindex.info
		* @param wvl wavelength (in um)
		* @return refractive index
		*/
		std::complex<double> n_fused_silica(double wvl);

		/**
		* @brief refractive index of air according to refractiveindex.info
		* @param wvl wavelength (in um)
		* @return refractive index
		*/
		std::complex<double> n_Air(double wvl);

		/**
		* @brief refractive index of vacuum (=1.0) for compatibility 
		*/
		std::complex<double> n_Vacuum(double wvl);

		/**
* @brief refractive index formula used for glasses from Schott
*/
		std::complex<double> n_Sellmeier(double lambda, double C[3], double B[3]);

		/**
		* @brief refractive index of LASF55-glass from Schott (coefficients according taken from www.schott.com)
		*/
		std::complex<double> n_LASF55(double wvl);
		                
                std::complex<double> n_lin(double wvl);

        /**
		 * @brief refractive index function of Poly(methyl mathacrylate), PMMA 
		 * refractive index function for PMMA taken from refractiveindex.info (Sultanova et al., Acta Physica Polonica A 116, 585-587 (2009))  
		*/
		std::complex<double> n_PMMA (double wvl);

		/** 
		* @brief refractive index function for highly absorbing medium with real part 1.0
		*/

		std::complex<double> n_ABS(double wvl);

		std::complex<double> n_test (double wvl);
	}
}
