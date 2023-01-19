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
		std::complex<double> n_BK7(double wvl)
		{
			double wvl2 = wvl * wvl;
			std::complex<double> n = sqrt(1.03961212 * wvl2 / (wvl2 - 0.00600069867) + 0.231792344 * wvl2 / (wvl2 - 0.0200179144) + 1.01046945 * wvl2 / (wvl2 - 103.560653) + 1);
			return n;
		}

		/**
		* @brief refractive index of air according to refractiveindex.info
		* @param wvl wavelength (in um)
		* @return refractive index
		*/
		std::complex<double> n_Air(double wvl)
		{
			double wvl_2 = 1.0 / (wvl * wvl);
			std::complex<double> n = 0.05792105 / (238.0185 - wvl_2) + 0.00167917 / (57.362 - wvl_2) + 1.0;
			return n;
		}

		/**
		* @brief refractive index of vacuum (=1.0) for compatibility 
		*/
		std::complex<double> n_Vacuum(double wvl)
		{
			return 1.0;
		}
	}
}