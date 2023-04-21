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

		/**
* @brief refractive index formula used for glasses from Schott
*/
		std::complex<double> n_Sellmeier(double lambda, double C[3], double B[3])
		{
			std::complex<double> n;
			double l2 = lambda * lambda;
			n = sqrt(B[1] * l2 / (l2 - C[1]) + B[2] * l2 / (l2 - C[2]) + B[3] * l2 / (l2 - C[3]) + 1.0);
			return n;
		}

		/**
		* @brief refractive index of LASF55-glass from Schott (coefficients according taken from www.schott.com)
		*/
		std::complex<double> n_LASF55(double wvl)
		{
			double B[3] = { 2.30861228,0.354736638,1.922271250 };
			double C[3] = { 0.01304469950,0.0557524221,133.1968690 };
			return n_Sellmeier(wvl, C, B);
		}
                
                std::complex<double> n_lin(double wvl)
		{
			double n0 = 1.5;
			double m = -0.1;
			double c = n0 - m;
			return m * wvl + c;
		}

		std::complex<double> n_test (double wvl)
		{
			if (fabs(wvl-1.0)<0.001) return 1.5;
			return 1.4;      
		}

	}
}
