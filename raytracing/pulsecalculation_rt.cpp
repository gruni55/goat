#include "pulsecalculation_rt.h"
#include <chrono>
namespace GOAT
{
	namespace raytracing
	{
		pulseCalculation_rt::pulseCalculation_rt() : pulseCalculationBase(Scene())
		{
		}

		pulseCalculation_rt::pulseCalculation_rt(Scene S) : pulseCalculationBase(S)
		{
		}

		void pulseCalculation_rt::initCalculation(double& omegaStart, double& domega)
		{
			pulseCalculationBase::initCalculation(omegaStart, domega);
			 rt = Raytrace_usp_rt(S, nn);
			 rt.setRefractiveIndexFunctions(trafoparms.nList);
			 rt.setNumReflex(numReflex);
			 sigma = trafoparms.dt / (2.0 * M_LN2);
		}

		void pulseCalculation_rt::oneFrequency(double t, double omega, double omega0)
		{
			double dw = omega - omega0;
			double wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega;
			std::complex<decltype(dw)> weight = exp(-dw * dw * sigma * sigma / 2.0);
			rt.trace(omega, weight);
		}

		

	}
}
