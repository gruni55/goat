#pragma once
#include "superarray.h"
#include "raytrace_usp_rt.h"
#include <vector>
#include "fft.h"
#include "pulsecalculation_base.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		 * @brief This class provides pulse calculation (only with ray tracing)
		 * The pulses from a mode-locked laser can be described by series of frequency modes. 
		 * The electric field at a certain location at time t can be described by
		 * \f$ \vec E(\vec r,t)=\sum\f$
		 */
		class pulseCalculation_rt : public pulseCalculationBase
		{
		public:
			pulseCalculation_rt();
			pulseCalculation_rt(Scene S);
			Raytrace_usp_rt rt;

		protected:		
			double sigma; 
			void initCalculation(double& omegaStart, double& domega) override;
			void oneFrequency(double t, double omega, double omega0);
		};
	}
}