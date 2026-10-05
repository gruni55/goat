#pragma once

#include "pulsecalculation_rt.h"
#include "kirchhoff.h"
namepace GOAT
{
	namespace raytracing
	{
		class pulseCalculation_Kirchhoff : public pulseCalculation_rt
		{
			public:
				pulseCalculation_Kirchhoff();
				pulseCalculation_Kirchhoff(Scene S);
				void oneFrequency(double omega, std::complex<double> weight) override;
		};
	}
}