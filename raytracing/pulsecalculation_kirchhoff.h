#pragma once
#include "pulsecalculation_base.h"
#include "raytrace.h"

namespace GOAT
{
  namespace raytracing
  {
	class pulseCalculationKirchhoff : public pulseCalculationBase
	{
	public:
		pulseCalculationKirchhoff (Scene S);
		void setSpatialResolution(double dx); ///< sets the spatial resolution of the calculation grid (in µm)
	protected:
	  void initCalculation(double& omegaStart, double& domega) override;
	  void oneFrequency(double t, double omega, double omega0) override;
	  Raytrace_pure rt;
	  double sigma;
	};
  }
}