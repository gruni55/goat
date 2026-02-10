#pragma once

#include "superarray.h"
#include "fft.h"
#include "raytrace_usp.h"
#include <vector>
namespace GOAT
{
  namespace raytracing
  {
	  class pulseCalculationBase
	  {
	  public:
		  pulseCalculationBase(Scene S);
		  void setDefaults();
		  void setPulseWidth(double dt);
		  void setSpatialResolution(double dx);
		  void setNumReflex(int numReflex);
		  void setCenterWavelength(double wvl);
		  void setBandwidth(double dWvl);
		  void setRepetitionRate(double rep);
		  void setSpectralRanges(int nI);
		  void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double)>> nList);
		  void setNumWavelengthsPerRange(int nS);
		  void field(double t);
		  INDEX_TYPE getNumCellsPerDirection() { return nn; }

	  protected:
		  
		  virtual void initCalculation(double& omegaStart, double& domega);
		  virtual void oneFrequency(double t, double omega, double omega0) = 0;
		  void calcTrafoParms();
		  int numReflex = INEL_MAX_NREFLEX;
		  TrafoParms trafoparms;
		  INDEX_TYPE  nn = 0;    ///< number of cells over the whole width of the calculation space (i.e. 2*r0).
		  double dWvl = 0.02;  ///< spectral width of the light (default 20nm)
		  double Domega = 0;  ///< spectral width in frequencies (unit: fs^-1)
		  Scene S;
	  };
   }
}