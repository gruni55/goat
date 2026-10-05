#pragma once

#include "superarray.h"
#include "fft.h"
#include "raytrace_usp.h"
#include <vector>
#include "computeSettings.h"
namespace GOAT
{
  namespace raytracing
  {
	  /**
	  * @brief This class provides the base functionality for the pulse calculation. 
	  */
	  class pulseCalculationBase
	  {
	  public:
		  pulseCalculationBase(Scene S); ///< constructor with Scene Object
		  /**
		  * This method sets all parameters to default values. The default values are as follows:
		  * | Parameter | Default value |
		  * |-----------|---------------|
		  * | Pulse width (trafoparms.dt) | 100 fs |
		  * | Center wavelength (trafoparms.wvl) | 1 µm |
		  * | number of spectral ranges (trafoparms.nI) | 250 |
		  * | number of subdivision per spectral range (trafoparms.nS) | 1 |
		  * | spatial resolution  | 1 µm |
		  */
		  void setComputeSettings(GOAT::computeSettings settings) { this->settings = settings; }
		  void setDefaults(); 
		  void setPulseWidth(double dt); ///< sets the pulse width (in femto seconds)
		  void setSpatialResolution(double dx); ///< sets the spatial resolution of the calculation grid (in µm)
		  void setNumReflex(int numReflex); ///< sets the number of reflections per ray considered in the raytracing part
		  void setCenterWavelength(double wvl); ///< sets the center wavelength of the pulse (in µm)
		  void setBandwidth(double dWvl); ///< sets the bandwidth of the pulse (in nm)
		  void setRepetitionRate(double rep); ///< sets the repetition rate of the pulse train (in femto seconds). If this value is >0, the pulse train will be considered in the calculation
		  void setSpectralRanges(int nI); ///< sets the number of spectral ranges. The spectral range is the range of frequencies, which are considered in the calculation. The spectral width is determined by the pulse width and the center wavelength. The spectral range is subdivided into nI spectral points, which are calculated separately.
		  void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double)>> nList); ///< sets the refractive index functions. The refractive index functions describe the wavelength dependence of the refractive index for each object and the surrounding medium. The order of the functions has to be as follows: first, the function for the surrounding medium, then one function for each object in the same order as they are stored in the scene.
		  void setNumWavelengthsPerRange(int nS); ///< sets the number of subdivision per spectral range. The spectral range is subdivided into nI spectral points, which are calculated separately. This parameter determines how many spectral points are calculated within each spectral range.
		  void field(double t); ///< calculates the electric field at time t. 
		  INDEX_TYPE getNumCellsPerDirection() { return nn; }  ///< returns the number of cells per direction, which is used for the calculation grid. The calculation grid is a square grid with nn x nn x nn cells, which covers the whole width of the calculation space (i.e. 2*r0).
		  Scene S;

	  protected:
		  
		  virtual void initCalculation(double& omegaStart, double& domega);
		  virtual void oneFrequency(double t, double omega, double omega0) = 0;
		  void calcTrafoParms();
		  int numReflex = INEL_MAX_NREFLEX;
		  TrafoParms trafoparms;
		  INDEX_TYPE  nn = 0;    ///< number of cells over the whole width of the calculation space (i.e. 2*r0).
		  double dWvl = 0.02;  ///< spectral width of the light (default 20nm)
		  double Domega = 0;  ///< spectral width in frequencies (unit: fs^-1)
		  
		  computeSettings settings;
	  };
   }
}