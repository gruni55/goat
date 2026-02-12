#include "pulsecalculation_kirchhoff.h"
#include "kirchhoff.h"

namespace GOAT
{
  namespace raytracing
  {
	pulseCalculationKirchhoff::pulseCalculationKirchhoff(Scene S) : pulseCalculationBase(S)
	{
	}

	void pulseCalculationKirchhoff::setSpatialResolution(double dx)
	{

		pulseCalculationBase::setSpatialResolution(dx);
	    S.setNumberOfCellsPerDirection(nn);
		
	}

	void pulseCalculationKirchhoff::initCalculation(double& omegaStart, double& domega)  
	{
		pulseCalculationBase::initCalculation(omegaStart, domega);
		rt = Raytrace_pure(S);
		rt.setNumReflex(numReflex);
		sigma = trafoparms.dt / (2.0 * M_LN2);
	}

	void pulseCalculationKirchhoff::oneFrequency(double t, double omega, double omega0)
	{
		/* set weight for the current frequency, which is determined by the Fourier transform of the pulse.
		The Fourier transform of a Gaussian pulse is also a Gaussian function, which is given by exp(-dw ^ 2 * sigma ^ 2 / 2),
		where dw is the difference between the current frequency and the center frequency, and sigma is the width of the pulse in the time domain.
		The width of the pulse in the time domain is related to the width of the pulse in the frequency domain by sigma = dt / (2 * sqrt(2 * ln(2))),
		where dt is the pulse width in femto seconds. */

		double dw = omega - omega0;
		double wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega;		
		std::complex<decltype(dw)> weight = exp(-dw * dw * sigma * sigma / 2.0);

		// prepare everything for the raytracing part, i.e. set the wavelength and the refractive index for each object and the surrounding medium
		rt.S.cleanAllDetectors();

		rt.S.setRaytype(LIGHTSRC_RAYTYPE_IRAY);
		for (int i = 0; i < S.nObj; i++)
			rt.S.Obj[i]->setn(trafoparms.nList[i](wvl));

		for (int i = 0; i < S.nLS; i++)
			rt.S.LS[i]->setWavelength(wvl);
		std::cout << "Raytracing for frequency " << omega << " (wavelength " << wvl << " um) with weight " << weight << std::endl;
		rt.trace();
		//rt.S.Det[0]->save("C:\\tmp\\test.dat");
		rt.S.multAllDetectors(weight);
		std::cout << "Raytracing finished for frequency " << omega << std::endl;
		for (int i = 0; i < rt.S.nK3D; i++)
		{
			rt.S.k3D[i]->calc(wvl);
		}
		

	}
  }
}