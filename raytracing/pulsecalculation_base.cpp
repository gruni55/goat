#include <chrono>
#include "pulsecalculation_base.h"

namespace GOAT
{
  namespace raytracing
  {
	pulseCalculationBase::pulseCalculationBase(Scene S)
	{
	  this->S = S;
	  setDefaults();
	}
	void pulseCalculationBase::setDefaults()
	{
	  trafoparms.dt = 100;
	  trafoparms.wvl = 1.0;
	  trafoparms.nI = 10;
	  trafoparms.nR = 1;
	  trafoparms.nS = 50;
	  setPulseWidth(trafoparms.dt);
	  setSpatialResolution(1.0);
	  setCenterWavelength(trafoparms.wvl);

	}
	void pulseCalculationBase::setPulseWidth(double dt)
	{
	  this->trafoparms.dt = dt;
	  calcTrafoParms();
	}
	void pulseCalculationBase::setSpatialResolution(double dx)
	{
		nn = 2.0 * S.r0 / dx;
	}

	void pulseCalculationBase::field(double t)
	{
		double omegaStart, domega;
		initCalculation(omegaStart, domega);
		for (int iOmega = 0; iOmega < trafoparms.nI; iOmega++)
		{
			auto start = std::chrono::high_resolution_clock::now();
			double omega = omegaStart + (double)iOmega * domega;
			oneFrequency(t, omega, trafoparms.omega0);
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << "% time for frequency " << iOmega << ": " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " s" << std::endl;
		}
	}

	void pulseCalculationBase::initCalculation(double &omegaStart, double &domega)
	{
		double omega0 = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.wvl;
		Domega = 5.0 * 4.0 * M_LN2 / trafoparms.dt;
		domega = Domega / (double)trafoparms.nS;
		omegaStart = omega0 - Domega / 2.0;
	}

	void pulseCalculationBase::calcTrafoParms()
	{
		double Sigma = sqrt(2.0 * M_LN2) / trafoparms.dt;
		double Domega = 8.0 * Sigma;
		trafoparms.omegaStart = trafoparms.omega0 - Domega / 2.0;
		trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
		double lambdaStart = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.omegaEnd;
		double lambdaEnd = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.omegaStart;
	}

	void pulseCalculationBase::setNumReflex(int numReflex)
	{
		this->numReflex = numReflex;
	}

	void pulseCalculationBase::setCenterWavelength(double wvl)
	{
		trafoparms.wvl = wvl;
		trafoparms.omega0 = C_LIGHT_MU_FS / wvl * 2.0 * M_PI;
		calcTrafoParms();
	}

	void pulseCalculationBase::setBandwidth(double dWvl)
	{
		this->dWvl = dWvl;
		Domega = 2.0 * M_PI * C_LIGHT_MU_FS * dWvl / (trafoparms.wvl * trafoparms.wvl);
		trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
		trafoparms.omegaStart = trafoparms.omega0 + Domega / 2.0;
	}

	void pulseCalculationBase::setRepetitionRate(double rep)
	{
		double Domega = trafoparms.omegaEnd - trafoparms.omegaStart;
		trafoparms.nS = ceil(Domega / (rep * (double)trafoparms.nI));
	}

	void pulseCalculationBase::setSpectralRanges(int nI)
	{
		trafoparms.nI = nI;		
	}

	void pulseCalculationBase::setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double)>> nList)
	{
		trafoparms.nList = nList;	
	}

	void pulseCalculationBase::setNumWavelengthsPerRange(int nS)
	{
		trafoparms.nS = nS;
	}
  }
}