#include "pulsecalculation_rt.h"
#include <chrono>
namespace GOAT
{
	namespace raytracing
	{
		pulseCalculation_rt::pulseCalculation_rt()
		{
			setDefaults();
		}

		pulseCalculation_rt::pulseCalculation_rt(Scene S)
		{
			this->S = S;
			setDefaults();
		}

		void pulseCalculation_rt::setDefaults()
		{
			trafoparms.dt = 100;
			trafoparms.wvl = 1.0;
			trafoparms.nI = 10;
			trafoparms.nR = 1;
			trafoparms.nS = 50;
			setPulseWidth(trafoparms.dt);
			setSpatialResolution(1.0);
		}
	

			void pulseCalculation_rt::calcTrafoParms()
			{
				double Sigma = sqrt(2.0 * M_LN2) / trafoparms.dt;
				Domega = 8.0 * Sigma;
				trafoparms.omegaStart = trafoparms.omega0 - Domega / 2.0;
				trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
				double lambdaStart = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.omegaEnd;
				double lambdaEnd = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.omegaStart;
				std::cout << "% wvl-range:" << lambdaStart << "\t" << lambdaEnd << std::endl;
			}

			void pulseCalculation_rt::field(double t)
			{
				double omega0 = 2.0 * M_PI * C_LIGHT_MU_FS / trafoparms.wvl;
				Domega = 5.0 * 4.0 * M_LN2 / trafoparms.dt;
				std::cout << "Domega=" << Domega << std::endl;
				double domega = Domega / (double)trafoparms.nI;
				double omegaStart = omega0 - Domega / 2.0;
				double omega;
				double wvl1, wvl2;
				rt = Raytrace_usp_rt(S, nn);
				rt.setRefractiveIndexFunctions(trafoparms.nList);
				rt.setNumReflex(numReflex);
				double wvl;
				double sigma = trafoparms.dt / (2.0 * M_LN2);
				std::complex<double> weight;
				double dw;

				// loop over the frequency ranges
				for (int iOmega = 0; iOmega < trafoparms.nI; iOmega++)
				{
					omega = omegaStart + (double)iOmega * domega;
					dw = omega - omega0;
					weight = exp(-I * omega * t) * exp(-dw * dw * sigma * sigma / 2.0);
					wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega; // center wavelength of the current range

					// ------ for output only ------
				//	wvl1 = 2.0 * M_PI * C_LIGHT_MU_FS / (omega - 0.5 * domega);
				//	wvl2 = 2.0 * M_PI * C_LIGHT_MU_FS / (omega + 0.5 * domega);
					std::cout << "%  " << iOmega  << std::endl << std::flush;
				auto starttime=std::chrono::high_resolution_clock::now();
					rt.trace(omega, weight);
				auto rttime=std::chrono::high_resolution_clock::now();
				std::cout << "% time for rt: " << std::chrono::duration_cast<std::chrono::microseconds>(rttime - starttime).count() / 1000000.0 << " s" << std::endl;
				}
			}

			void pulseCalculation_rt::setPulseWidth(double dt)
			{
				trafoparms.dt = dt;
				calcTrafoParms();
			}

			void pulseCalculation_rt::setCenterWavelength(double wvl)
			{
				trafoparms.wvl = wvl;
				trafoparms.omega0 = C_LIGHT_MU_FS / wvl * 2.0 * M_PI;
				calcTrafoParms();
			}

			void pulseCalculation_rt::setBandwidth(double dWvl)
			{
				this->dWvl = dWvl;
				 Domega = 2.0 * M_PI * C_LIGHT_MU_FS * dWvl / (trafoparms.wvl * trafoparms.wvl);
				trafoparms.omegaEnd = trafoparms.omega0 + Domega / 2.0;
				trafoparms.omegaStart = trafoparms.omega0 + Domega / 2.0;
			}

			void pulseCalculation_rt::setRepetitionRate(double rep)
			{
				 Domega = trafoparms.omegaEnd - trafoparms.omegaStart;
				trafoparms.nS = ceil(Domega / (rep * (double)trafoparms.nI));
			}
			

			void pulseCalculation_rt::setSpectralRanges(int nI)
			{
				trafoparms.nI = nI;
			}
			


			void pulseCalculation_rt::setSpatialResolution(double dx)
			{
				nn = 2.0 * S.r0 / dx;
			}

			void pulseCalculation_rt::setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList)
			{
				trafoparms.nList = nList;
				rt.setRefractiveIndexFunctions(nList);
			}

			void pulseCalculation_rt::setNumReflex(int numReflex)
			{
				this->numReflex = numReflex;
				rt.setNumReflex(numReflex);
			}				
	}
}
